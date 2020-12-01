import torch
import numpy as np
from . import quadrics
from external.ldif.util import camera_util
import torch.nn.functional as F
from external.ldif.util import file_util


def homogenize(m):
  # Adds homogeneous coordinates to a [..., N,N] matrix.
  m = F.pad(m, [0, 1, 0, 1], "constant", 0)
  m[..., -1, -1] = 1
  return m

def _unflatten(config, vector):
    return torch.split(vector, [1, 3, 6, config['model']['mesh_reconstruction']['implicit_parameter_length']], -1)

class StructuredImplicit(object):
    def __init__(self, config, constant, center, radius, iparam, net=None):
        # (ldif.representation.structured_implicit_function.StructuredImplicit.from_activation)
        self.config = config
        self.implicit_parameter_length = config['model']['mesh_reconstruction']['implicit_parameter_length']
        self.element_count = config['model']['mesh_reconstruction']['element_count']
        self.sym_element_count = config['model']['mesh_reconstruction']['sym_element_count']

        self.constants = constant
        self.radii = radius
        self.centers = center
        self.iparams = iparam
        self.effective_element_count = self.element_count + self.sym_element_count
        self.device = constant.device
        self.batch_size = constant.size(0)
        self.net = net
        self._packed_vector = None
        self._analytic_code = None
        self._all_centers = None

    @classmethod
    def from_packed_vector(cls, config, packed_vector, net):
        """Parse an already packed vector (NOT a network activation)."""
        constant, center, radius, iparam = _unflatten(config, packed_vector)
        return cls(config, constant, center, radius, iparam, net)

    @classmethod
    def from_activation(cls, config, activation, net):
        constant, center, radius, iparam = _unflatten(config, activation)
        constant = -torch.abs(constant)
        radius_var = torch.sigmoid(radius[..., :3])
        radius_var = 0.15 * radius_var
        radius_var = radius_var * radius_var
        max_euler_angle = np.pi / 4.0
        radius_rot = torch.clamp(radius[..., 3:], -max_euler_angle, max_euler_angle)
        radius = torch.cat([radius_var, radius_rot], -1)
        center /= 2
        return cls(config, constant, center, radius, iparam, net)

    def _tile_for_symgroups(self, elements):
        # Tiles an input tensor along its element dimension based on symmetry
        # (ldif.representation.structured_implicit_function._tile_for_symgroups)
        sym_elements = elements[:, :self.sym_element_count, ...]
        elements = torch.cat([elements, sym_elements], 1)
        return elements

    def _generate_symgroup_samples(self, samples):
        samples = samples.unsqueeze(1).expand(-1, self.element_count, -1, -1)
        sym_samples = samples[:, :self.sym_element_count].clone()
        # sym_samples *= torch.tensor([1, 1, -1], dtype=torch.float32, device=self.device)  # reflect across the XY plane
        sym_samples *= torch.tensor([-1, 1, 1], dtype=torch.float32, device=self.device)  # reflect across the YZ plane
        effective_samples = torch.cat([samples, sym_samples], 1)
        return effective_samples

    def compute_world2local(self):
        tx = torch.eye(3, device=self.device).expand(self.batch_size, self.element_count, -1, -1)
        centers = self.centers.unsqueeze(-1)
        tx = torch.cat([tx, -centers], -1)
        lower_row = torch.tensor([0., 0., 0., 1.], device=self.device).expand(self.batch_size, self.element_count, 1, -1)
        tx = torch.cat([tx, lower_row], -2)

        # Compute a rotation transformation
        rotation = camera_util.roll_pitch_yaw_to_rotation_matrices(self.radii[..., 3:6]).inverse()
        diag = 1.0 / (torch.sqrt(self.radii[..., :3] + 1e-8) + 1e-8)
        scale = torch.diag_embed(diag)

        # Apply both transformations and return the transformed points.
        tx3x3 = torch.matmul(scale, rotation)
        return torch.matmul(homogenize(tx3x3), tx)

    def implicit_values(self, local_samples):
        # Computes the implicit values given local input locations.
        iparams = self._tile_for_symgroups(self.iparams)
        values = self.net.eval_implicit_parameters(iparams, local_samples)
        return values

    @property
    def all_centers(self):
        if self._all_centers is None:
            sym_centers = self.centers[:, :self.sym_element_count].clone()
            sym_centers[:, :, 0] *= -1  # reflect across the YZ plane
            self._all_centers = torch.cat([self.centers, sym_centers], 1)
        return self._all_centers

    def class_at_samples(self, samples, apply_class_transfer=True):
        # (ldif.representation.structured_implicit_function.StructuredImplicit.class_at_samples)
        effective_constants = self._tile_for_symgroups(self.constants)
        effective_centers = self._tile_for_symgroups(self.centers)
        effective_radii = self._tile_for_symgroups(self.radii)

        effective_samples = self._generate_symgroup_samples(samples)
        constants_quadrics = torch.zeros(self.batch_size, self.effective_element_count, 4, 4, device=self.device)
        constants_quadrics[:, :, -1:, -1] = effective_constants

        per_element_constants, per_element_weights = quadrics.compute_shape_element_influences(
            constants_quadrics, effective_centers, effective_radii, effective_samples
        )

        # We currently have constants, weights with shape:
        # [batch_size, element_count, sample_count, 1].
        # We need to use the net to get a same-size grid of offsets.
        # The input samples to the net must have shape
        # [batch_size, element_count, sample_count, 3], while the current samples
        # have shape [batch_size, sample_count, 3]. This is because each sample
        # should be evaluated in the relative coordinate system of the
        # The world2local transformations for each element. Shape [B, EC, 4, 4].
        effective_world2local = self._tile_for_symgroups(self.compute_world2local())
        local_samples = torch.matmul(F.pad(effective_samples, [0, 1], "constant", 1),
                                     effective_world2local.transpose(-1, -2))[..., :3]
        implicit_values = self.implicit_values(local_samples)

        residuals = 1 + implicit_values
        local_decisions = per_element_constants * per_element_weights * residuals
        local_weights = per_element_weights
        sdf = torch.sum(local_decisions, 1)
        if apply_class_transfer:
            sdf = torch.sigmoid(100 * (sdf + 0.07))

        return sdf, (local_decisions, local_weights)

    @property
    def vector(self):
        if self._packed_vector is None:
            self._packed_vector = torch.cat([self.constants, self.centers, self.radii, self.iparams], -1)
        return self._packed_vector

    @property
    def analytic_code(self):
        if self._analytic_code is None:
            self._analytic_code = torch.cat([self.constants, self.centers, self.radii], -1)
        return self._analytic_code

    def savetxt(self, path):
        assert self.vector.shape[0] == 1
        sif_vector = self.vector.squeeze().cpu().numpy()
        sif_vector[:, 4:7] = np.sqrt(np.maximum(sif_vector[:, 4:7], 0))
        out = 'SIF\n%i %i %i\n' % (self.element_count, 0, self.implicit_parameter_length)
        for row_idx in range(self.element_count):
            row = ' '.join(10 * ['%.9g']) % tuple(sif_vector[row_idx, :10].tolist())
            symmetry = int(row_idx < self.sym_element_count)
            row += ' %i' % symmetry
            implicit_params = ' '.join(self.implicit_parameter_length * ['%.9g']) % (
                tuple(sif_vector[row_idx, 10:].tolist()))
            row += ' ' + implicit_params
            row += '\n'
            out += row
        file_util.writetxt(path, out)

    def unbind(self):
        return [StructuredImplicit.from_packed_vector(self.config, self.vector[i:i+1], self.net)
                for i in range(self.vector.size(0))]

    def __getitem__(self, item):
        return StructuredImplicit.from_packed_vector(self.config, self.vector[item], self.net)

    def dict(self):
        return {'constant': self.constants, 'radius': self.radii, 'center': self.centers, 'iparam': self.iparams}


