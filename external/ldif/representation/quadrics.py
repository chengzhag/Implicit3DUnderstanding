import torch
import torch.nn.functional as F
from external.ldif.util import camera_util


def sample_quadric_surface(quadric, center, samples):
    # Sample the quadric surfaces and the RBFs in world space, and composite them.
    # (ldif.representation.quadrics.sample_quadric_surface)
    samples = samples - center.unsqueeze(2)
    homogeneous_sample_coords = F.pad(samples, [0, 1], "constant", 1)
    half_distance = torch.matmul(quadric, homogeneous_sample_coords.transpose(-1, -2))
    half_distance = half_distance.transpose(-1, -2)
    algebraic_distance = torch.sum(homogeneous_sample_coords * half_distance, -1, keepdim=True)
    return algebraic_distance

def decode_covariance_roll_pitch_yaw(radius, invert=False):
    # Converts 6-D radus vectors to the corresponding covariance matrices.
    # (ldif.representation.quadrics.decode_covariance_roll_pitch_yaw
    d = 1.0 / (radius[..., :3] + 1e-8) if invert else radius[..., :3]
    diag = torch.diag_embed(d)
    rotation = camera_util.roll_pitch_yaw_to_rotation_matrices(radius[..., 3:6])
    return torch.matmul(torch.matmul(rotation, diag), rotation.transpose(-1, -2))

def sample_cov_bf(center, radius, samples):
    # Samples gaussian radial basis functions at specified coordinates.
    # (ldif.representation.quadrics.sample_cov_bf)
    diff = samples - center.unsqueeze(2)
    x, y, z = diff.unbind(-1)

    inv_cov = decode_covariance_roll_pitch_yaw(radius, invert=True)
    inv_cov = torch.reshape(inv_cov, [inv_cov.shape[0], -1, 1, 9])
    c00, c01, c02, _, c11, c12, _, _, c22 = inv_cov.unbind(-1)
    dist = (x * (c00 * x + c01 * y + c02 * z)
            + y * (c01 * x + c11 * y + c12 * z)
            + z * (c02 * x + c12 * y + c22 * z))
    dist = torch.exp(-0.5 * dist)
    return dist.unsqueeze(-1)

def compute_shape_element_influences(quadrics, centers, radii, samples):
    # compute shape element influences (ldif.representation.quadrics.compute_shape_element_influences)
    sampled_quadrics = sample_quadric_surface(quadrics, centers, samples)
    sampled_rbfs = sample_cov_bf(centers, radii, samples)
    return sampled_quadrics, sampled_rbfs