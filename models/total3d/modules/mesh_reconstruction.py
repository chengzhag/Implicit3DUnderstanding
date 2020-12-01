# Definition of PoseNet
# author: ynie
# date: March, 2020

import torch
import torch.nn as nn
from models.registers import MODULES
from configs.data_config import number_pnts_on_template, pix3d_n_classes
from models.modules import resnet
from models.modules.resnet import model_urls
import torch.utils.model_zoo as model_zoo
import torch.nn.functional as F
from net_utils.misc import weights_init, sphere_edges, sphere_faces, sphere_edge2face, sphere_adjacency, sphere_points_normals, sample_points_on_edges
from external.ldif.representation.structured_implicit_function import StructuredImplicit
import numpy as np
from external.ldif.util import np_util
from external.ldif.inference import extract_mesh
import trimesh
from external.PIFu.lib import mesh_util
from external.ldif.util import file_util
import os
import struct
import tempfile
import shutil
import subprocess


class PointGenCon(nn.Module):
    def __init__(self, bottleneck_size = 2500, output_dim = 3):
        super(PointGenCon, self).__init__()
        self.conv1 = torch.nn.Conv1d(bottleneck_size, bottleneck_size, 1)
        self.conv2 = torch.nn.Conv1d(bottleneck_size, bottleneck_size//2, 1)
        self.conv3 = torch.nn.Conv1d(bottleneck_size//2, bottleneck_size//4, 1)
        self.conv4 = torch.nn.Conv1d(bottleneck_size//4, output_dim, 1)

        self.th = nn.Tanh()
        self.bn1 = torch.nn.BatchNorm1d(bottleneck_size)
        self.bn2 = torch.nn.BatchNorm1d(bottleneck_size//2)
        self.bn3 = torch.nn.BatchNorm1d(bottleneck_size//4)

    def forward(self, x):
        x = F.relu(self.bn1(self.conv1(x)))
        x = F.relu(self.bn2(self.conv2(x)))
        x = F.relu(self.bn3(self.conv3(x)))
        x = self.th(self.conv4(x))
        return x

class EREstimate(nn.Module):
    def __init__(self, bottleneck_size=2500, output_dim = 3):
        super(EREstimate, self).__init__()
        self.conv1 = torch.nn.Conv1d(bottleneck_size, bottleneck_size, 1)
        self.conv2 = torch.nn.Conv1d(bottleneck_size, bottleneck_size//2, 1)
        self.conv3 = torch.nn.Conv1d(bottleneck_size//2, bottleneck_size//4, 1)
        self.conv4 = torch.nn.Conv1d(bottleneck_size//4, output_dim, 1)

        self.bn1 = torch.nn.BatchNorm1d(bottleneck_size)
        self.bn2 = torch.nn.BatchNorm1d(bottleneck_size//2)
        self.bn3 = torch.nn.BatchNorm1d(bottleneck_size//4)

    def forward(self, x):
        x = F.relu(self.bn1(self.conv1(x)))
        x = F.relu(self.bn2(self.conv2(x)))
        x = F.relu(self.bn3(self.conv3(x)))
        x = self.conv4(x)
        return x


@MODULES.register_module
class DensTMNet(nn.Module):
    def __init__(self, cfg, optim_spec=None, bottleneck_size=1024, n_classes=pix3d_n_classes, pretrained_encoder=True):
        super(DensTMNet, self).__init__()

        '''Optimizer parameters used in training'''
        self.optim_spec = optim_spec

        '''Module parameters'''
        self.num_points = number_pnts_on_template
        self.subnetworks = cfg.config['data']['tmn_subnetworks']
        self.train_e_e = cfg.config['data']['with_edge_classifier']

        '''Modules'''
        self.encoder = resnet.resnet18_full(pretrained=False, num_classes=1024,
                                            input_channels=4 if cfg.config['data'].get('mask', False) else 3)
        self.decoders = nn.ModuleList(
            [PointGenCon(bottleneck_size=3 + bottleneck_size + n_classes) for i in range(0, self.subnetworks)])

        if self.train_e_e:
            self.error_estimators = nn.ModuleList(
                [EREstimate(bottleneck_size=3 + bottleneck_size + n_classes, output_dim=1) for i in range(0, max(self.subnetworks-1, 1))])

            self.face_samples = cfg.config['data']['face_samples']

        # initialize weight
        self.apply(weights_init)

        # initialize resnet
        if pretrained_encoder:
            pretrained_dict = model_zoo.load_url(model_urls['resnet18'])
            model_dict = self.encoder.state_dict()
            if pretrained_dict['conv1.weight'].shape != model_dict['conv1.weight'].shape:
                model_dict['conv1.weight'][:,:3,...] = pretrained_dict['conv1.weight']
                pretrained_dict.pop('conv1.weight')
            pretrained_dict = {k: v for k, v in pretrained_dict.items() if
                               k in model_dict and not k.startswith('fc.')}
            model_dict.update(pretrained_dict)
            self.encoder.load_state_dict(model_dict)

    def unfreeze_parts(self, loose_parts):
        # freeze all
        for param in self.parameters():
            param.requires_grad = False
        print('All layers freezed.')
        # unfreeze parts
        if 'encoder' in loose_parts:
            for param in self.encoder.parameters():
                param.requires_grad = True
            print('Encoder unfrozen.')

    def freeze_encoder(self):
        for param in self.encoder.parameters():
            param.requires_grad = False
        print('Encoder freezed.')

    def freeze_by_stage(self, stage, loose_parts):
        if stage >= 1:
            # freeze all
            for param in self.parameters():
                param.requires_grad = False
            print('All layers freezed.')

            if 'decoder' in loose_parts:
                # unfreeze the last sub-network of decoders.
                for param in self.decoders[-1].parameters():
                    param.requires_grad = True
                print('Decoder unfrozen.')

            if 'ee' in loose_parts and hasattr(self, 'error_estimators'):
                # unfreeze the last sub-network of error estimators.
                for param in self.error_estimators[-1].parameters():
                    param.requires_grad = True
                print('EE unfrozen.')

    def forward(self, image, size_cls, threshold = 0.1, factor = 1., mask_status = None, reconstruction = 'mesh'):
        mode = 'train' if self.training else 'test'
        device = image.device

        n_batch = image.size(0)
        n_edges = sphere_edges.shape[0]

        # image encoding
        image = image.contiguous()
        afeature = self.encoder(image)
        code = torch.cat([afeature, size_cls], 1)

        if mask_status is not None:
            code4recon = code[mask_status.nonzero()]
            n_batch = code4recon.size(0)
            if n_batch == 0:
                return {'mgn_afeature':afeature}
        else:
            code4recon = code

        if reconstruction is None:
            return {'mgn_afeature':afeature}

        if mode == 'test':
            current_faces = sphere_faces.clone().unsqueeze(0).to(device)
            current_faces = current_faces.repeat(n_batch, 1, 1)
        else:
            current_faces = None

        current_edges = sphere_edges.clone().unsqueeze(0).to(device)
        current_edges = current_edges.repeat(n_batch, 1, 1)

        current_shape_grid = sphere_points_normals[:, :3].t().expand(n_batch, 3, self.num_points).to(device)

        # outputs for saving
        out_shape_points = []
        out_sampled_mesh_points = []
        out_indicators = []

        # boundary faces for boundary refinement
        boundary_point_ids = torch.zeros(size=(n_batch, self.num_points), dtype=torch.uint8).to(device)
        remove_edges_list = []

        # AtlasNet deformation + topoly modification
        for i in range(self.subnetworks):
            current_image_grid = code4recon.unsqueeze(2).expand(code4recon.size(0), code4recon.size(1),
                                                           current_shape_grid.size(2)).contiguous()
            current_image_grid = torch.cat((current_shape_grid, current_image_grid), 1).contiguous()
            current_shape_grid = current_shape_grid + self.decoders[i](current_image_grid)

            # save deformed point cloud
            out_shape_points.append(current_shape_grid)

            if i == self.subnetworks - 1 and self.subnetworks > 1:
                remove_edges_list = [item for item in remove_edges_list if len(item)]
                if remove_edges_list:
                    remove_edges_list = torch.unique(torch.cat(remove_edges_list), dim=0)
                    for batch_id in range(n_batch):
                        rm_edges = remove_edges_list[remove_edges_list[:, 0] == batch_id, 1]
                        if len(rm_edges) > 0:
                            rm_candidates, counts = torch.unique(sphere_edges[rm_edges], return_counts=True)
                            boundary_ids = counts < sphere_adjacency[rm_candidates - 1].sum(1)
                            boundary_point_ids[batch_id][rm_candidates[boundary_ids] - 1] = 1

                return {'mesh_coordinates_results': out_shape_points, 'points_from_edges': out_sampled_mesh_points,
                        'point_indicators': out_indicators, 'output_edges': current_edges,
                        'boundary_point_ids': boundary_point_ids, 'faces': current_faces,
                        'mgn_afeature':afeature}

            if self.train_e_e:
                # sampling from deformed mesh
                sampled_points = sample_points_on_edges(current_shape_grid, current_edges, quantity=self.face_samples, mode=mode)

                # save sampled points from deformed mesh
                out_sampled_mesh_points.append(sampled_points)

                # preprare for face error estimation
                current_image_grid = code4recon.unsqueeze(2).expand(code4recon.size(0), code4recon.size(1), sampled_points.size(2)).contiguous()
                current_image_grid = torch.cat((sampled_points, current_image_grid), 1).contiguous()

                # estimate the distance from deformed points to gt mesh.
                indicators = self.error_estimators[i](current_image_grid)
                indicators = indicators.view(n_batch, 1, n_edges, self.face_samples)
                indicators = indicators.squeeze(1)
                indicators = torch.mean(indicators, dim=2)

                # save estimated distance values from deformed points to gt mesh.
                out_indicators.append(indicators)
                # remove faces and modify the topology
                remove_edges = torch.nonzero(torch.sigmoid(indicators) < threshold)
                remove_edges_list.append(remove_edges)

                for batch_id in range(n_batch):
                    rm_edges = remove_edges[remove_edges[:, 0] == batch_id, 1]
                    if len(rm_edges)>0:
                        # cutting edges in training
                        current_edges[batch_id][rm_edges, :] = 1
                        if mode == 'test':
                            current_faces[batch_id][sphere_edge2face[rm_edges].sum(0).type(torch.bool), :] = 1

                threshold *= factor

        return {'mesh_coordinates_results':out_shape_points, 'points_from_edges':out_sampled_mesh_points,
                'point_indicators':out_indicators, 'output_edges':current_edges,
                'boundary_point_ids':boundary_point_ids, 'faces':current_faces,
                'mgn_afeature':afeature}


class BatchedCBNLayer(nn.Module):
    def __init__(self, f_dim=32):
        super(BatchedCBNLayer, self).__init__()
        self.fc_beta = nn.Linear(f_dim, f_dim)
        self.fc_gamma = nn.Linear(f_dim, f_dim)
        self.register_buffer('running_mean', torch.zeros(1))
        self.register_buffer('running_var', torch.ones(1))

    def forward(self, shape_embedding, sample_embeddings):
        beta = self.fc_beta(shape_embedding)
        gamma = self.fc_gamma(shape_embedding)
        if self.training:
            batch_mean, batch_variance = sample_embeddings.mean().detach(), sample_embeddings.var().detach()
            self.running_mean = 0.995 * self.running_mean + 0.005 * batch_mean
            self.running_var = 0.995 * self.running_var + 0.005 * batch_variance
        sample_embeddings = (sample_embeddings - self.running_mean) / torch.sqrt(self.running_var + 1e-5)

        out = gamma.unsqueeze(1) * sample_embeddings + beta.unsqueeze(1)

        return out


class BatchedOccNetResnetLayer(nn.Module):
    def __init__(self, f_dim=32):
        super(BatchedOccNetResnetLayer, self).__init__()
        self.bn1 = BatchedCBNLayer(f_dim=f_dim)
        self.fc1 = nn.Linear(f_dim, f_dim)
        self.bn2 = BatchedCBNLayer(f_dim=f_dim)
        self.fc2 = nn.Linear(f_dim, f_dim)

    def forward(self, shape_embedding, sample_embeddings):
        sample_embeddings = self.bn1(shape_embedding, sample_embeddings)
        init_sample_embeddings = sample_embeddings

        sample_embeddings = torch.relu(sample_embeddings)
        sample_embeddings = self.fc1(sample_embeddings)
        sample_embeddings = self.bn2(shape_embedding, sample_embeddings)

        sample_embeddings = torch.relu(sample_embeddings)
        sample_embeddings = self.fc2(sample_embeddings)

        return init_sample_embeddings + sample_embeddings


class OccNetDecoder(nn.Module):
    def __init__(self, f_dim=32):
        super(OccNetDecoder, self).__init__()
        self.fc1 = nn.Linear(3, f_dim)
        self.resnet = BatchedOccNetResnetLayer(f_dim=f_dim)
        self.bn = BatchedCBNLayer(f_dim=f_dim)
        self.fc2 = nn.Linear(f_dim, 1)

    def write_occnet_file(self, path):
        """Serializes an occnet network and writes it to disk."""
        f = file_util.open_file(path, 'wb')

        def write_fc_layer(layer):
            weights = layer.weight.t().cpu().numpy()
            biases = layer.bias.cpu().numpy()
            f.write(weights.astype('f').tostring())
            f.write(biases.astype('f').tostring())

        def write_cbn_layer(layer):
            write_fc_layer(layer.fc_beta)
            write_fc_layer(layer.fc_gamma)
            running_mean = layer.running_mean.item()
            running_var = layer.running_var.item()
            f.write(struct.pack('ff', running_mean, running_var))

        # write_header
        f.write(struct.pack('ii', 1, self.fc1.out_features))
        # write_input_layer
        write_fc_layer(self.fc1)
        # write_resnet
        write_cbn_layer(self.resnet.bn1)
        write_fc_layer(self.resnet.fc1)
        write_cbn_layer(self.resnet.bn2)
        write_fc_layer(self.resnet.fc2)
        # write_cbn_layer
        write_cbn_layer(self.bn)
        # write_activation_layer
        weights = self.fc2.weight.t().cpu().numpy()
        bias = self.fc2.bias.data.item()
        f.write(weights.astype('f').tostring())
        f.write(struct.pack('f', bias))
        f.close()

    def forward(self, embedding, samples):
        sample_embeddings = self.fc1(samples)
        sample_embeddings = self.resnet(embedding, sample_embeddings)
        sample_embeddings = self.bn(embedding, sample_embeddings)
        vals = self.fc2(sample_embeddings)
        return vals


@MODULES.register_module
class LDIF(nn.Module):
    def __init__(self, cfg, optim_spec=None, n_classes=pix3d_n_classes,
                 pretrained_encoder=True):
        super(LDIF, self).__init__()

        '''Optimizer parameters used in training'''
        self.optim_spec = optim_spec

        '''Module parameters'''
        self.cfg = cfg
        self.bottleneck_size = cfg.config['model']['mesh_reconstruction'].get('bottleneck_size', 2048)
        cfg.config['model']['mesh_reconstruction']['bottleneck_size'] = self.bottleneck_size
        self.element_count = cfg.config['model']['mesh_reconstruction']['element_count']
        self.sym_element_count = cfg.config['model']['mesh_reconstruction']['sym_element_count']
        self.effective_element_count = self.element_count + self.sym_element_count
        cfg.config['model']['mesh_reconstruction']['effective_element_count'] = self.effective_element_count
        self.implicit_parameter_length = cfg.config['model']['mesh_reconstruction']['implicit_parameter_length']
        self.element_embedding_length = 10 + self.implicit_parameter_length
        cfg.config['model']['mesh_reconstruction']['analytic_code_len'] = 10 * self.element_count
        cfg.config['model']['mesh_reconstruction']['structured_implicit_vector_len'] = \
            self.element_embedding_length * self.element_count
        self._temp_folder = None

        '''Modules'''
        self.encoder = resnet.resnet18_full(pretrained=False, num_classes=self.bottleneck_size,
                                            input_channels=4 if cfg.config['data'].get('mask', False) else 3)
        self.mlp = nn.Sequential(
            nn.Linear(self.bottleneck_size + n_classes, self.bottleneck_size), nn.LeakyReLU(0.2, True),
            nn.Linear(self.bottleneck_size, self.bottleneck_size), nn.LeakyReLU(0.2, True),
            nn.Linear(self.bottleneck_size, self.element_count * self.element_embedding_length)
        )
        self.decoder = OccNetDecoder(f_dim=self.implicit_parameter_length)

        # initialize weight
        self.apply(weights_init)

        # initialize resnet
        if pretrained_encoder:
            pretrained_dict = model_zoo.load_url(model_urls['resnet18'])
            model_dict = self.encoder.state_dict()
            if pretrained_dict['conv1.weight'].shape != model_dict['conv1.weight'].shape:
                model_dict['conv1.weight'][:,:3,...] = pretrained_dict['conv1.weight']
                pretrained_dict.pop('conv1.weight')
            pretrained_dict = {k: v for k, v in pretrained_dict.items() if
                               k in model_dict and not k.startswith('fc.')}
            model_dict.update(pretrained_dict)
            self.encoder.load_state_dict(model_dict)

    def eval_implicit_parameters(self, implicit_parameters, samples):
        batch_size, element_count, element_embedding_length = list(implicit_parameters.shape)
        sample_count = samples.shape[-2]
        batched_parameters = torch.reshape(implicit_parameters, [batch_size * element_count, element_embedding_length])
        batched_samples = torch.reshape(samples, [batch_size * element_count, sample_count, -1])
        batched_vals = self.decoder(batched_parameters, batched_samples)
        vals = torch.reshape(batched_vals, [batch_size, element_count, sample_count, 1])
        return vals

    def extract_mesh(self, structured_implicit, resolution=64, extent=0.75, num_samples=10000,
                     cuda=True, marching_cube=True):
        if cuda:
            mesh = []
            for s in structured_implicit.unbind():
                if self._temp_folder is None:
                    self._temp_folder = tempfile.mktemp(dir='/dev/shm')
                    os.makedirs(self._temp_folder)
                    self.decoder.write_occnet_file(os.path.join(self._temp_folder, 'serialized.occnet'))
                    shutil.copy('./external/ldif/ldif2mesh/ldif2mesh', self._temp_folder)
                si_path = os.path.join(self._temp_folder, 'ldif.txt')
                grd_path = os.path.join(self._temp_folder, 'grid.grd')

                s.savetxt(si_path)
                cmd = (f"{os.path.join(self._temp_folder, 'ldif2mesh')} {si_path}" 
                       f" {os.path.join(self._temp_folder, 'serialized.occnet')}"
                       f' {grd_path} -resolution {resolution} -extent {extent}')
                subprocess.check_output(cmd, shell=True)
                _, volume = file_util.read_grd(grd_path)
                _, m = extract_mesh.marching_cubes(volume, extent)
                mesh.append(m)
        else:
            mesh = mesh_util.reconstruction(structured_implicit=structured_implicit, resolution=resolution,
                                            b_min=np.array([-extent] * 3), b_max=np.array([extent] * 3),
                                            use_octree=True, num_samples=num_samples, marching_cube=marching_cube)
        return mesh

    def forward(self, image=None, size_cls=None, samples=None, occnet2gaps=None, structured_implicit=None,
                resolution=None, cuda=True, reconstruction='mesh', apply_class_transfer=True):
        return_dict = {}
        # predict structured_implicit
        return_structured_implicit = structured_implicit
        if isinstance(structured_implicit, dict):
            structured_implicit = StructuredImplicit(config=self.cfg.config, **structured_implicit, net=self)
        elif structured_implicit is None or isinstance(structured_implicit, bool):
            # encoder (ldif.model.model.StructuredImplicitModel.forward)
            # image encoding (ldif.nets.cnn.early_fusion_cnn)
            embedding = self.encoder(image)
            return_dict['ldif_afeature'] = embedding
            embedding = torch.cat([embedding, size_cls], 1)
            structured_implicit_activations = self.mlp(embedding)
            structured_implicit_activations = torch.reshape(
                structured_implicit_activations, [-1, self.element_count, self.element_embedding_length])
            return_dict['structured_implicit_activations'] = structured_implicit_activations

            # SIF decoder
            structured_implicit = StructuredImplicit.from_activation(
                self.cfg.config, structured_implicit_activations, self)
        else:
            raise NotImplementedError

        return_dict['structured_implicit'] = structured_implicit.dict()

        # if only want structured_implicit
        if return_structured_implicit is True:
            return return_dict

        # predict class or mesh
        if samples is not None:
            global_decisions, local_outputs = structured_implicit.class_at_samples(samples, apply_class_transfer)
            return_dict.update({'global_decisions': global_decisions,
                                'element_centers': structured_implicit.centers})
            return return_dict
        elif reconstruction is not None:
            if resolution is None:
                resolution =  self.cfg.config['data'].get('marching_cube_resolution', 128)
            mesh = self.extract_mesh(structured_implicit, extent=self.cfg.config['data']['bounding_box'],
                                     resolution=resolution, cuda=cuda, marching_cube=reconstruction == 'mesh')
            if reconstruction == 'mesh':
                if occnet2gaps is not None:
                    mesh = [m.apply_transform(t.inverse().cpu().numpy()) if not isinstance(m, trimesh.primitives.Sphere) else m
                            for m, t in zip(mesh, occnet2gaps)]

                mesh_coordinates_results = []
                faces = []
                for m in mesh:
                    mesh_coordinates_results.append(
                        torch.from_numpy(m.vertices).type(torch.float32).transpose(-1, -2).to(structured_implicit.device))
                    faces.append(torch.from_numpy(m.faces).to(structured_implicit.device) + 1)
                return_dict.update({'mesh': mesh, 'mesh_coordinates_results': [mesh_coordinates_results, ],
                                    'faces': faces, 'element_centers': structured_implicit.centers})
            elif reconstruction == 'sdf':
                return_dict.update({'sdf': mesh[0], 'mat': mesh[1], 'element_centers': structured_implicit.centers})
            else:
                raise NotImplementedError
            return return_dict
        else:
            return return_dict

    def __del__(self):
        if self._temp_folder is not None:
            shutil.rmtree(self._temp_folder)
