# Total3D: model loader
# author: ynie
# date: Feb, 2020

from models.registers import METHODS, MODULES, LOSSES
from models.network import BaseNetwork
import torch
from torch import nn
from configs.data_config import obj_cam_ratio
from external.ldif.representation.structured_implicit_function import StructuredImplicit
import numpy as np
from models.loss import get_phy_loss_samples


@METHODS.register_module
class TOTAL3D(BaseNetwork):

    def __init__(self, cfg):
        '''
        load submodules for the network.
        :param config: customized configurations.
        '''
        super(BaseNetwork, self).__init__()
        self.cfg = cfg

        phase_names = []
        if cfg.config[cfg.config['mode']]['phase'] in ['layout_estimation', 'joint']:
            phase_names += ['layout_estimation']
        if cfg.config[cfg.config['mode']]['phase'] in ['object_detection', 'joint']:
            phase_names += ['object_detection']
        if cfg.config[cfg.config['mode']]['phase'] in ['joint']:
            phase_names += ['mesh_reconstruction']
            if 'output_adjust' in cfg.config['model'].keys():
                phase_names += ['output_adjust']

        if (not cfg.config['model']) or (not phase_names):
            cfg.log_string('No submodule found. Please check the phase name and model definition.')
            raise ModuleNotFoundError('No submodule found. Please check the phase name and model definition.')

        '''load network blocks'''
        for phase_name in phase_names:
            if phase_name not in cfg.config['model'].keys():
                continue
            net_spec = cfg.config['model'][phase_name]
            method_name = net_spec['method']
            # load specific optimizer parameters
            optim_spec = self.load_optim_spec(cfg.config, net_spec)
            subnet = MODULES.get(method_name)(cfg, optim_spec)
            self.add_module(phase_name, subnet)

            '''load corresponding loss functions'''
            setattr(self, phase_name + '_loss', LOSSES.get(self.cfg.config['model'][phase_name]['loss'], 'Null')(
                self.cfg.config['model'][phase_name].get('weight', 1), cfg.config))

        '''Add joint loss'''
        setattr(self, 'joint_loss', LOSSES.get('JointLoss', 'Null')(1))

        '''Multi-GPU setting'''
        # Note that for object_detection, we should extract relational features, thus it does not support parallel training.
        if cfg.config[cfg.config['mode']]['phase'] in ['layout_estimation', 'joint']:
            self.layout_estimation = nn.DataParallel(self.layout_estimation)
        if cfg.config[cfg.config['mode']]['phase'] in ['joint']:
            self.mesh_reconstruction = nn.DataParallel(self.mesh_reconstruction)

        '''freeze submodules or not'''
        self.freeze_modules(cfg)

    def get_extra_results(self, all_output):
        # center and coef of the reconstructed mesh
        extra_results = {}
        if isinstance(self.mesh_reconstruction.module, MODULES.get('LDIF')):
            structured_implicit = all_output['structured_implicit']
            in_coor_min = structured_implicit.all_centers.min(dim=1)[0]
            in_coor_max = structured_implicit.all_centers.max(dim=1)[0]

            obj_center = (in_coor_max + in_coor_min) / 2.
            obj_center[:, 2] *= -1
            obj_coef = (in_coor_max - in_coor_min) / 2.
        else:
            if all_output['meshes'] is None:
                return {}
            obj_center = (all_output['meshes'].max(dim=0)[0] + all_output['meshes'].min(dim=0)[0]) / 2.
            obj_center = obj_center.detach()
            obj_coef = (all_output['meshes'].max(dim=0)[0] - all_output['meshes'].min(dim=0)[0]) / 2.
            obj_coef = obj_coef.detach()

        extra_results.update({'obj_center': obj_center, 'obj_coef': obj_coef})
        return extra_results

    def forward(self, data):
        all_output = {}

        if self.cfg.config[self.cfg.config['mode']]['phase'] in ['layout_estimation', 'joint']:
            pitch_reg_result, roll_reg_result, \
            pitch_cls_result, roll_cls_result, \
            lo_ori_reg_result, lo_ori_cls_result, \
            lo_centroid_result, lo_coeffs_result, a_features = self.layout_estimation(data['image'])

            layout_output = {'pitch_reg_result':pitch_reg_result, 'roll_reg_result':roll_reg_result,
                             'pitch_cls_result':pitch_cls_result, 'roll_cls_result':roll_cls_result,
                             'lo_ori_reg_result':lo_ori_reg_result, 'lo_ori_cls_result':lo_ori_cls_result,
                             'lo_centroid_result':lo_centroid_result, 'lo_coeffs_result':lo_coeffs_result,
                             'lo_afeatures': a_features}
            all_output.update(layout_output)

        if self.cfg.config[self.cfg.config['mode']]['phase'] in ['object_detection', 'joint']:
            size_reg_result, \
            ori_reg_result, ori_cls_result, \
            centroid_reg_result, centroid_cls_result, \
            offset_2D_result, a_features, \
            r_features, a_r_features = self.object_detection(data['patch'], data['size_cls'], data['g_features'],
                                                             data['split'], data['rel_pair_counts'])
            object_output = {'size_reg_result':size_reg_result, 'ori_reg_result':ori_reg_result,
                             'ori_cls_result':ori_cls_result, 'centroid_reg_result':centroid_reg_result,
                             'centroid_cls_result':centroid_cls_result, 'offset_2D_result':offset_2D_result,
                             'odn_afeature': a_features, 'odn_rfeatures': r_features, 'odn_arfeatures': a_r_features}
            all_output.update(object_output)

        if self.cfg.config[self.cfg.config['mode']]['phase'] in ['joint']:
            # predict meshes
            if self.cfg.config['mode']=='train':
                if isinstance(self.mesh_reconstruction.module, MODULES.get('LDIF')):
                    # for LDIF
                    mesh_output = self.mesh_reconstruction(
                        data['patch_for_mesh'], data['cls_codes'], structured_implicit=True)
                else:
                    # for MGNet
                    mesh_output = self.mesh_reconstruction(
                        data['patch_for_mesh'], data['cls_codes'], mask_status=data['mask_status'])
                    if data['mask_flag'] == 1:
                        meshes = mesh_output['mesh_coordinates_results'][-1]
                        mesh_output['meshes'] = meshes
                    else:
                        mesh_output['meshes'] = None
            else:
                # for test
                reconstruction = 'mesh' if self.cfg.config['full'] else None
                mesh_output = self.mesh_reconstruction(data['patch_for_mesh'], data['cls_codes'],
                                                       reconstruction=reconstruction)
                out_points = mesh_output.get('mesh_coordinates_results', [None])
                out_faces = mesh_output.get('faces', None)
                mesh_output.update({'meshes': out_points[-1], 'out_faces': out_faces})

            # construct structured_implicti from dict
            if 'structured_implicit' in mesh_output:
                mesh_output['structured_implicit'] = StructuredImplicit(
                    config=self.cfg.config, **mesh_output['structured_implicit'])

            # convert to SUNRGBD coordinates
            if mesh_output.get('meshes') is not None:
                if isinstance(mesh_output['meshes'], list):
                    for m in mesh_output['meshes']:
                        m[2, :] *= -1
                elif mesh_output['meshes'] is not None:
                    mesh_output['meshes'][:, 2, :] *= -1
            mesh_output['mgn'] = self.mesh_reconstruction
            all_output.update(mesh_output)

            # get extra_results
            all_output.update(self.get_extra_results(all_output))

            if hasattr(self, 'output_adjust'):
                input = all_output.copy()
                input['size_cls'] = data['size_cls']
                input['cls_codes'] = data['cls_codes']
                input['g_features'] = data['g_features']
                input['bdb2D_pos'] = data['bdb2D_pos']
                input['K'] = data['K']
                input['split'] = data['split']
                input['rel_pair_counts'] = data['rel_pair_counts']
                refined_output = self.output_adjust(input)
                all_output.update(refined_output)

        if all_output:
            return all_output
        else:
            raise NotImplementedError

    def loss(self, est_data, gt_data):
        '''
        calculate loss of est_out given gt_out.
        '''
        loss_weights = self.cfg.config.get('loss_weights', {})
        if self.cfg.config[self.cfg.config['mode']]['phase'] in ['layout_estimation', 'joint']:
            layout_loss, layout_results = self.layout_estimation_loss(est_data, gt_data, self.cfg.bins_tensor)
            layout_loss_weighted = {k: v * loss_weights.get(k, 1.0) for k, v in layout_loss.items()}
            total_layout_loss = sum(layout_loss_weighted.values())
            total_layout_loss_unweighted = sum([v.detach() for v in layout_loss.values()])
            for key, value in layout_loss.items():
                layout_loss[key] = value.item()
        if self.cfg.config[self.cfg.config['mode']]['phase'] in ['object_detection', 'joint']:
            object_loss = self.object_detection_loss(est_data, gt_data)
            object_loss_weighted = {k: v * loss_weights.get(k, 1.0) for k, v in object_loss.items()}
            total_object_loss = sum(object_loss_weighted.values())
            total_object_loss_unweighted = sum([v.detach() for v in object_loss.values()])
            for key, value in object_loss.items():
                object_loss[key] = value.item()
        if self.cfg.config[self.cfg.config['mode']]['phase'] in ['joint']:
            joint_loss, extra_results = self.joint_loss(est_data, gt_data, self.cfg.bins_tensor, layout_results)
            joint_loss_weighted = {k: v * loss_weights.get(k, 1.0) for k, v in joint_loss.items()}
            mesh_loss = self.mesh_reconstruction_loss(est_data, gt_data, extra_results)
            mesh_loss_weighted = {k: v * loss_weights.get(k, 1.0) for k, v in mesh_loss.items()}

            total_joint_loss = sum(joint_loss_weighted.values()) + sum(mesh_loss_weighted.values())
            total_joint_loss_unweighted = \
                sum([v.detach() for v in joint_loss.values()]) \
                + sum([v.detach() if isinstance(v, torch.Tensor) else v for v in mesh_loss.values()])
            for key, value in mesh_loss.items():
                mesh_loss[key] = float(value)
            for key, value in joint_loss.items():
                joint_loss[key] = value.item()

        if self.cfg.config[self.cfg.config['mode']]['phase'] == 'layout_estimation':
            return {'total':total_layout_loss, **layout_loss, 'total_unweighted': total_layout_loss_unweighted}
        if self.cfg.config[self.cfg.config['mode']]['phase'] == 'object_detection':
            return {'total':total_object_loss, **object_loss, 'total_unweighted': total_object_loss_unweighted}
        if self.cfg.config[self.cfg.config['mode']]['phase'] == 'joint':
            total3d_loss = total_object_loss + total_joint_loss + obj_cam_ratio * total_layout_loss
            total3d_loss_unweighted = total_object_loss_unweighted + total_joint_loss_unweighted\
                                      + obj_cam_ratio * total_layout_loss_unweighted
            return {'total':total3d_loss, **layout_loss, **object_loss, **mesh_loss, **joint_loss,
                    'total_unweighted': total3d_loss_unweighted}
        else:
            raise NotImplementedError