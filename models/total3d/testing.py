# Tester for Total3D
# author: ynie
# date: April, 2020
import os
from models.testing import BaseTester
from .training import Trainer
import torch
from net_utils.libs import get_rotation_matrix_gt
from configs.data_config import NYU37_TO_PIX3D_CLS_MAPPING
from models.eval_metrics import get_iou_cuboid
from net_utils.libs import get_layout_bdb_sunrgbd, get_rotation_matix_result, get_bdb_evaluation, get_bdb_2d_result, get_corners_of_bb3d_no_index, get_iou
from scipy.io import savemat
from libs.tools import write_obj
import numpy as np
from ..loss import PoseLoss, JointLoss, ReconLoss
from utils.sunrgbd_utils import get_layout_bdb_from_corners, transform_to_world, cvt_R_ex_to_cam_R, unprocess_bdb3d


class Tester(BaseTester, Trainer):
    '''
    Tester object for SCNet.
    '''
    def __init__(self, cfg, net, device=None):
        super(Tester, self).__init__(cfg, net, device)
        self.layout_estimation_loss = PoseLoss()
        self.joint_loss = JointLoss()
        self.mesh_reconstruction_loss = ReconLoss()

    def to_device(self, data):
        data_output = super(Tester, self).to_device(data)
        data_output['sequence_id'] = data['sequence_id']
        if 'world_R_inv' in data['camera'].keys():
            data_output['world_R_inv'] = data['camera']['world_R_inv']
            data_output['lo_inv'] = data['layout']['lo_inv']
            data_output['bdb3d_inv'] = data['boxes_batch']['bdb3d_inv']
        return data_output

    def get_metric_values(self, est_data, gt_data):
        ''' Performs a evaluation step.
        '''
        metrics = {}

        # Layout IoU
        if 'lo_ori_reg_result' in est_data.keys():
            lo_bdb3D_out = get_layout_bdb_sunrgbd(self.cfg.bins_tensor, est_data['lo_ori_reg_result'],
                                                  torch.argmax(est_data['lo_ori_cls_result'], 1), est_data['lo_centroid_result'],
                                                  est_data['lo_coeffs_result'])

            layout_iou = []
            for index, sequence_id in enumerate(gt_data['sequence_id']):
                lo_iou = get_iou_cuboid(lo_bdb3D_out[index, :, :].cpu().numpy(), gt_data['lo_bdb3D'][index, :, :].cpu().numpy())
                layout_iou.append(lo_iou)

            metrics['layout_iou'] = np.mean(layout_iou)

            # camera orientation for evaluation
            cam_R_out, pitch, roll = get_rotation_matix_result(
                self.cfg.bins_tensor,
                torch.argmax(est_data['pitch_cls_result'], 1), est_data['pitch_reg_result'],
                torch.argmax(est_data['roll_cls_result'], 1), est_data['roll_reg_result'],
                return_degrees=True)
            cam_pitch_err = torch.abs(pitch - gt_data['cam_pitch_gt']).mean().item()
            cam_roll_err = torch.abs(roll - gt_data['cam_roll_gt']).mean().item()
            metrics['cam_pitch_err'] = cam_pitch_err / np.pi * 180
            metrics['cam_roll_err'] = cam_roll_err / np.pi * 180

        # projected center
        if 'offset_2D_result' in est_data.keys():
            P_result = torch.stack(((gt_data['bdb2D_pos'][:, 0] + gt_data['bdb2D_pos'][:, 2]) / 2 - (
                    gt_data['bdb2D_pos'][:, 2] - gt_data['bdb2D_pos'][:, 0]) * est_data['offset_2D_result'][:, 0],
                                    (gt_data['bdb2D_pos'][:, 1] + gt_data['bdb2D_pos'][:, 3]) / 2 - (
                                            gt_data['bdb2D_pos'][:, 3] - gt_data['bdb2D_pos'][:, 1]) * est_data['offset_2D_result'][:, 1]), 1)

            bdb3D_out_form_cpu, bdb3D_out = get_bdb_evaluation(self.cfg.bins_tensor, torch.argmax(est_data['ori_cls_result'], 1), est_data['ori_reg_result'],
                                                           torch.argmax(est_data['centroid_cls_result'], 1), est_data['centroid_reg_result'],
                                                           gt_data['size_cls'], est_data['size_reg_result'], P_result, gt_data['K'], cam_R_out, gt_data['split'], return_bdb=True)

            bdb2D_out = get_bdb_2d_result(bdb3D_out, cam_R_out, gt_data['K'], gt_data['split'])

            nyu40class_ids = []
            IoU3D = []
            IoU2D = []
            for index, evaluate_bdb in enumerate(bdb3D_out_form_cpu):
                NYU40CLASS_ID = int(evaluate_bdb['classid'])
                iou_3D = get_iou_cuboid(get_corners_of_bb3d_no_index(evaluate_bdb['basis'], evaluate_bdb['coeffs'],
                                                                     evaluate_bdb['centroid']),
                                        gt_data['bdb3D'][index, :, :].cpu().numpy())

                box1 = bdb2D_out[index, :].cpu().numpy()
                box2 = gt_data['bdb2D_from_3D_gt'][index, :].cpu().numpy()

                box1 = {'u1': box1[0], 'v1': box1[1], 'u2': box1[2], 'v2': box1[3]}
                box2 = {'u1': box2[0], 'v1': box2[1], 'u2': box2[2], 'v2': box2[3]}

                iou_2D = get_iou(box1, box2)

                nyu40class_ids.append(NYU40CLASS_ID)
                IoU3D.append(iou_3D)
                IoU2D.append(iou_2D)
            metrics['iou_3d'] = IoU3D
            metrics['iou_2d'] = IoU2D

        # Lg
        if self.cfg.config['full'] and 'meshes' in est_data.keys():
            mesh_data = est_data.copy()
            if isinstance(est_data['meshes'], list):
                mesh_data['meshes'] = [est_data['meshes'][i] for i in gt_data['mask_status'].nonzero()[0]]
            else:
                mesh_data['meshes'] = est_data['meshes'][gt_data['mask_status'].nonzero()]
            mesh_loss = self.mesh_reconstruction_loss(mesh_data, gt_data, {'cam_R_result': cam_R_out, 'bdb3D_result': bdb3D_out})
            metrics['Lg'] = getattr(mesh_loss['mesh_loss'], 'item', lambda :mesh_loss['mesh_loss'])()

        '''Save results'''
        if self.cfg.config['log']['save_results'] \
                and 'lo_ori_reg_result' in est_data.keys() \
                and 'offset_2D_result' in est_data.keys():

            save_path = self.cfg.config['log']['vis_path']

            for index, sequence_id in enumerate(gt_data['sequence_id']):
                save_path_per_img = os.path.join(save_path, str(sequence_id.item()))
                if not os.path.exists(save_path_per_img):
                    os.mkdir(save_path_per_img)

                # save layout results
                savemat(os.path.join(save_path_per_img, 'layout.mat'),
                        mdict={'layout': lo_bdb3D_out[index, :, :].cpu().numpy()})

                # save bounding boxes and camera poses
                interval = gt_data['split'][index].cpu().tolist()
                current_cls = nyu40class_ids[interval[0]:interval[1]]

                savemat(os.path.join(save_path_per_img, 'bdb_3d.mat'),
                        mdict={'bdb': bdb3D_out_form_cpu[interval[0]:interval[1]], 'class_id': current_cls})
                savemat(os.path.join(save_path_per_img, 'r_ex.mat'),
                        mdict={'cam_R': cam_R_out[index, :, :].cpu().numpy()})

                # save results in cooperative_scene_parsing coordinate
                if 'world_R_inv' in gt_data:
                    world_R_inv = gt_data['world_R_inv'][index].cpu().numpy()
                    lo_inv = gt_data['lo_inv'][index].cpu().numpy()
                    bdb3d_inv = gt_data['bdb3d_inv'][interval[0]:interval[1]].cpu().numpy()
                    layout_3D = get_layout_bdb_from_corners(lo_bdb3D_out[index, :, :].cpu().numpy())
                    bdb3ds_ws = bdb3D_out_form_cpu[interval[0]:interval[1]]
                    cam_R = cam_R_out[index, :, :].cpu().numpy()

                    layout_3D_co, bdb3ds_ws_co, cam_R_co = transform_to_world(layout_3D, bdb3ds_ws, cam_R, world_R_inv)

                    cam_R_co = cvt_R_ex_to_cam_R(cam_R_co)

                    layout_3D_co['basis'] = np.matmul(lo_inv, layout_3D_co['basis'])
                    layout_3D_co['coeffs'] = np.matmul(lo_inv, layout_3D_co['coeffs'])
                    layout_3D_co = get_corners_of_bb3d_no_index(
                        layout_3D_co['basis'], layout_3D_co['coeffs'], layout_3D_co['centroid'])
                    trans_mat = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
                    layout_3D_co = (trans_mat.dot(layout_3D_co.T)).T

                    bdb3ds_ws_co = unprocess_bdb3d(bdb3ds_ws_co, bdb3d_inv)

                    # save layout results
                    savemat(os.path.join(save_path_per_img, 'co_layout.mat'),
                            mdict={'layout': layout_3D_co})
                    savemat(os.path.join(save_path_per_img, 'co_bdb_3d.mat'),
                            mdict={'bdb': bdb3ds_ws_co, 'class_id': current_cls})
                    savemat(os.path.join(save_path_per_img, 'co_r_ex.mat'),
                            mdict={'cam_R': cam_R_co})

                if est_data['meshes'] is not None:
                    current_faces = est_data['out_faces'][interval[0]:interval[1]]
                    current_coordinates = est_data['meshes'][interval[0]:interval[1]]

                    for obj_id, obj_cls in enumerate(current_cls):
                        file_path = os.path.join(save_path_per_img, '%s_%s.obj' % (obj_id, obj_cls))

                        mesh_obj = {'v': current_coordinates[obj_id].transpose(-1,-2).cpu().numpy(),
                                    'f': current_faces[obj_id].cpu().numpy()}

                        write_obj(file_path, mesh_obj)

        return metrics

    def test_step(self, data):
        '''
        test by epoch
        '''
        '''load input and ground-truth data'''
        data = self.to_device(data)

        '''network forwarding'''
        est_data = self.net(data)

        loss = self.get_metric_values(est_data, data)
        return loss

    def visualize_step(self, epoch, phase, iter, data):
        ''' Performs a visualization step.
        '''
        pass
