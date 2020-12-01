# loss function library.
# author: ynie
# date: Feb, 2020
import torch
import torch.nn as nn
import torch.nn.functional as F
from configs.data_config import cls_reg_ratio
from external.pyTorchChamferDistance.chamfer_distance import ChamferDistance
from models.registers import LOSSES
from net_utils.libs import get_layout_bdb_sunrgbd, get_bdb_form_from_corners, \
    recover_points_to_world_sys, get_rotation_matix_result, get_bdb_3d_result, \
    get_bdb_2d_result, physical_violation, recover_points_to_obj_sys
import numpy as np
dist_chamfer = ChamferDistance()

cls_criterion = nn.CrossEntropyLoss(reduction='mean')
reg_criterion = nn.SmoothL1Loss(reduction='mean')
mse_criterion = nn.MSELoss(reduction='mean')
binary_cls_criterion = nn.BCEWithLogitsLoss(reduction='mean')


def cls_reg_loss(cls_result, cls_gt, reg_result, reg_gt):
    cls_loss = cls_criterion(cls_result, cls_gt)
    if len(reg_result.size()) == 3:
        reg_result = torch.gather(reg_result, 1, cls_gt.view(reg_gt.size(0), 1, 1).expand(reg_gt.size(0), 1, reg_gt.size(1)))
    else:
        reg_result = torch.gather(reg_result, 1, cls_gt.view(reg_gt.size(0), 1).expand(reg_gt.size(0), 1))
    reg_result = reg_result.squeeze(1)
    reg_loss = reg_criterion(reg_result, reg_gt)
    return cls_loss, cls_reg_ratio * reg_loss

def get_point_loss(points_in_world_sys, cam_R, cam_K, depth_maps, bdb3D_form, split, obj_masks, mask_status):
    '''
    get the depth loss for each mesh.
    :param points_in_world_sys: Number_of_objects x Number_of_points x 3
    :param cam_R: Number_of_scenes x 3 x 3
    :param cam_K: Number_of_scenes x 3 x 3
    :param depth_maps: Number_of_scenes depth maps in a list
    :param split: Number_of_scenes x 2 matrix
    :return: depth loss
    '''
    depth_loss = 0.
    n_objects = 0
    masked_object_id = -1

    device = cam_R.device

    for scene_id, obj_interval in enumerate(split):
        # map depth to 3d points in camera system.
        u, v = torch.meshgrid(torch.arange(0, depth_maps[scene_id].size(1)), torch.arange(0, depth_maps[scene_id].size(0)))
        u = u.t().to(device)
        v = v.t().to(device)
        u = u.reshape(-1)
        v = v.reshape(-1)
        z_cam = depth_maps[scene_id][v, u]
        u = u.float()
        v = v.float()

        # non_zero_indices = torch.nonzero(z_cam).t()[0]
        # z_cam = z_cam[non_zero_indices]
        # u = u[non_zero_indices]
        # v = v[non_zero_indices]

        # calculate coordinates
        x_cam = (u - cam_K[scene_id][0][2])*z_cam/cam_K[scene_id][0][0]
        y_cam = (v - cam_K[scene_id][1][2])*z_cam/cam_K[scene_id][1][1]

        # transform to toward-up-right coordinate system
        points_world = torch.cat([z_cam.unsqueeze(-1), -y_cam.unsqueeze(-1), x_cam.unsqueeze(-1)], -1)
        # transform from camera system to world system
        points_world = torch.mm(points_world, cam_R[scene_id].t())

        n_columns = depth_maps[scene_id].size(1)

        for loc_id, obj_id in enumerate(range(*obj_interval)):
            if mask_status[obj_id] == 0:
                continue
            masked_object_id += 1

            bdb2d = obj_masks[scene_id][loc_id]['msk_bdb']
            obj_msk = obj_masks[scene_id][loc_id]['msk']

            u_s, v_s = torch.meshgrid(torch.arange(bdb2d[0], bdb2d[2] + 1), torch.arange(bdb2d[1], bdb2d[3] + 1))
            u_s = u_s.t().long()
            v_s = v_s.t().long()
            index_dep = u_s + n_columns * v_s
            index_dep = index_dep.reshape(-1)
            in_object_indices = obj_msk.reshape(-1).nonzero()[0]

            # remove holes in depth maps
            if len(in_object_indices) == 0:
                continue

            object_pnts = points_world[index_dep,:][in_object_indices,:]
            # remove noisy points that out of bounding boxes
            inner_idx = torch.sum(torch.abs(
                torch.mm(object_pnts - bdb3D_form['centroid'][masked_object_id].view(1, 3), bdb3D_form['basis'][masked_object_id].t())) >
                                  bdb3D_form['coeffs'][masked_object_id], dim=1)

            inner_idx = torch.nonzero(inner_idx == 0).t()[0]

            if inner_idx.nelement() == 0:
                continue

            object_pnts = object_pnts[inner_idx, :]

            dist_1 = dist_chamfer(object_pnts.unsqueeze(0), points_in_world_sys[masked_object_id].unsqueeze(0))[0]
            depth_loss += torch.mean(dist_1)
            n_objects += 1
    return depth_loss/n_objects if n_objects > 0 else torch.tensor(0.).to(device)

class BaseLoss(object):
    '''base loss class'''
    def __init__(self, weight=1, config=None):
        '''initialize loss module'''
        self.weight = weight
        self.config = config

@LOSSES.register_module
class Null(BaseLoss):
    '''This loss function is for modules where a loss preliminary calculated.'''
    def __call__(self, loss):
        return self.weight * torch.mean(loss)

@LOSSES.register_module
class SVRLoss(BaseLoss):
    def __call__(self, est_data, gt_data, subnetworks, face_sampling_rate):
        device = est_data['mesh_coordinates_results'][0].device
        # chamfer losses
        chamfer_loss = torch.tensor(0.).to(device)
        edge_loss = torch.tensor(0.).to(device)
        boundary_loss = torch.tensor(0.).to(device)

        for stage_id, mesh_coordinates_result in enumerate(est_data['mesh_coordinates_results']):
            mesh_coordinates_result = mesh_coordinates_result.transpose(1, 2)
            # points to points chamfer loss
            dist1, dist2 = dist_chamfer(gt_data['mesh_points'], mesh_coordinates_result)[:2]
            chamfer_loss += (torch.mean(dist1)) + (torch.mean(dist2))

            # boundary loss
            if stage_id == subnetworks - 1:
                if 1 in est_data['boundary_point_ids']:
                    boundary_loss = torch.mean(dist2[est_data['boundary_point_ids']])

            # edge loss
            edge_vec = torch.gather(mesh_coordinates_result, 1,
                                    (est_data['output_edges'][:, :, 0] - 1).unsqueeze(-1).expand(est_data['output_edges'].size(0),
                                                                                     est_data['output_edges'].size(1), 3)) \
                       - torch.gather(mesh_coordinates_result, 1,
                                      (est_data['output_edges'][:, :, 1] - 1).unsqueeze(-1).expand(est_data['output_edges'].size(0),
                                                                                       est_data['output_edges'].size(1), 3))

            edge_vec = edge_vec.view(edge_vec.size(0) * edge_vec.size(1), edge_vec.size(2))
            edge_loss += torch.mean(torch.pow(torch.norm(edge_vec, p=2, dim=1), 2))

        chamfer_loss = 100 * chamfer_loss / len(est_data['mesh_coordinates_results'])
        edge_loss = 100 * edge_loss / len(est_data['mesh_coordinates_results'])
        boundary_loss = 100 * boundary_loss

        # face distance losses
        face_loss = torch.tensor(0.).to(device)
        for points_from_edges_by_step, points_indicator_by_step in zip(est_data['points_from_edges'], est_data['point_indicators']):
            points_from_edges_by_step = points_from_edges_by_step.transpose(1, 2).contiguous()
            _, dist2_face, _, idx2 = dist_chamfer(gt_data['mesh_points'], points_from_edges_by_step)
            idx2 = idx2.long()
            dist2_face = dist2_face.view(dist2_face.shape[0], dist2_face.shape[1] // face_sampling_rate,
                                         face_sampling_rate)

            # average distance to nearest face.
            dist2_face = torch.mean(dist2_face, dim=2)
            local_dens = gt_data['densities'][:, idx2[:]][range(gt_data['densities'].size(0)), range(gt_data['densities'].size(0)), :]
            in_mesh = (dist2_face <= local_dens).float()
            face_loss += binary_cls_criterion(points_indicator_by_step, in_mesh)

        if est_data['points_from_edges']:
            face_loss = face_loss / len(est_data['points_from_edges'])

        return {'chamfer_loss': chamfer_loss, 'face_loss': 0.01 * face_loss,
                'edge_loss': 0.1 * edge_loss, 'boundary_loss': 0.5 * boundary_loss}

@LOSSES.register_module
class PoseLoss(BaseLoss):
    def __call__(self, est_data, gt_data, bins_tensor):
        pitch_cls_loss, pitch_reg_loss = cls_reg_loss(est_data['pitch_cls_result'], gt_data['pitch_cls'], est_data['pitch_reg_result'], gt_data['pitch_reg'])
        roll_cls_loss, roll_reg_loss = cls_reg_loss(est_data['roll_cls_result'], gt_data['roll_cls'], est_data['roll_reg_result'], gt_data['roll_reg'])
        lo_ori_cls_loss, lo_ori_reg_loss = cls_reg_loss(est_data['lo_ori_cls_result'], gt_data['lo_ori_cls'], est_data['lo_ori_reg_result'], gt_data['lo_ori_reg'])
        lo_centroid_loss = reg_criterion(est_data['lo_centroid_result'], gt_data['lo_centroid']) * cls_reg_ratio
        lo_coeffs_loss = reg_criterion(est_data['lo_coeffs_result'], gt_data['lo_coeffs']) * cls_reg_ratio

        lo_bdb3D_result = get_layout_bdb_sunrgbd(bins_tensor, est_data['lo_ori_reg_result'], gt_data['lo_ori_cls'], est_data['lo_centroid_result'],
                                                 est_data['lo_coeffs_result'])
        # layout bounding box corner loss
        lo_corner_loss = cls_reg_ratio * reg_criterion(lo_bdb3D_result, gt_data['lo_bdb3D'])

        return {'pitch_cls_loss':pitch_cls_loss, 'pitch_reg_loss':pitch_reg_loss,
                'roll_cls_loss':roll_cls_loss, 'roll_reg_loss':roll_reg_loss,
                'lo_ori_cls_loss':lo_ori_cls_loss, 'lo_ori_reg_loss':lo_ori_reg_loss,
                'lo_centroid_loss':lo_centroid_loss, 'lo_coeffs_loss':lo_coeffs_loss,
                'lo_corner_loss':lo_corner_loss}, {'lo_bdb3D_result':lo_bdb3D_result}

@LOSSES.register_module
class DetLoss(BaseLoss):
    def __call__(self, est_data, gt_data):
        # calculate loss
        size_reg_loss = reg_criterion(est_data['size_reg_result'], gt_data['size_reg']) * cls_reg_ratio
        ori_cls_loss, ori_reg_loss = cls_reg_loss(est_data['ori_cls_result'], gt_data['ori_cls'], est_data['ori_reg_result'], gt_data['ori_reg'])
        centroid_cls_loss, centroid_reg_loss = cls_reg_loss(est_data['centroid_cls_result'], gt_data['centroid_cls'],
                                                          est_data['centroid_reg_result'], gt_data['centroid_reg'])
        offset_2D_loss = reg_criterion(est_data['offset_2D_result'], gt_data['offset_2D'])

        return {'size_reg_loss':size_reg_loss, 'ori_cls_loss':ori_cls_loss, 'ori_reg_loss':ori_reg_loss,
                'centroid_cls_loss':centroid_cls_loss, 'centroid_reg_loss':centroid_reg_loss,
                'offset_2D_loss':offset_2D_loss}

@LOSSES.register_module
class ReconLoss(BaseLoss):
    def __call__(self, est_data, gt_data, extra_results):
        if gt_data['mask_flag'] == 0:
            point_loss = 0.
        else:
            # get the world coordinates for each 3d object.
            bdb3D_form = get_bdb_form_from_corners(extra_results['bdb3D_result'], gt_data['mask_status'])
            obj_points_in_world_sys = recover_points_to_world_sys(bdb3D_form, est_data['meshes'])
            point_loss = 100 * get_point_loss(obj_points_in_world_sys, extra_results['cam_R_result'],
                                              gt_data['K'], gt_data['depth_maps'], bdb3D_form, gt_data['split'],
                                              gt_data['obj_masks'], gt_data['mask_status'])

            # remove samples without depth map
            if torch.isnan(point_loss):
                point_loss = 0.

        return {'mesh_loss':point_loss}

@LOSSES.register_module
class JointLoss(BaseLoss):
    def __call__(self, est_data, gt_data, bins_tensor, layout_results):
        # predicted camera rotation
        cam_R_result = get_rotation_matix_result(bins_tensor,
                                                 gt_data['pitch_cls'], est_data['pitch_reg_result'],
                                                 gt_data['roll_cls'], est_data['roll_reg_result'])

        # projected center
        P_result = torch.stack(
            ((gt_data['bdb2D_pos'][:, 0] + gt_data['bdb2D_pos'][:, 2]) / 2 - (gt_data['bdb2D_pos'][:, 2] - gt_data['bdb2D_pos'][:, 0]) * est_data['offset_2D_result'][:, 0],
             (gt_data['bdb2D_pos'][:, 1] + gt_data['bdb2D_pos'][:, 3]) / 2 - (gt_data['bdb2D_pos'][:, 3] - gt_data['bdb2D_pos'][:, 1]) * est_data['offset_2D_result'][:, 1]), 1)

        # retrieved 3D bounding box
        bdb3D_result, _ = get_bdb_3d_result(bins_tensor,
                                            gt_data['ori_cls'],
                                            est_data['ori_reg_result'],
                                            gt_data['centroid_cls'],
                                            est_data['centroid_reg_result'],
                                            gt_data['size_cls'],
                                            est_data['size_reg_result'],
                                            P_result,
                                            gt_data['K'],
                                            cam_R_result,
                                            gt_data['split'])

        # 3D bounding box corner loss
        corner_loss = 5 * cls_reg_ratio * reg_criterion(bdb3D_result, gt_data['bdb3D'])

        # 2D bdb loss
        bdb2D_result = get_bdb_2d_result(bdb3D_result, cam_R_result, gt_data['K'], gt_data['split'])
        bdb2D_loss = 20 * cls_reg_ratio * reg_criterion(bdb2D_result, gt_data['bdb2D_from_3D_gt'])

        # physical violation loss
        phy_violation, phy_gt = physical_violation(layout_results['lo_bdb3D_result'], bdb3D_result, gt_data['split'])
        phy_loss = 20 * mse_criterion(phy_violation, phy_gt)

        return {'phy_loss':phy_loss, 'bdb2D_loss':bdb2D_loss, 'corner_loss':corner_loss},\
               {'cam_R_result':cam_R_result, 'bdb3D_result':bdb3D_result}

@LOSSES.register_module
class LDIFLoss(BaseLoss):
    def __call__(self, est_data, gt_data):
        # calculate loss (ldif.training.loss.compute_loss)
        uniform_sample_loss = nn.MSELoss()(est_data['uniform_class'], gt_data['uniform_class'])
        uniform_sample_loss *= self.config['model']['mesh_reconstruction']['uniform_loss_weight']

        near_surface_sample_loss = nn.MSELoss()(est_data['near_surface_class'], gt_data['near_surface_class'])
        near_surface_sample_loss *= self.config['model']['mesh_reconstruction']['near_surface_loss_weight']

        element_centers = est_data['element_centers']
        xyzw_samples = F.pad(element_centers, [0, 1], "constant", 1) # 维度为 [batch_size, sample_count, 4]
        xyzw_samples = torch.matmul(xyzw_samples, gt_data['world2grid'])[..., :3] # 维度为 [batch_size, sample_count, 3]
        # xyzw_samples = torch.matmul(torch.Tensor([[[1.2] * 3 + [1]], ] * xyzw_samples.shape[0]).to(element_centers.device), gt_data['world2grid'])[..., :3] # 测试边界
        grid = gt_data['grid']
        scale_fac = torch.Tensor(list(grid.shape)[1:]).to(element_centers.device) / 2 - 0.5
        xyzw_samples /= scale_fac
        xyzw_samples = xyzw_samples.unsqueeze(1).unsqueeze(1) # 维度为 [batch_size, 1, 1, sample_count, 3]
        grid = grid.unsqueeze(1)
        gt_sdf_at_centers = F.grid_sample(grid, xyzw_samples, mode='bilinear', padding_mode='zeros') # https://pytorch.org/docs/master/nn.functional.html#torch.nn.functional.grid_sample
        gt_sdf_at_centers = torch.where(gt_sdf_at_centers > self.config['data']['coarse_grid_spacing'] / 1.1,
                                        gt_sdf_at_centers, torch.zeros(1).to(gt_sdf_at_centers.device))
        gt_sdf_at_centers *= self.config['model']['mesh_reconstruction']['lowres_grid_inside_loss_weight']
        element_center_lowres_grid_inside_loss = torch.mean((gt_sdf_at_centers + 1e-04) ** 2) + 1e-05

        bounding_box = self.config['data']['bounding_box']
        lower, upper = -bounding_box, bounding_box
        lower_error = torch.max(lower - element_centers, torch.zeros(1).cuda())
        upper_error = torch.max(element_centers - upper, torch.zeros(1).cuda())
        bounding_box_constraint_error = lower_error * lower_error + upper_error * upper_error
        bounding_box_error = torch.mean(bounding_box_constraint_error)
        inside_box_loss = self.config['model']['mesh_reconstruction']['inside_box_loss_weight'] * bounding_box_error

        return {'uniform_sample_loss': uniform_sample_loss,
                'near_surface_sample_loss': near_surface_sample_loss,
                'fixed_bounding_box_loss': inside_box_loss,
                'lowres_grid_inside_loss':element_center_lowres_grid_inside_loss}


def get_phy_loss_samples(ldif, structured_implicit, ldif_center, ldif_coef, phy_loss_samples,
                         return_range=False, surface_optimize=False):
    # get inside points from blob centers
    centers = structured_implicit.all_centers.clone()
    sample_points = centers

    # get inside points with random sampling
    bbox_samples = (torch.rand([len(centers), phy_loss_samples * (3 if surface_optimize else 10), 3],
                               device=centers.device) - 0.5) * 2 * ldif_coef.unsqueeze(1) + ldif_center.unsqueeze(1)
    sample_points = torch.cat([sample_points, bbox_samples], 1)

    # optimize to get surface points
    if surface_optimize:
        surface_samples = sample_points.clone()
        surface_samples.requires_grad = True
        ldif_grad = []
        for param in ldif.parameters():
            ldif_grad.append(param.requires_grad)
            param.requires_grad = True
        optimizer = torch.optim.SGD([surface_samples], 200)
        with torch.enable_grad():
            for i in range(10):
                optimizer.zero_grad()
                est_sdf = ldif(
                    samples=surface_samples,
                    structured_implicit=structured_implicit.dict(),
                    apply_class_transfer=False,
                )['global_decisions'] + 0.07
                # from external.PIFu.lib import sample_util
                # sample_util.save_samples_truncted_prob(f'out/inside_sample_{i}.ply', surface_samples[0].detach().cpu(),
                #                                        (est_sdf[0] < 0).detach().cpu())
                error = torch.mean(abs(est_sdf))
                error.backward()
                optimizer.step()
                # print(f"({i}) sdf error: {error:.4f}")
        for param, requires_grad in zip(ldif.parameters(), ldif_grad):
            param.requires_grad = requires_grad
        surface_samples.requires_grad = False
        sample_points = torch.cat([sample_points, surface_samples], 1)

    # remove outside points
    est_sdf = ldif(
        samples=sample_points,
        structured_implicit=structured_implicit.dict(),
        apply_class_transfer=True,
    )['global_decisions']
    inside_samples = []
    in_coor_min = []
    in_coor_max = []
    for i, (s, mask) in enumerate(zip(sample_points, est_sdf < 0.5)):
        inside_sample = s[mask.squeeze(), :]
        if return_range:
            if len(inside_sample) > 0:
                inside_min = inside_sample.min(0)[0]
                inside_max = inside_sample.max(0)[0]
            if len(inside_sample) <= 0 or (inside_min == inside_max).sum() > 0:
                inside_min = centers[i].min(0)[0]
                inside_max = centers[i].max(0)[0]
            in_coor_min.append(inside_min)
            in_coor_max.append(inside_max)
        # from external.PIFu.lib import sample_util
        # sample_util.save_samples_truncted_prob('out/inside_sample.ply', s.detach().cpu(),
        #                                        mask.detach().cpu())
        if len(inside_sample) <= 0:
            inside_sample = centers[i]
        # from external.PIFu.lib import sample_util
        # sample_util.save_samples_truncted_prob('out/inside_sample.ply', inside_sample.detach().cpu(),
        #                                        np.zeros(inside_sample.shape[0]))
        p_ids = np.random.choice(len(inside_sample), phy_loss_samples, replace=True)
        inside_sample = inside_sample[p_ids]
        inside_samples.append(inside_sample)
    inside_samples = torch.stack(inside_samples)

    if return_range:
        in_coor_min = torch.stack(in_coor_min)
        in_coor_max = torch.stack(in_coor_max)
        return inside_samples, in_coor_min, in_coor_max
    return inside_samples


@LOSSES.register_module
class LDIFReconLoss(BaseLoss):
    def __call__(self, est_data, gt_data, extra_results):
        ldif = est_data['mgn']
        get_phy_loss = self.config.get('loss_weights', {}).get('ldif_phy_loss', 0.0) > 0
        device = gt_data['image'].device

        loss_settings = self.config['model']['mesh_reconstruction'].get('loss_settings', {})
        loss_type = loss_settings.get('type', 'classmse')
        scale_before_func = loss_settings.get('scale_before_func', 100)
        phy_loss_samples = loss_settings.get('phy_loss_samples', 128)
        phy_loss_objects = loss_settings.get('phy_loss_objects', 128)
        surface_optimize = self.config['model']['mesh_reconstruction']['loss_settings'].get('surface_optimize', False)
        sdf_data = {}

        if get_phy_loss:
            bdb3D_form = get_bdb_form_from_corners(extra_results['bdb3D_result'])
            structured_implicit = est_data['structured_implicit']
            ldif_center, ldif_coef = est_data['obj_center'], est_data['obj_coef']
            if 'ldif_sampling_points' in est_data:
                inside_samples = est_data['ldif_sampling_points']
            else:
                obj_center = ldif_center.clone()
                obj_center[:, 2] *= -1
                inside_samples = get_phy_loss_samples(ldif, structured_implicit, obj_center, ldif_coef,
                                                      phy_loss_samples, surface_optimize=surface_optimize)

            # put points to other objects' coor
            inside_samples[:, :, 2] *= -1
            obj_samples = recover_points_to_world_sys(bdb3D_form, inside_samples, ldif_center, ldif_coef)
            max_sample_points = (gt_data['split'][:, 1] - gt_data['split'][:, 0] - 1).max() * obj_samples.shape[1]
            if max_sample_points == 0:
                sdf_data['ldif_phy_loss'] = None
            else:
                est_sdf = []
                for start, end in gt_data['split']:
                    assert end > start
                    if end > start + 1:
                        centroids = bdb3D_form['centroid'][start:end]
                        centroids = centroids.unsqueeze(0).expand(len(centroids), -1, -1)
                        distances = F.pairwise_distance(
                            centroids.reshape(-1, 3), centroids.transpose(0, 1).reshape(-1, 3), 2
                        ).reshape(len(centroids), len(centroids))
                        for obj_ind in range(start, end):
                            other_obj_dis = distances[obj_ind - start]
                            _, nearest = torch.sort(other_obj_dis)
                            other_obj_sample = obj_samples[start:end].index_select(
                                0, nearest[1:phy_loss_objects + 1]).reshape(-1, 3)
                            other_obj_sample = recover_points_to_obj_sys(
                                {k: v[obj_ind:obj_ind + 1] for k, v in bdb3D_form.items()},
                                other_obj_sample.unsqueeze(0),
                                ldif_center[obj_ind:obj_ind + 1],
                                ldif_coef[obj_ind:obj_ind + 1]
                            )
                            other_obj_sample[:, :, 2] *= -1
                            sdf = ldif(
                                samples=other_obj_sample,
                                structured_implicit=structured_implicit[obj_ind:obj_ind + 1].dict(),
                                apply_class_transfer=False,
                            )['global_decisions']
                            est_sdf.append(sdf.squeeze())
                if len(est_sdf) == 0:
                    sdf_data['ldif_phy_loss'] = None
                else:
                    est_sdf = torch.cat(est_sdf) + 0.07
                    est_sdf[est_sdf > 0] = 0
                    gt_sdf = torch.full(est_sdf.shape, 0., device=device, dtype=torch.float32)
                    sdf_data['ldif_phy_loss'] = (est_sdf, gt_sdf)

        # compute final loss
        loss = {}
        if not isinstance(loss_type, list):
            loss_type = [loss_type] * len(sdf_data)
        for lt, (k, sdf) in zip(loss_type, sdf_data.items()):
            if sdf is None:
                loss[k] = 0.
            else:
                est_sdf, gt_sdf = sdf

                if 'class' in lt:
                    est_sdf = torch.sigmoid(scale_before_func * est_sdf) - 0.5
                    gt_sdf[gt_sdf > 0] = 0.5
                    gt_sdf[gt_sdf < 0] = -0.5
                elif 'sdf' in lt:
                    est_sdf = scale_before_func * est_sdf
                else:
                    raise NotImplementedError

                if 'mse' in lt:
                    point_loss = nn.MSELoss()(est_sdf, gt_sdf)
                elif 'l1' in lt:
                    point_loss = nn.L1Loss()(est_sdf, gt_sdf)
                elif 'sl1' in lt:
                    point_loss = nn.SmoothL1Loss()(est_sdf, gt_sdf)
                else:
                    raise NotImplementedError

                loss[k] = point_loss

        return loss
