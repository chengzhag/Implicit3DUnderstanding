# Tester for LIEN+LDIF.
# author: chengzhang

import os
from models.testing import BaseTester
import torch
from .training import Trainer
from external.pyTorchChamferDistance.chamfer_distance import ChamferDistance
from libs.tools import write_obj
dist_chamfer = ChamferDistance()
from models.loss import LDIFLoss
from external.ldif.inference.metrics import mesh_chamfer_via_points
from external.ldif.util.file_util import read_mesh
from libs.tools import read_obj, sample_pnts_from_obj, normalize_to_unit_square
import trimesh
import numpy as np
import tempfile
import shutil
import subprocess
from collections import defaultdict
from .modules import LDIF
from ..mgnet.modules import MGNet
import trimesh


class Tester(BaseTester, Trainer):
    '''
    Tester object for SCNet.
    '''
    def __init__(self, cfg, net, device=None):
        super(Tester, self).__init__(cfg, net, device)
        self._temp_folder = None

    def get_metric_values(self, est_data, data):
        losses = defaultdict(list)
        if isinstance(self.net, MGNet):
            final_mesh = est_data['mesh_coordinates_results'][-1].transpose(1, 2)
            mesh_points = data['mesh_points']
        elif isinstance(self.net, LDIF):
            final_mesh = [p.transpose(-1, -2) for p in est_data['mesh_coordinates_results'][-1]]
            mesh_points = np.stack([trimesh.sample.sample_surface(m, 10000)[0] for m in data['mesh']])
            mesh_points = torch.from_numpy(mesh_points).type(torch.float32).to(self.device)  # [10000, 3]
        else:
            raise NotImplementedError

        for index in range(len(final_mesh)):
            dist1, dist2 = dist_chamfer(mesh_points[index].unsqueeze(0), final_mesh[index].unsqueeze(0))[:2]
            losses['Avg_Chamfer'].append(((torch.mean(dist1)) + (torch.mean(dist2))).item())

            if 'mesh' in est_data.keys():
                mesh = est_data['mesh'][index]
            else:
                mesh = trimesh.Trimesh(final_mesh[index].cpu().numpy(), est_data['faces'][index].cpu().numpy() - 1)

            if self.cfg.config['log']['save_results']:
                img_path = data['img_path'][index]
                rel_dir = os.path.join(*img_path.split('/')[-3:])
                rel_dir = os.path.splitext(rel_dir)[0] + '.ply'
                rel_dir = f"{self.cfg.config['log']['vis_path']}/{rel_dir}"
                if not os.path.isdir(os.path.dirname(rel_dir)):
                    os.makedirs(os.path.dirname(rel_dir))
                mesh.export(rel_dir)

            if self._temp_folder is None:
                self._temp_folder = tempfile.mktemp(dir='/dev/shm')
                os.makedirs(self._temp_folder)
                shutil.copy('./external/ldif/gaps/bin/x86_64/mshalign', self._temp_folder)

            if self.cfg.config['full']:
                # ICP mesh alignment
                output_file = os.path.join(self._temp_folder, 'output.ply')
                mesh.export(output_file)
                align_file = os.path.join(self._temp_folder, 'align.ply')
                gt_file = data['mesh_path'][index]
                cmd = f"{os.path.join(self._temp_folder, 'mshalign')} {output_file} {gt_file} {align_file}"
                subprocess.check_output(cmd, shell=True)

                # Chamfer distance
                align_mesh = read_mesh(align_file)
                gt = mesh_points[index].cpu().numpy()
                for ext, output in zip(('woICP', 'wICP'), (mesh, align_mesh)):
                    losses[f'chamfer_{ext}'].append(mesh_chamfer_via_points(points1=output.sample(10000), points2=gt))

        return losses

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

    def __del__(self):
        if self._temp_folder is not None:
            shutil.rmtree(self._temp_folder)
