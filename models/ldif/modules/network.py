# LIEN+LDIF: model loader
# author: chengzhang

from models.registers import METHODS, MODULES, LOSSES
from models.network import BaseNetwork
import torch
from torch import nn


@METHODS.register_module
class LDIF(BaseNetwork):

    def __init__(self, cfg):
        '''
        load submodules for the network.
        :param config: customized configurations.
        '''
        super(BaseNetwork, self).__init__()
        self.cfg = cfg

        '''load network blocks'''
        for phase_name, net_spec in cfg.config['model'].items():
            method_name = net_spec['method']
            # load specific optimizer parameters
            optim_spec = self.load_optim_spec(cfg.config, net_spec)
            subnet = MODULES.get(method_name)(cfg, optim_spec)
            self.add_module(phase_name, subnet)

            '''load corresponding loss functions'''
            setattr(self, phase_name + '_loss', LOSSES.get(self.cfg.config['model'][phase_name]['loss'], 'Null')(
                self.cfg.config['model'][phase_name].get('weight', 1), cfg.config))

        '''Multi-GPU setting'''
        self.mesh_reconstruction = nn.DataParallel(self.mesh_reconstruction)

        '''freeze submodules or not'''
        self.freeze_modules(cfg)

    def forward(self, data):
        if 'uniform_samples' in data.keys():
            samples = torch.cat([data['near_surface_samples'], data['uniform_samples']], 1)
            len_near_surface = data['near_surface_class'].shape[1]
            est_data = self.mesh_reconstruction(data['img'], data['cls'], samples=samples)
            est_data['near_surface_class'] = est_data['global_decisions'][:, :len_near_surface, ...]
            est_data['uniform_class'] = est_data['global_decisions'][:, len_near_surface:, ...]
        else:
            est_data = self.mesh_reconstruction(data['img'], data['cls'], occnet2gaps=data.get('occnet2gaps'))

        return est_data

    def loss(self, est_data, gt_data):
        '''
        calculate loss of est_out given gt_out.
        '''
        loss = self.mesh_reconstruction_loss(est_data, gt_data)
        total_loss = sum(loss.values())
        for key, item in loss.items():
            loss[key] = item.item()
        return {'total':total_loss, **loss}
