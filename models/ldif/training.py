# Trainer for LIEN+LDIF.
# author: chengzhang

from models.training import BaseTrainer
import torch


class Trainer(BaseTrainer):
    '''
    Trainer object for total3d.
    '''
    def eval_step(self, data):
        '''
        performs a step in evaluation
        :param data (dict): data dictionary
        :return:
        '''
        loss = self.compute_loss(data)
        loss['total'] = loss['total'].item()
        return loss

    def visualize_step(self, epoch, phase, iter, data):
        ''' Performs a visualization step.
        '''
        pass

    def to_device(self, data):
        device = self.device
        ndata = {}
        for k, v in data.items():
            if type(v) is torch.Tensor and v.dtype is torch.float32:
                ndata[k] = v.to(device)
            else:
                ndata[k] = v
        return ndata

    def compute_loss(self, data):
        '''
        compute the overall loss.
        :param data (dict): data dictionary
        :return:
        '''
        '''load input and ground-truth data'''
        data = self.to_device(data)

        '''network forwarding'''
        est_data = self.net(data)

        '''computer losses'''
        loss = self.net.loss(est_data, data)
        return loss
