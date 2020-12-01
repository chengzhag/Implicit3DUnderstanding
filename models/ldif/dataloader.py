# Dataloader of LIEN+LDIF.
# author: chengzhang

import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torchvision import transforms
import torch.utils.data
from models.datasets import PIX3DLDIF
import collections
from PIL import Image
from configs.data_config import pix3d_n_classes
import numpy as np
from external.ldif.util import gaps_util, file_util
from random import random

default_collate = torch.utils.data.dataloader.default_collate


neighbors = 30

HEIGHT_PATCH = 256
WIDTH_PATCH = 256

class LDIF_Dataset(PIX3DLDIF):
    def __init__(self, config, mode):
        super(LDIF_Dataset, self).__init__(config, mode)
        mean = [0.485, 0.456, 0.406]
        std = [0.229, 0.224, 0.225]
        if self.mode=='train':
            self.data_transforms = transforms.Compose([
                transforms.Resize((280, 280)),
                transforms.RandomCrop((HEIGHT_PATCH, WIDTH_PATCH)),
                transforms.RandomHorizontalFlip(),
                transforms.ToTensor(),
                transforms.Normalize(mean, std)
            ])
        else:
            self.data_transforms = transforms.Compose([
                transforms.Resize((HEIGHT_PATCH, WIDTH_PATCH)),
                transforms.ToTensor(),
                transforms.Normalize(mean, std)
            ])

    def __getitem__(self, index):
        sample_info = self.split[index]
        sample = {}

        image = np.load(sample_info['img_path'])
        image = Image.fromarray(image)
        sample['img'] = self.data_transforms(image)

        cls_codes = torch.zeros(pix3d_n_classes)
        cls_codes[sample_info['class_id']] = 1
        sample['cls'] = cls_codes

        if self.config['model']['mesh_reconstruction']['method'] == 'LDIF':
            if self.mode == 'test':
                sample['mesh'] = file_util.read_mesh(sample_info['mesh_path'])
                occnet2gaps = file_util.read_txt_to_np(sample_info['occnet2gaps_path'])
                sample['occnet2gaps'] = np.reshape(occnet2gaps, [4, 4])
            else:
                near_surface_samples = gaps_util.read_pts_file(sample_info['nss_points_path'])
                p_ids = np.random.choice(near_surface_samples.shape[0],
                                         self.config['data']['near_surface_samples'],
                                         replace=False)
                near_surface_samples = near_surface_samples[p_ids, :]
                sample['near_surface_class'] = (near_surface_samples[:, 3:] > 0).astype(np.float32)
                sample['near_surface_samples'] = near_surface_samples[:, :3]

                uniform_samples = gaps_util.read_pts_file(sample_info['uniform_points_path'])
                p_ids = np.random.choice(uniform_samples.shape[0],
                                         self.config['data']['uniform_samples'],
                                         replace=False)
                uniform_samples = uniform_samples[p_ids, :]
                sample['uniform_class'] = (uniform_samples[:, 3:] > 0).astype(np.float32)
                sample['uniform_samples'] = uniform_samples[:, :3]

                sample['world2grid'], sample['grid'] = file_util.read_grd(sample_info['coarse_grid_path'])
                # from external.PIFu.lib import sample_util
                # sample_util.save_samples_truncted_prob('near_surface_samples.ply', sample['near_surface_samples'], sample['near_surface_class'])
                # sample_util.save_samples_truncted_prob('uniform_samples.ply', sample['uniform_samples'], sample['uniform_class'])
        elif self.config['model']['mesh_reconstruction']['method'] == 'DensTMNet':
            sample['sequence_id'] = sample_info['sample_id']
            sample['mesh_points'] = np.fromfile(
                sample_info['gt_3dpoints_path'], dtype=np.float).reshape(-1, 3).astype(np.float32)
            sample['densities'] = np.fromfile(
                sample_info['densities_path'], dtype=np.float).astype(np.float32)
            if self.mode == 'train':
                p_ids = np.random.choice(sample['mesh_points'].shape[0], 5000, replace=False)
                sample['mesh_points'] = sample['mesh_points'][p_ids, :]
                sample['densities'] = sample['densities'][p_ids]
        else:
            raise NotImplementedError
        sample.update(sample_info)
        return sample

def recursive_convert_to_torch(elem):
    if torch.is_tensor(elem):
        return elem
    elif type(elem).__module__ == 'numpy':
        if elem.size == 0:
            return torch.zeros(elem.shape).type(torch.DoubleTensor)
        else:
            return torch.from_numpy(elem)
    elif isinstance(elem, int):
        return torch.LongTensor([elem])
    elif isinstance(elem, float):
        return torch.DoubleTensor([elem])
    elif isinstance(elem, collections.Mapping):
        return {key: recursive_convert_to_torch(elem[key]) for key in elem}
    elif isinstance(elem, collections.Sequence):
        return [recursive_convert_to_torch(samples) for samples in elem]
    else:
        return elem

def collate_fn(batch):
    """
    Data collater.

    Assumes each instance is a dict.
    Applies different collation rules for each field.
    Args:
        batches: List of loaded elements via Dataset.__getitem__
    """
    collated_batch = {}
    # iterate over keys
    for key in batch[0]:
        try:
            collated_batch[key] = default_collate([elem[key] for elem in batch])
        except TypeError:
            collated_batch[key] = [elem[key] for elem in batch]

    return collated_batch

def LDIF_dataloader(config, mode='train'):
    dataloader = DataLoader(dataset=LDIF_Dataset(config, mode),
                            num_workers=config['device']['num_workers'],
                            batch_size=config[mode]['batch_size'],
                            shuffle=(mode == 'train'),
                            collate_fn=collate_fn)
    return dataloader
