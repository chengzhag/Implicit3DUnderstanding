# Base data of networks
# author: ynie
# date: Feb, 2020
import os
from torch.utils.data import Dataset
import json
from configs.pix3d_config import Config
from tqdm import tqdm


class SUNRGBD(Dataset):
    def __init__(self, config, mode):
        '''
        initiate SUNRGBD dataset for data loading
        :param config: config file
        :param mode: train/val/test mode
        '''
        self.config = config
        if mode == 'val':
            mode = 'test'
        self.mode = mode
        split_file = os.path.join(config['data']['split'], mode + '.json')
        with open(split_file) as file:
            split = json.load(file)
        self.split = []
        skipped = 0
        for s in tqdm(split):
            if os.path.exists(s):
                self.split.append(s)
            else:
                skipped += 1
        print(f'{skipped}/{len(split)} missing samples')

    def __len__(self):
        return len(self.split)


class PIX3D(Dataset):
    def __init__(self, config, mode):
        '''
        initiate PIX3D dataset for data loading
        :param config: config file
        :param mode: train/val/test mode
        '''
        self.config = config
        if mode == 'val':
            mode = 'test'
        self.mode = mode
        split_file = os.path.join(config['data']['split'], mode + '.json')
        with open(split_file) as file:
            self.split = json.load(file)

    def __len__(self):
        return len(self.split)


class PIX3DLDIF(Dataset):
    def __init__(self, config, mode):
        '''
        initiate PIX3DLDIF dataset for data loading
        :param config: config file
        :param mode: train/val/test mode
        '''
        self.config = config
        self.mode = mode
        if mode == 'val':
            mode = 'test'
        split_file = os.path.join(config['data']['split'], mode + '.json')
        with open(split_file) as file:
            split = json.load(file)

        config_data = Config('pix3d')
        with open(config_data.metadata_file, 'r') as file:
            metadatas = json.load(file)
        ids = [int(os.path.basename(file).split('.')[0]) for file in split if 'flipped' not in file]
        sample_info = []
        skipped = 0
        for id in tqdm(ids):
            metadata = metadatas[id]
            info = {}

            rel_img = metadata['img'].replace('img/', '').split('/')
            model_folder = '.'.join(os.path.splitext(metadata['model'])[0].split('/')[-2:])
            rel_folder = os.path.join(config_data.root_path, 'ldif', *rel_img[:-1], model_folder)
            img_name = os.path.splitext(rel_img[-1])[0]
            info['img_path'] = os.path.join(rel_folder, img_name + '.npy')
            info['nss_points_path'] = os.path.join(rel_folder, 'nss_points.sdf')
            info['uniform_points_path'] = os.path.join(rel_folder, 'uniform_points.sdf')
            info['coarse_grid_path'] = os.path.join(rel_folder, 'coarse_grid.grd')
            info['occnet2gaps_path'] = os.path.join(rel_folder, 'orig_to_gaps.txt')
            watertight = self.config['data'].get('watertight',
                                                 self.config['model']['mesh_reconstruction']['method'] == 'LDIF')
            info['mesh_path'] = os.path.join(rel_folder, 'mesh_orig.ply' if watertight else 'mesh_normalized.ply')
            ext_mgnet = 'mgn' if watertight else 'org'
            info['gt_3dpoints_path'] = os.path.join(rel_folder, f'gt_3dpoints.{ext_mgnet}')
            info['densities_path'] = os.path.join(rel_folder, f'densities.{ext_mgnet}')
            if not all([os.path.exists(path) for path in info.values()]) :
                skipped += 1
                continue

            info['sample_id'] = id
            info['class_id'] = config_data.classnames.index(metadata['category'])
            info['class_name'] = metadata['category']
            sample_info.append(info)
        print(f'{skipped}/{len(ids)} missing samples')
        self.split = sample_info

    def __len__(self):
        return len(self.split)