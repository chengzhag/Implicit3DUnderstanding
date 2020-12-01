# Utility functions in training and testing
# author: ynie
# date: Feb, 2020

import os
import yaml
import logging
from datetime import datetime
from configs.data_config import Config as Data_Config
from net_utils.libs import to_dict_tensor
import argparse
from functools import reduce
from operator import getitem
import ast


def update_recursive(dict1, dict2):
    ''' Update two config dictionaries recursively.

    Args:
        dict1 (dict): first dictionary to be updated
        dict2 (dict): second dictionary which entries should be used

    '''
    for k, v in dict2.items():
        if k not in dict1:
            dict1[k] = dict()
        if isinstance(v, dict):
            update_recursive(dict1[k], v)
        else:
            dict1[k] = v

def set_by_path(tree, keys, value):
    """Set a value in a nested object in tree by sequence of keys."""
    keys = keys.split('.')
    get_by_path(tree, keys[:-1])[keys[-1]] = value

def get_by_path(tree, keys):
    """Access a nested object in tree by sequence of keys."""
    return reduce(getitem, keys, tree)

class CONFIG(object):
    '''
    Stores all configures
    '''
    def __init__(self, parser):
        '''
        Loads config file
        :param path (str): path to config file
        :return:
        '''
        args, modified_args = parser.parse_known_args()
        input = args.config
        mode = args.mode
        if mode in ['qtrain', 'qtest']:
            mode = mode[1:]
            full = False
        else:
            full = True

        # read config file
        self.config = self.read_to_dict(input)
        # for compatibility with newly added val mode
        if 'val' not in self.config.keys():
            self.config['val'] = self.config['test'].copy()

        # modify config from modified_args
        parser = argparse.ArgumentParser('Modified config.')
        modified_args = ' '.join(modified_args).replace('=', ' ').split(' ') if modified_args else []
        modified_keys = [modified_args[i * 2] for i in range(len(modified_args) // 2)]

        def add_argument(d, super_name=None):
            for k, v in d.items():
                n = f'{super_name}.{k}' if super_name else k
                if isinstance(v, dict):
                    add_argument(v, n)
                else:
                    flag = f'--{n}'
                    if flag in modified_keys:
                        if type(v) == list:
                            type_v = type(v[0])
                            nargs = '+'
                        else:
                            type_v = type(v)
                            nargs = None
                        if type_v is bool:
                            type_v = ast.literal_eval
                        parser.add_argument(flag, type=type_v, nargs=nargs)

        add_argument(self.config)
        modified_args = vars(parser.parse_args(modified_args))
        if mode in ['test', 'demo']:
            log_path = modified_args.get('log.path', self.config['log']['path'])
            weight = os.path.join(log_path, 'model_best.pth')
            if os.path.exists(weight):
                self.config['weight'] = weight
                modified_args['resume'] = True
            else:
                modified_args['resume'] = False
        resume = modified_args.get('resume', self.config['resume'])
        updated_args = []
        for k, v in modified_args.items():
            if resume and k in ['resume']:
                print(f'Warning: when resuming, config change "{k}: {v}" will not be saved')
                continue
            updated_args.append(k)
            set_by_path(self.config, k, v)

        # initialize save_path, logger, vis_path
        if resume:
            self._save_path = self.config['log']['path']
            self._logger = self.load_logger()
        else:
            # update save_path to config file
            self._save_path = os.path.join(self.config['log']['path'], datetime.now().strftime('%y%m%d%H%M%S%f')[:-4])
            self._logger = self.load_logger()
            self.update_config(log={'path': self._save_path})

            # update visualization path
            vis_path = os.path.join(self._save_path, 'visualization')
            if not os.path.exists(vis_path):
                os.mkdir(vis_path)
            self.update_config(log={'vis_path': vis_path})

        # write new config file or update existing config file
        self.write_config()

        # update configs not saved to config file
        for k, v in modified_args.items():
            if k not in updated_args:
                set_by_path(self.config, k, v)
        args_dict = args.__dict__
        args_dict['mode'] = mode
        self.update_config(args_dict)
        self.config['full'] = full

        # initiate environments
        if "CUDA_VISIBLE_DEVICES" not in os.environ:
            os.environ["CUDA_VISIBLE_DEVICES"] = self.config['device']['gpu_ids']
        else:
            self.config['device']['gpu_ids'] = os.environ["CUDA_VISIBLE_DEVICES"]
        from net_utils.utils import initiate_environment
        initiate_environment(self.config)

    @property
    def logger(self):
        return self._logger

    @property
    def save_path(self):
        return self._save_path

    def load_logger(self):
        # set file handler
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        logfile = os.path.join(self.save_path, 'log.txt')
        file_handler = logging.FileHandler(logfile)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        self.__file_handler = file_handler

        # configure logger
        logger = logging.getLogger('Empty')
        logger.setLevel(logging.INFO)
        logger.addHandler(file_handler)

        return logger

    def log_string(self, content):
        self._logger.info(content)
        print(content)

    def read_to_dict(self, input):
        if not input:
            return dict()
        if isinstance(input, str) and os.path.isfile(input):
            if input.endswith('yaml'):
                with open(input, 'r') as f:
                    config = yaml.load(f, Loader=yaml.FullLoader)
            else:
                ValueError('Config file should be with the format of *.yaml')
        elif isinstance(input, dict):
            config = input
        else:
            raise ValueError('Unrecognized input type (i.e. not *.yaml file nor dict).')

        return config

    def update_config(self, *args, **kwargs):
        '''
        update config and corresponding logger setting
        :param input: dict settings add to config file
        :return:
        '''
        cfg1 = dict()
        for item in args:
            cfg1.update(self.read_to_dict(item))

        cfg2 = self.read_to_dict(kwargs)

        new_cfg = {**cfg1, **cfg2}

        update_recursive(self.config, new_cfg)
        # when update config file, the corresponding logger should also be updated.
        self.__update_logger()

    def write_config(self):
        output_file = os.path.join(self._save_path, 'out_config.yaml')

        with open(output_file, 'w') as file:
            yaml.dump(self.config, file, default_flow_style = False)

    def __update_logger(self):
        # configure logger
        name = self.config['mode'] if 'mode' in self.config else self._logger.name
        logger = logging.getLogger(name)
        logger.setLevel(logging.INFO)
        logger.addHandler(self.__file_handler)
        self._logger = logger

def mount_external_config(cfg):
    if cfg.config['data']['dataset'] == 'sunrgbd':
        dataset_config = Data_Config('sunrgbd')
        bins_tensor = to_dict_tensor(dataset_config.bins, if_cuda=cfg.config['device']['use_gpu'])
        setattr(cfg, 'dataset_config', dataset_config)
        setattr(cfg, 'bins_tensor', bins_tensor)
    return cfg
