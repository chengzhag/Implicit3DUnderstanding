# Main file for training and testing
# author: ynie
# date: Feb, 2020

import argparse
from configs.config_utils import CONFIG
import os
import train, test,demo


def parse_args():
    '''PARAMETERS'''
    parser = argparse.ArgumentParser('Total 3D Understanding.')
    parser.add_argument('config', type=str, default='configs/total3d_mgnet.yaml',
                        help='configure file for training or testing.')
    parser.add_argument('--mode', type=str, default='qtrain', help='train, test, demo or qtrain, qtest')
    parser.add_argument('--demo_path', type=str, default='demo/inputs/1', help='Please specify the demo path.')
    parser.add_argument('--name', type=str, default=None, help='wandb exp name.')
    parser.add_argument('--sweep', action='store_true')
    return parser

if __name__ == '__main__':
    parser = parse_args()
    cfg = CONFIG(parser)

    '''Configuration'''
    cfg.log_string('Loading configurations.')
    cfg.log_string(cfg.config)

    '''Run'''
    if cfg.config['mode'] == 'train':
        try:
            train.run(cfg)
        except KeyboardInterrupt:
            pass
        except:
            raise
        cfg.update_config(mode='test', resume=True, weight=os.path.join(cfg.save_path, 'model_best.pth'))
    if cfg.config['mode'] == 'test':
        test.run(cfg)
    if cfg.config['mode'] == 'demo':
        demo.run(cfg)

