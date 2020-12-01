# Testing functions.
# author: ynie
# date: April, 2020
from net_utils.utils import LossRecorder, ETA
from time import time
import numpy as np
import torch
import wandb


def test_func(cfg, tester, test_loader):
    '''
    test function.
    :param cfg: configuration file
    :param tester: specific tester for networks
    :param test_loader: dataloader for testing
    :return:
    '''
    batch_size = cfg.config[cfg.config['mode']]['batch_size']
    loss_recorder = LossRecorder(batch_size)
    cfg.log_string('-'*100)
    eta_calc = ETA(smooth=0.99, ignore_first=True)
    for iter, data in enumerate(test_loader):
        loss = tester.test_step(data)

        # visualize intermediate results.
        tester.visualize_step(0, cfg.config['mode'], iter, data)

        loss_recorder.update_loss(loss, data.get('class_name', None))

        eta = eta_calc(len(test_loader) - iter - 1)
        if ((iter + 1) % cfg.config['log']['print_step']) == 0:
            pretty_loss = []
            for k, v in loss.items():
                if isinstance(v, list):
                    pretty_loss.append(str(k) + ': [' + ', '.join([f'{i:.3f}' for i in v]) + ']')
                else:
                    pretty_loss.append(f"{k}: {v:.3f}")
            pretty_loss = '{' + ', '.join(pretty_loss) + '}'
            cfg.log_string('Process: Phase: %s. Epoch %d: %d/%d. ETA: %s. Current loss: %s.'
                           % (cfg.config['mode'], 0, iter + 1, len(test_loader), eta, pretty_loss))
            wandb.summary['ETA'] = str(eta)
            for key, test_loss in loss_recorder.loss_recorder.items():
                cfg.log_string('Test loss (%s): %f' % (key, test_loss.avg))

    return loss_recorder.loss_recorder

def test(cfg, tester, test_loader):
    '''
    train epochs for network
    :param cfg: configuration file
    :param tester: specific tester for networks
    :param test_loader: dataloader for testing
    :return:
    '''
    cfg.log_string('-' * 100)
    # set mode
    tester.net.train(cfg.config['mode'] == 'train')
    start = time()
    with torch.no_grad():
        test_loss_recoder = test_func(cfg, tester, test_loader)
    cfg.log_string('Test time elapsed: (%f).' % (time()-start))
    table = None
    for key, test_loss in test_loss_recoder.items():
        cfg.log_string('Test loss (%s): %f' % (key, test_loss.avg))
        wandb.summary.update({f"{key}_avg": test_loss.avg})
        wandb.summary.update({f"{key}_hist": wandb.Histogram(test_loss.val)})
        if len(test_loss.cls) > 0:
            if table is None:
                table = wandb.Table(columns=['metric'] + [k for k in test_loss.cls.keys()] + ['mean'])
            cfg.log_string({k: v.avg for k, v in test_loss.cls.items()})
            table.add_data(key, *[f"{v.avg:.2f}" for v in test_loss.cls.values()], f"{test_loss.avg:.5f}")
    if table is not None:
        wandb.summary['metrics_table'] = table
