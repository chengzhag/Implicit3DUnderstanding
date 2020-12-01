# Training script
# author: ynie
# date: Feb, 2020
from models.optimizers import load_optimizer, load_scheduler
from net_utils.utils import load_device, load_model, load_trainer, load_dataloader
from net_utils.utils import CheckpointIO
from train_epoch import train
from configs.config_utils import mount_external_config
import os
import wandb


def run(cfg):
    resume = cfg.config['resume']
    name = cfg.config['name']
    id = cfg.config['log']['path'].split('/')[-1]
    if not resume and not name:
        name = id
    if cfg.config['sweep']:
        name = None
        id = None
    wandb.init(project="implicit3dunderstanding", config=cfg.config, dir=cfg.config['log']['path'],
               name=name, id=id, resume=resume)
    wandb.summary['pid'] = os.getpid()
    wandb.summary['ppid'] = os.getppid()
    # if cfg.config != dict(wandb.config):
    #     cfg.log_string('Updating configurations from wandb.')
    #     cfg.config.update(wandb.config)
    #     cfg.write_config()
    #     cfg.log_string(cfg.config)

    if resume:
        cfg.update_config(weight=os.path.join(cfg.config['log']['path'], 'model_last.pth'))

    '''Begin to run network.'''
    checkpoint = CheckpointIO(cfg)

    '''Mount external config data'''
    cfg = mount_external_config(cfg)

    '''Load save path'''
    cfg.log_string('Data save path: %s' % (cfg.save_path))

    '''Load device'''
    cfg.log_string('Loading device settings.')
    device = load_device(cfg)

    '''Load data'''
    cfg.log_string('Loading dataset.')
    train_loader = load_dataloader(cfg.config, mode='train')
    test_loader = load_dataloader(cfg.config, mode='val')

    '''Load net'''
    cfg.log_string('Loading model.')
    net = load_model(cfg, device=device)
    checkpoint.register_modules(net=net)
    cfg.log_string(net)
    wandb.watch(net, log=None)

    '''Load optimizer'''
    cfg.log_string('Loading optimizer.')
    optimizer = load_optimizer(config=cfg.config, net=net)
    checkpoint.register_modules(optimizer=optimizer)

    '''Load scheduler'''
    cfg.log_string('Loading optimizer scheduler.')
    scheduler = load_scheduler(config=cfg.config, optimizer=optimizer)
    checkpoint.register_modules(scheduler=scheduler)

    '''Check existing checkpoint (resume or finetune)'''
    checkpoint.parse_checkpoint()

    '''Load trainer'''
    cfg.log_string('Loading trainer.')
    trainer = load_trainer(cfg=cfg, net=net, optimizer=optimizer, device=device)

    '''Start to train'''
    cfg.log_string('Start to train.')
    num_params = sum(p.numel() for p in net.parameters())
    wandb.summary['num_params'] = num_params
    cfg.log_string('Total number of parameters in {0:s}: {1:d}.'.format(cfg.config['method'], num_params))

    train(cfg=cfg, trainer=trainer, scheduler=scheduler, checkpoint=checkpoint, train_loader=train_loader, val_loader=test_loader)

    cfg.log_string('Training finished.')