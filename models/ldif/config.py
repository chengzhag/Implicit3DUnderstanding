# Configure trainer and tester
# author: chengzhang

from .training import Trainer
from .testing import Tester
from .dataloader import LDIF_dataloader

def get_trainer(cfg, net, optimizer, device=None):
    return Trainer(cfg=cfg, net=net, optimizer=optimizer, device=device)

def get_tester(cfg, net, device=None):
    return Tester(cfg=cfg, net=net, device=device)

def get_dataloader(config, mode):
    return LDIF_dataloader(config=config, mode=mode)