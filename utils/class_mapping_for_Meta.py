from utils.sunrgbd_utils import get_NYU37_class_id
import scipy.io as scio
from tqdm import tqdm
from configs.data_config import NYU40CLASSES
import numpy as np

meta_path = 'external/cooperative_scene_parsing/evaluation/SUNRGBDtoolbox/Metadata/SUNRGBDMeta3DBB_v2.mat'
save_path = 'external/cooperative_scene_parsing/evaluation/vis/SUNRGBDMeta37.mat'
meta_data = scio.loadmat(meta_path)['SUNRGBDMeta']

for i in tqdm(range(len(meta_data[0]))):
    if len(meta_data[0, i]['groundtruth3DBB']) > 0:
        for j in range(len(meta_data[0, i]['groundtruth3DBB'][0])):
            classname = meta_data[0, i]['groundtruth3DBB']['classname'][0, j][0]
            meta_data[0, i]['groundtruth3DBB']['classname'][0, j] = \
                np.array([NYU40CLASSES[get_NYU37_class_id([classname])[1]]])

scio.savemat(save_path, {'SUNRGBDMeta': meta_data}, do_compression=True)