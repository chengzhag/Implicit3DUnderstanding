from utils.sunrgbd_config import SUNRGBD_CONFIG
from configs.data_config import NYU40CLASSES
from utils.sunrgbd_utils import readsunrgbdframe, process_sunrgbd_frame, check_bdb
from utils.vis_tools_sunrgbd import Scene3D_SUNRGBD
import os
import numpy as np
import pickle
from configs.data_config import Config
from utils.sunrgbd_utils import proj_from_point_to_2d, get_corners_of_bb3d_no_index
from net_utils.libs import get_iou
from libs.tools import camera_cls_reg_sunrgbd, layout_size_avg_residual, ori_cls_reg, obj_size_avg_residual, bin_cls_reg, list_of_dict_to_dict_of_list
import json
from multiprocessing import Pool
import shutil
from pathlib import Path
from utils.sunrgbd_utils import get_NYU37_class_id
from configs.data_config import NYU40CLASSES


if __name__ == '__main__':
    config = SUNRGBD_CONFIG()

    # visualize SUNRGBD samples
    np.set_printoptions(suppress=True)
    for i, id in enumerate((276, 724, 765, 845, 1109, 1140, 1149, 2435)):
        input_folder = f"demo/inputs/{i+1}/"
        os.makedirs(input_folder, exist_ok=True)
        with open(os.path.join(config.clean_data_root, 'data_all', str(id) + '.pickle'), 'rb') as f:
            img_info = pickle.load(f, encoding='latin1')
        # Path(input_folder + "cam_K.txt").touch()
        np.savetxt(input_folder + "cam_K.txt", img_info['K'], fmt='%.4f')
        with open(input_folder + "detections.json", 'w') as f:
            bdb2d = [{'bbox': [b['x1'], b['y1'], b['x2'], b['y2']],
                      'class': NYU40CLASSES[get_NYU37_class_id([b['classname']])[1]]} for b in img_info['bdb2d']]
            json.dump(bdb2d, f)
        imgrgb_path = img_info['imgrgb_path'].replace('/home/siyuan/Documents/Dataset/SUNRGBD_ALL', config.data_root)
        shutil.copy(imgrgb_path, input_folder + "img.jpg")
