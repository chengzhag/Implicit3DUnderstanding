# Holistic 3D Scene Understanding from a Single Image with Implicit Representation

Cheng Zhang, Zhaopeng Cui, Yinda Zhang, Shuaicheng Liu, Bing Zeng, Marc Pollefeys

<img src="demo/inputs/1/img.jpg" alt="img.jpg" width="20%" /> <img src="demo/outputs/1/3dbbox.png" alt="3dbbox.png" width="20%" /> <img src="demo/outputs/1/recon.png" alt="recon.png" width="20%" /> <br>
<img src="demo/inputs/2/img.jpg" alt="img.jpg" width="20%" /> <img src="demo/outputs/2/3dbbox.png" alt="3dbbox.png" width="20%" /> <img src="demo/outputs/2/recon.png" alt="recon.png" width="20%" /> <br>
<img src="demo/inputs/3/img.jpg" alt="img.jpg" width="20%" /> <img src="demo/outputs/3/3dbbox.png" alt="3dbbox.png" width="20%" /> <img src="demo/outputs/3/recon.png" alt="recon.png" width="20%" />

![pipeline](figures/pipeline.png)

## Install

```
conda env create -f environment.yml
conda activate Im3D
python project.py build
```


## Demo
1.Download the pretrained checkpoint and unzip it into
```
out/total3d/20110611514267/
```
2.Change current directory to ```Implicit3DUnderstanding/``` and run the demo, which will generate 3D detection result and rendered scene mesh to ```demo/output/1/```
```
CUDA_VISIBLE_DEVICES=0 python main.py out/total3d/20110611514267/out_config.yaml --mode demo --demo_path demo/inputs/1
```
3.In case you want to run it off screen (for example, with SSH)
```
sudo apt install xvfb
CUDA_VISIBLE_DEVICES=0 xvfb-run -a -s "-screen 0 800x600x24" python main.py out/total3d/20110611514267/out_config.yaml --mode demo --demo_path demo/inputs/1
```
4.If you want to run it interactively, change the last line of demo.py
```
scene_box.draw3D(if_save=True, save_path = '%s/recon.png' % (save_path))
```
to
```
scene_box.draw3D(if_save=False, save_path = '%s/recon.png' % (save_path))
```


## Data preparation
We follow [Total3DUnderstanding](https://github.com/yinyunie/Total3DUnderstanding) to use [SUN-RGBD](https://rgbd.cs.princeton.edu/) to train our Scene Graph Convolutional Network (SGCN), and use [Pix3D](http://pix3d.csail.mit.edu/) to train our Local Implicit Embedding Network
(LIEN) with [Local Deep Implicit Functions](https://github.com/google/ldif) (LDIF) decoder.

##### Preprocess SUN-RGBD data

Please follow [Total3DUnderstanding](https://github.com/yinyunie/Total3DUnderstanding) to directly download or generate the processed train/test data.

In case you prefer processing by yourself,
according to [issue #6](https://github.com/yinyunie/Total3DUnderstanding/issues/6) of [Total3DUnderstanding](https://github.com/yinyunie/Total3DUnderstanding),
there are a few typos in json files of SUNRGBD dataset, which is mostly solve by the json loader.
However, one typo still needs to be fixed by hand.
Please find ```{"name":""propulsion"tool"}``` in ```data/sunrgbd/Dataset/SUNRGBD/kv2/kinect2data/002922_2014-06-26_15-43-16_094959634447_rgbf000089-resize/annotation2Dfinal/index.json``` and remove ```""propulsion```.

##### Preprocess Pix3D data
We use a different data process pipeline with [Total3DUnderstanding](https://github.com/yinyunie/Total3DUnderstanding). Please follow these steps to generate the train/test data:

1. Download the [Pix3D dataset](http://pix3d.csail.mit.edu/) to 
```
data/pix3d/metadata
```
2. Run below to generate the train/test data into 'data/pix3d/ldif'
```
python utils/preprocess_pix3d4ldif.py
```


## Training and Testing
We use [wandb](https://www.wandb.com/) for logging and visualization.
You can register a wandb account and login before training by ```wandb login```.
In case you don't need to visualize the training process, you can put ```WANDB_MODE=dryrun``` before the commands bellow.

Thanks to the well-structured code of [Total3DUnderstanding](https://github.com/yinyunie/Total3DUnderstanding), we use the same method to manage parameters of each experiment with configuration files (```configs/****.yaml```).
We first follow [Total3DUnderstanding](https://github.com/yinyunie/Total3DUnderstanding) to pretrain each individual module, then jointly finetune the full model with additional physical violation loss.

##### Pretraining
We use the [pretrained checkpoint](https://livebournemouthac-my.sharepoint.com/:u:/g/personal/ynie_bournemouth_ac_uk/EcA66Nb1aI1KitzX7avbE10BiHGzovf3rqQebeJHmFB4QA?e=4hE8zv) of [Total3DUnderstanding](https://github.com/yinyunie/Total3DUnderstanding) to load weights for ODN.
Please download and rename the checkpoint to ```out/pretrained_models/total3d/model_best.pth```.
Other modules can be trained then tested with the following steps:

1. Train LEN by:
```
python main.py configs/layout_estimation.yaml
```
The pretrained checkpoint can be found at 
```
out/layout_estimation/[start_time]/model_best.pth
```
2. Train LIEN + LDIF by:
```
python main.py configs/ldif.yaml
```
The pretrained checkpoint can be found at
```
out/ldif/[start_time]/model_best.pth
```
The training process is followed with a quick test without scene mesh generated. In case you want to generate the object mesh and evaluate the Chamfer distance during testing:
```
python main.py configs/ldif.yaml --mode train
```
3. Replace the checkpoint directories of LEN and LIEN in ```configs/total3d_ldif_gcnn.yaml``` with the checkpoints trained above, then train SGCN by:
```
python main.py configs/total3d_ldif_gcnn.yaml
```
The pretrained checkpoint can be found at
```
out/total3d/[start_time]/model_best.pth
```

##### Joint finetune

Replace the checkpoint directory in ```configs/total3d_ldif_gcnn_joint.yaml``` with the one trained in the last step above, then train the full model by:
```
python main.py configs/total3d_ldif_gcnn_joint.yaml
```
The trained model can be found at
```
out/total3d/[start_time]/model_best.pth
```
The training process is followed with a quick test without scene mesh generated. In case you want to generate the scene mesh during testing (which will cost a day on 1080ti due to the unoptimized interface of LDIF CUDA kernel):
```
python main.py configs/total3d_ldif_gcnn_joint.yaml --mode train
```
The testing resaults can be found at
```
out/total3d/[start_time]/visualization
```

##### Testing

1. The training process above already include a testing process. In case you want to test LIEN+LDIF or full model by yourself:
```
python main.py out/[ldif/total3d]/[start_time]/model_best.pth --mode test
```
The results will be saved to evaluation metrics will be logged to wandb as run summary.

2. Visualize the i-th 3D scene interacively by
```
python utils/visualize.py --result_path out/total3d/[start_time]/visualization --sequence_id [i]
```
or save the 3D detection result and rendered scene mesh by
```
python utils/visualize.py --result_path out/total3d/[start_time]/visualization --sequence_id [i] --save_path []
```
In case you do not have a screen:
```
python utils/visualize.py --result_path out/total3d/[start_time]/visualization --sequence_id [i] --save_path [] --offscreen
```

##### About the testing speed

When using marching cube to reconstruct mesh with LDIF, the CUDA kernel is used to optimize the speed.
However, we ported the file-based interface of the kernel in [Tensorflow implementation](https://github.com/google/ldif) which is bottlenecked by the file IO.
This is the main reason why our method takes so much time in object and scene mesh reconstruction.
To improve the speed, you can lower the parameter ```data.marching_cube_resolution``` in the configuration file.

## Citation

We thank the following great works:
- [Total3DUnderstanding](https://github.com/yinyunie/Total3DUnderstanding) for their well-structured code. We construct our network based on their well-structured code.
- [Coop](https://github.com/thusiyuan/cooperative_scene_parsing) for their dataset. We used their processed dataset with 2D detector prediction.
- [LDIF](https://github.com/google/ldif) for their novel representation method. We ported their LDIF decoder from Tensorflow to PyTorch.
- [Graph R-CNN](https://github.com/jwyang/graph-rcnn.pytorch/blob/master/README.md) for their scene graph design. We adopted their GCN implemention to construct our SGCN.

If you find them helpful, please cite:
```
@InProceedings{Nie_2020_CVPR,
author = {Nie, Yinyu and Han, Xiaoguang and Guo, Shihui and Zheng, Yujian and Chang, Jian and Zhang, Jian Jun},
title = {Total3DUnderstanding: Joint Layout, Object Pose and Mesh Reconstruction for Indoor Scenes From a Single Image},
booktitle = {IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)},
month = {June},
year = {2020}
}
@inproceedings{huang2018cooperative,
  title={Cooperative Holistic Scene Understanding: Unifying 3D Object, Layout, and Camera Pose Estimation},
  author={Huang, Siyuan and Qi, Siyuan and Xiao, Yinxue and Zhu, Yixin and Wu, Ying Nian and Zhu, Song-Chun},
  booktitle={Advances in Neural Information Processing Systems},
  pages={206--217},
  year={2018}
}	
@inproceedings{genova2020local,
    title={Local Deep Implicit Functions for 3D Shape},
    author={Genova, Kyle and Cole, Forrester and Sud, Avneesh and Sarna, Aaron and Funkhouser, Thomas},
    booktitle={Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition},
    pages={4857--4866},
    year={2020}
}
@inproceedings{yang2018graph,
    title={Graph r-cnn for scene graph generation},
    author={Yang, Jianwei and Lu, Jiasen and Lee, Stefan and Batra, Dhruv and Parikh, Devi},
    booktitle={Proceedings of the European Conference on Computer Vision (ECCV)},
    pages={670--685},
    year={2018}
}
```



