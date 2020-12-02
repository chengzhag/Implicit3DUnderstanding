clc;
clear all;

old_SUNRGBD_path = '/home/siyuan/Documents/Dataset/SUNRGBD_ALL/SUNRGBD';
SUNRGBD_path = '/home/zhangcheng/datasets/SUNRGBD/SUNRGBD';
toolboxpath = '../../evaluation/SUNRGBDtoolbox';
addpath(genpath(toolboxpath));
load(fullfile(toolboxpath,'/traintestSUNRGBD/allsplit.mat'));
load(fullfile(toolboxpath,'/Metadata/SUNRGBDMeta.mat'));
new_3d = load(fullfile(toolboxpath, '/Metadata/SUNRGBDMeta3DBB_v2.mat'));
cls_mapping = readtable('../../../../data/sunrgbd/class_mapping_from_toolbox.csv');
cls = {'void', 'wall', 'floor', 'cabinet', 'bed', 'chair', ...
    'sofa', 'table', 'door', 'window', 'bookshelf', ...
    'picture', 'counter', 'blinds', 'desk', 'shelves', ...
    'curtain', 'dresser', 'pillow', 'mirror', 'floor_mat', ...
    'clothes', 'ceiling', 'books', 'refridgerator', 'television', ...
    'paper', 'towel', 'shower_curtain', 'box', 'whiteboard', ...
    'person', 'night_stand', 'toilet', 'sink', 'lamp', ...
    'bathtub', 'bag', 'otherstructure', 'otherfurniture', 'otherprop'};
disp_cls = {'cabinet', 'bed', 'chair', 'sofa', 'table', 'door', ...
    'window', 'bookshelf', 'picture', 'counter', 'blinds', ...
    'desk', 'shelves', 'curtain', 'dresser', 'pillow', ...
    'mirror', 'clothes', 'books', 'refridgerator', 'television', ...
    'paper', 'towel', 'shower_curtain', 'box', 'whiteboard', ...
    'person', 'night_stand', 'toilet', 'sink', 'lamp', 'bathtub', ...
    'bag'};
thresh_iou = [0.01,0.05:0.05:0.6,0.8:0.2:1];
result_path = '/home/zhangcheng/projects/sdf_seg/src/Total3DUnderstanding/out/total3d/20111223070510/visualization';
vis = false;

%% evaluation
iou_precision =[];
iou_recall =[];
label_recall =[];
label_precision = [];
iouHolistic_all = []; 
A = [];
testset = zeros(5050, 1);
index = 1;
for imageId = 1:5050
    data = SUNRGBDMeta(imageId);
    if ~exist(fullfile(result_path, int2str(imageId), 'co_bdb_3d.mat'), 'file');
        continue;
    end
    if size(data.gtCorner3D, 2) == 0
        continue;
    end
    sequenceMapping{imageId} = index;
    index = index + 1;
end

for imageId = 1:5050
    data = SUNRGBDMeta(imageId);
    if ~exist(fullfile(result_path, int2str(imageId), 'co_bdb_3d.mat'), 'file');
        continue;
    end
    if size(data.gtCorner3D, 2) == 0
        continue;
    end
    testset(imageId) = 1;
    disp(imageId);
    bdb_3d = load(fullfile(result_path, int2str(imageId), 'co_bdb_3d.mat'));
    bdb = bdb_3d.bdb;
    layout = load(fullfile(result_path, int2str(imageId), 'co_layout.mat'));
    predictRoom = layout.layout;
    gt_bdb = new_3d.SUNRGBDMeta(imageId).groundtruth3DBB;
    % evaluate
    for i = 1:size(bdb, 2)
    	bdb_temp = bdb{i};
        bdb_temp.sceneid = sequenceMapping{imageId};
        bdb_temp.confidence = 1;
        predictedBbs(i) = bdb_temp;
        A = [A,predictedBbs(i)];
    end
end
%%
% segment

bdb_all = A';
num_class = length(cls);
bdb_total{num_class} = [];
for i = 1:size(bdb_all, 1)
    classid = bdb_all(i, 1).classid+1;
    sceneid = bdb_all(i, 1).sceneid;
    bdb_info = bdb_all(i, 1);
    if isempty(bdb_total{classid})
        bdb_total{classid}.allBb3dtight = bdb_info;
        bdb_total{classid}.allTestImgIds = sceneid;
    else
        bdb_total{classid}.allBb3dtight = [bdb_total{classid}.allBb3dtight; bdb_info];
        bdb_total{classid}.allTestImgIds = [bdb_total{classid}.allTestImgIds; sceneid];
    end
end
tempfolder = tempname;
mkdir(tempfolder);
for i = 1:num_class
    if ~isempty(bdb_total{i})
        allBb3dtight = bdb_total{i}.allBb3dtight;
        allTestImgIds = bdb_total{i}.allTestImgIds;
        save([tempfolder '/' cls{i} '.mat'], 'allBb3dtight', 'allTestImgIds');
    end
end

%%
% compute 
for i = 1:length(disp_cls)
className = disp_cls{i};
obj_path = fullfile(tempfolder, [className '.mat']);
if ~exist(obj_path, 'file')
    continue
end
load(obj_path,'allTestImgIds','allBb3dtight')
split = load(fullfile(toolboxpath,'/traintestSUNRGBD/allsplit.mat'));
testset_path = split.alltest;
testset_path = testset_path(testset==1);
[groundTruthBbs,all_sequenceName] = benchmark_groundtruth(className,fullfile(toolboxpath,'Metadata/'),testset_path,cls_mapping);
if size(groundTruthBbs, 1) == 0
    continue
end
[apScore,precision,recall,isTp,isFp,isMissed,gtAssignment,maxOverlaps] = computePRCurve3D(className,allBb3dtight,allTestImgIds,groundTruthBbs,zeros(length(groundTruthBbs),1));
result_all = struct('apScore',apScore,'precision',precision,'recall',recall,'isTp',isTp,'isFp',isFp,'isMissed',isMissed,'gtAssignment',gtAssignment);
disp(className);
disp(apScore);
score{i} = apScore;
%% 
% figure,
% plot(recall,precision)
% title(className)
end
disp(mean(cell2mat(score(:))));
%%
rmdir(tempfolder, 's');
