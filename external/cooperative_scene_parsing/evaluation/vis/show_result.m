%% show our result
clc;
close all;
clear all;
addpath(genpath('.'))
SUNRGBDMeta = load('SUNRGBDMetaUS.mat', 'SUNRGBDMeta');
SUNRGBDMeta = SUNRGBDMeta.SUNRGBDMeta;
old_SUNRGBD_path = '/home/siyuan/Documents/Dataset/SUNRGBD_ALL/SUNRGBD';
SUNRGBDMetav2 = load('SUNRGBDMeta37.mat', 'SUNRGBDMeta');
SUNRGBDMetav2 = SUNRGBDMetav2.SUNRGBDMeta;
SUNRGBD_path = 'Z:/datasets/SUNRGBD/SUNRGBD';
cls = {'void', 'wall', 'floor', 'cabinet', 'bed', 'chair', ...
    'sofa', 'table', 'door', 'window', 'bookshelf', ...
    'picture', 'counter', 'blinds', 'desk', 'shelves', ...
    'curtain', 'dresser', 'pillow', 'mirror', 'floor_mat', ...
    'clothes', 'ceiling', 'books', 'refridgerator', 'television', ...
    'paper', 'towel', 'shower_curtain', 'box', 'whiteboard', ...
    'person', 'night_stand', 'toilet', 'sink', 'lamp', ...
    'bathtub', 'bag', 'otherstructure', 'otherfurniture', 'otherprop'};
load('nyu_colorbox.mat');
nyu_colorbox = nyu_colorbox * 8/9;

save_root = 'Z:\projects\sdf_seg\src\paper\scnrecon\';
paths = {
    {fullfile(save_root, 'ours_coop'), ...
    'Z:/projects/sdf_seg/src/Total3DUnderstanding/out/total3d/20111223070510/visualization'}, ...
    {fullfile(save_root, 'gt_coop')}, ...
    {fullfile(save_root, 'Total3D_coop'), ...
    'Z:/projects/sdf_seg/src/Total3DUnderstanding/out/total3d/20110212532979/visualization'}
    };
% search_folder = 'Z:\projects\sdf_seg\src\paper\scnrecon\ours_20110611514267';
vis_pc = false;
views3d = {'oblique'}; % top or oblique
dosave = true;
legend = table('Size', [0, 2], 'VariableNames', {'class', 'color'}, 'VariableTypes', {'string', 'string'});

if exist('search_folder', 'var')
    image_files = dir(search_folder);
    image_files = {image_files(3:end).name};
    imageIds = cellfun(@(x)sscanf(x, '%d'), image_files);
    imageIds = unique(imageIds);
else
    imageIds = [276, 765, 1149]; % set id of scene you want to visualize
end

%% change path
for imageId = imageIds
    try
        viewsCam = {};
        for ipath = 1:length(paths)
            path = paths{ipath};
            save_path = path{1};
            gt = true;
            if length(path) == 2
                result_path = path{2};
                gt = false;
            end
            if dosave
                if ~exist(save_path, 'dir')
                    mkdir(save_path);
                end
            end

            bdb_2d = {};
            disp(imageId);
            data = SUNRGBDMeta(imageId);
            data.depthpath = strrep(data.depthpath, old_SUNRGBD_path, SUNRGBD_path);
            data.rgbpath = strrep(data.rgbpath, old_SUNRGBD_path, SUNRGBD_path);
            [rgb,points3d,depthInpaint,imsize]=read3dPoints(data);
            if ~gt
                bdb_result_path = fullfile(result_path, int2str(imageId), 'co_bdb_3d.mat');
                layout_path = fullfile(result_path, int2str(imageId), 'co_layout.mat');
                if ~exist(bdb_result_path, 'file')
                    continue;
                end
                r_path = fullfile(result_path, int2str(imageId), 'co_r_ex.mat');
                load(r_path);
                load(bdb_result_path);
                bdb_3d = {};
                for i = 1:size(bdb, 2)
                    bdb_3d{i} = bdb{1, i};
                end
                roomLayout = load(layout_path);
                roomLayout = roomLayout.layout';
                class_id = class_id + 1;
            else
                cam_R = data.Rtilt;
                bdb_3d = num2cell(SUNRGBDMetav2(imageId).groundtruth3DBB);
                roomLayout = data.gtCorner3D;
                class_id = zeros(1, length(bdb_3d));
                for i = 1:length(bdb_3d)
                    class_name = bdb_3d{i}.classname;
                    id = find(ismember(cls, class_name));
                    class_id(i) = id;
                    if id > 1
                        color = nyu_colorbox(id, :);
                        color = reshape(dec2hex(uint8(255.*color)).', 1, 6);
                        legend = [legend; {class_name, color}];
                    end
                end
            end
            %%

            % get 2d corners
            for kk = 1:length(bdb_3d)
                bdb_2d_temp = get_corners_of_bb3d(bdb_3d{kk});
                bdb_2d_temp = bdb_2d_temp(:, [1, 3, 2]);
                bdb_2d_temp(:, 2) = - bdb_2d_temp(:, 2);
                bdb_2d_temp = (data.K*inv(cam_R)*bdb_2d_temp')';
                bdb_2d_corner(:, 1) = bdb_2d_temp(:, 1) ./ bdb_2d_temp(:, 3);
                bdb_2d_corner(:, 2) = bdb_2d_temp(:, 2) ./ bdb_2d_temp(:, 3);
                bdb_2d{kk} = bdb_2d_corner;
            end
            
            iminfo = imfinfo(data.rgbpath);
%             %% draw 2D
%             fig = figure;
%             imshow(data.rgbpath);
%             hold on;
%             for kk =1:length(bdb_2d)
%                 if class_id(kk) == 1
%                     continue
%                 end
%                 draw_square_3d(bdb_2d{kk}, nyu_colorbox(class_id(kk), :), 5);
%                 x_max = min(bdb_2d{kk});
%         %         text(double(x_max(1)), double(x_max(2)), num2str(kk), 'Color','red','FontSize',14);
%             end
%             hold off;
%             set(gca,'xtick',[])
%             set(gca,'ytick',[])
%             axis off;
%             set(gca, 'position', [0 0 1 1 ]);
%             set(gcf, 'Position', [0,0,iminfo.Width,iminfo.Height]);
%             if dosave
%                 saveas(fig, fullfile(save_path, [num2str(imageId) '_bbox.png']));
%                 close(fig);
%             end

            %% draw 3D
            if ipath == 1
                f3d = figure;
                set(gcf, 'Position', [0,0,iminfo.Width,iminfo.Height]);
            else
                set(0, 'currentfigure', f3d);
                cla;
            end
            grid on
    %         set(gca,'xtick',[])
    %         set(gca,'ytick',[])
    %         set(gca, 'ztick', [])
    %         axis off;
            axis equal;
            if vis_pc
                vis_point_cloud(points3d,rgb, 100000);
            end
            hold on;
            for kk = 1:length(bdb_3d)
                if class_id(kk) == 1
                    continue
                end
                vis_cube(bdb_3d{kk}, nyu_colorbox(class_id(kk), :), 2);
            end
            Linewidth = 3;
            maxhight = 1.2;
            drawRoom(roomLayout,'b',Linewidth,maxhight);
            hold off;
    %         set(gca,'CameraViewAngleMode','manual')
            for iview = 1:length(views3d)
                view3d = views3d{iview};
                if ipath == 1
                    if strcmp(view3d, 'top')
                        view(0,90);
                    elseif strcmp(view3d, 'oblique')
                        view(-15,45);
                    end
                    viewsCam{iview} = {
                        get(gca, 'CameraPosition'), get(gca, 'CameraTarget'), ...
                        xlim, ylim, zlim};
                else
    %                 set(gca,'CameraViewAngleMode','manual')
                    newcp = viewsCam{iview}{1};
                    set(gca,'CameraPosition',newcp);
                    newct = viewsCam{iview}{2};
                    set(gca,'CameraTarget',newct);
                    xlim(viewsCam{iview}{3});
                    ylim(viewsCam{iview}{4});
                    zlim(viewsCam{iview}{5});
                end
    %             zoom(1)
    %             set(gca, 'position', [0 0 1 1 ]);
                if vis_pc
                    ext = '_3dpc.png';
                else
                    ext = '_3d.png';
                end
                if dosave
                    saveas(f3d, fullfile(save_path, [num2str(imageId) '_' view3d ext]));
                end
            end
        end
        if dosave
            close(f3d);
        end
        close all;
    catch
        warning([num2str(imageId) ' visualization failed']);
    end
end
legend = unique(legend)

