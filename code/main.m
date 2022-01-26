%% Author: Arpit Aggarwal
clc
close all
clear 
addpath(genpath('pwd'))


% hard-coded paths for masks and images
patches_dir = "../../ovarian_cancer_results/patches_final/";
patches = dir(fullfile(patches_dir, '*.png'));
epi_stroma_masks_dir = "../../ovarian_cancer_results/epi_stroma_masks_final/";
nuclei_masks_dir = "../../ovarian_cancer_results/nuclei_masks_final/";
collagen_masks_dir = "../../ovarian_cancer_results/collagen_fiber_masks_final/";

%% get collagen mask for each patch
for index = 1:1
    filename = patches(index).name
    current_patch = imread(patches_dir + filename);
    epi_stroma_mask = imread(epi_stroma_masks_dir + filename);
    nuclei_mask = imread(nuclei_masks_dir + filename);

    % hyperparameters for calculating collagen features
    win_size = 200;
    filter_scale = 3;
    orient_cooccur_scheme = 1;
    feature_descriptor = 6;
    orient_bin_interval = 10;
    orient_num = 180 / orient_bin_interval;

    %% extract collagen fiber mask
    frag_thresh = filter_scale*5; %remove detected collagen fragments with an area lower than the predefined threshold
    [bifs] = compute_bifs(current_patch, filter_scale, .1, 1); %use BIF based model to extract the linear structures which was used as the representative of collagen fibers
    collagen_mask = bifs == feature_descriptor;
    [height, width] = size(collagen_mask);       
    collagen_mask = (collagen_mask & (1 - epi_stroma_mask));
    collagen_mask = (collagen_mask & (1 - nuclei_mask));
    collagen_mask = bwareaopen(collagen_mask, frag_thresh);
    %patch_collagen_mask = labeloverlay(current_patch, collagen_mask, 'transparency', 0, 'Colormap', [0,0,1]);
    %imwrite(patch_collagen_mask, collagen_masks_dir + filename);

    %% collagen centroid and orientation information extraction
    collogen_props = regionprops('table', collagen_mask, 'Centroid', 'Orientation', 'Area');
    colg_center = collogen_props.Centroid;
    colg_area = collogen_props.Area;
    colg_orient = collogen_props.Orientation;

    %% feature extraction
    colg_orient_bin = fix(colg_orient / orient_bin_interval);
    colg_orient_bin = colg_orient_bin + 9;
    win_size_ind = 0;
    step_size = win_size / 2;
    win_x_ind = 0;
    for win_x = 1:step_size:width-win_size+1
        win_x_ind = win_x_ind+1;
        win_y_ind = 0;
        for win_y = 1:step_size:height-win_size+1
            win_y_ind = win_y_ind+1;
            p_orient_occur = [];
            inwin_colg_ind = find(colg_center(:,1)>=win_x & colg_center(:,1)<win_x+win_size-1 & colg_center(:,2)>=win_y & colg_center(:,2)<win_y+win_size-1);
            inwin_epi_stroma_mask = epi_stroma_mask(win_y:win_y+win_size-1, win_x:win_x+win_size-1);
            number_of_zeros = sum(inwin_epi_stroma_mask(:) == 0);
            if number_of_zeros >= 12000
                inwin_colg_orient = colg_orient_bin(inwin_colg_ind); 
                inwin_colg_area = colg_area(inwin_colg_ind);
                if length(inwin_colg_orient)>=5
                    [orient_occur_feats] = disorder_feat_extract(inwin_colg_orient, inwin_colg_area, orient_num, orient_cooccur_scheme);
                    if isfield(orient_occur_feats, 'val')
                        cfod_map(win_y_ind, win_x_ind, :) = orient_occur_feats.val;
                    end 
                end
            end
        end
    end

    cfod_map_final = [];
    for feat_ind = 1:13
        cfod_map_move_avg = movmean(cfod_map(:,:,feat_ind), 2, 1, 'omitnan', 'Endpoint', 'discard');
        cfod_map_final(:, :, feat_ind) = movmean(cfod_map_move_avg, 2, 2, 'omitnan', 'Endpoint', 'discard');
    end
    size(current_patch)
    size(cfod_map_final(:,:,5))
    writematrix(cfod_map_final(:,:,5), "sample.txt");

    % plot heatmap
    figure();
    heat_map = heatmap(cfod_map_final(:,:,5));
    heat_map.GridVisible = 'off';
    colormap(jet);
    heat_map.ColorbarVisible = 'off';
    heat_map.FontColor = 'none';
    exportgraphics(gcf, filename);
end
%% after acquisation of the feature map, a set of statistics e.g. mean, std, sknewness could be calculated. 