%% Author: Arpit Aggarwal
clc
close all
clear 
addpath(genpath('pwd'))

% hard-coded paths
patches_dir = "../../tcga_ovarian_cancer_results/example_patches/";
feature_map_dir = "../../tcga_ovarian_cancer_results/example_collagen_results/";
patches = dir(fullfile(patches_dir, '*.png'));
for index = 1:1
    filename = patches(index).name;
    filename
    current_patch = imread(patches_dir + filename);

    filename = extractBefore(filename, ".png");
    matrix = load(feature_map_dir + filename + ".mat");
    file_feature_map = matrix.matrix;

    figure
    roi = imresize(isnan(file_feature_map), [3000, 3000], 'nearest');
    imshow(current_patch);
    hold on
    handle = imagesc(imresize(file_feature_map, [3000, 3000], 'nearest'), [1.5 2.5]);
    set(handle, 'AlphaData', 0.5*(1-double(roi)));
    colorbar
end