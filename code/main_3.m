%% Author: Arpit Aggarwal
clc
close all
clear 
addpath(genpath('pwd'))

% hard-coded paths for masks and images
files_dir = "../../ovarian_cancer_results/collagen_feature_maps_final/high_risk/";
files = dir(fullfile(files_dir, '*.mat'));

%% get feature map for a slide
for index = 1:1
    filename = files(index).name;
    filename = extractBefore(filename, ".mat");
    filename = "TCGA-42-2593"
    filename

    % load matrix and find superpixels
    matrix = load(files_dir + filename + ".mat");
    file_feature_map = matrix.file_feature_map;
    file_feature_map(isnan(file_feature_map)) = 0;
    [L, N] = superpixels(file_feature_map, 100);

    file_feature_map(file_feature_map == 0) = NaN;
    output_image = zeros(size(file_feature_map), 'like', file_feature_map);
    idx = label2idx(L);
    num_rows = size(file_feature_map, 1);
    num_cols = size(file_feature_map, 2);
    for labelVal = 1:N
        output_image(idx{labelVal}) = mean(file_feature_map(idx{labelVal}), 'omitnan');
    end

    figure
    imshow(output_image)

    % plot heatmap
    figure
    BW = boundarymask(L);
    imshow(BW);
    %heat_map = heatmap(file_feature_map, 'CellLabelColor','none');
    %heat_map.GridVisible = 'off';
    %heat_map.FontColor = 'none';
end