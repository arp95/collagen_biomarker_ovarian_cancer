%% Author: Arpit Aggarwal
clc
close all
clear 
addpath(genpath('pwd'))

% HPC Paths
%files_dir = "/scratch/users/axa1399/tcga_ovarian_cancer/patches/";
%files = dir(fullfile(files_dir, '*.svs'));

% hard-coded paths for masks and images
files_dir = "../../ovarian_cancer_files/";
feature_maps_dir = "../../ovarian_cancer_results/collagen_feature_maps_final/sample_collagen_feature_maps/";
files = dir(fullfile(files_dir, '*.txt'));
feature_maps = dir(fullfile(feature_maps_dir, '*.mat'));

%% get feature map for a slide
for index = 1:length(files)
    filename = files(index).name;
    filename = extractBefore(filename, ".txt");
    filename

    file_feature_map = [];
    for index1 = 1:length(feature_maps)
        file_feature_map_index1 = feature_maps(index1).name;
        file_feature_map_index1 = extractBefore(file_feature_map_index1, ".mat");
        file_feature_map_index1_split = split(file_feature_map_index1, "_");
        row = file_feature_map_index1_split(2);
        col = file_feature_map_index1_split(3);
        row = cellfun(@str2num, row);
        col = cellfun(@str2num, col);
        matrix = load(feature_maps_dir + file_feature_map_index1 + ".mat");
        size_matrix = size(matrix.matrix);
        row = ((row / 3000) * size_matrix(1)) + 1;
        col = ((col / 3000) * size_matrix(2)) + 1;
        size_matrix
        file_feature_map(row:row+size_matrix(1)-1, col:col+size_matrix(2)-1) = matrix.matrix;
    end

    % update zero values in feature map to nan
    file_feature_map(file_feature_map == 0) = NaN;
    file_feature_map_avg = movmean(file_feature_map, 2, 1, 'omitnan', 'Endpoint', 'discard');
    file_feature_map_avg = movmean(file_feature_map_avg, 2, 2, 'omitnan', 'Endpoint', 'discard');

    % plot heatmap
    figure
    heat_map = heatmap(file_feature_map_avg, 'CellLabelColor','none');
    heat_map.GridVisible = 'off';
    heat_map.FontColor = 'none';
end