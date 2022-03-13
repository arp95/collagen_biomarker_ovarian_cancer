%% Author: Arpit Aggarwal
clc
close all
clear 
addpath(genpath('pwd'))


% HPC Paths
files_dir = "/mnt/rstor/CSE_BME_AXM788/data/TCGA_Ovarian Cancer/TCGA_Ovarian_Diagnostic_Path/";
feature_maps_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/collagen_feature_maps_250/";
files = dir(fullfile(files_dir, '*.svs'));
feature_maps = dir(fullfile(feature_maps_dir, '*.mat'));
collagen_masks_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/sample_2/";

% hard-coded paths for masks and images
%files_dir = "../../ovarian_cancer_files/";
%feature_maps_dir = "../../ovarian_cancer_results/sample_collagen_feature_maps/";
%collagen_masks_dir = "../../ovarian_cancer_results/final_collagen_feature_maps/";
%files = dir(fullfile(files_dir, '*.txt'));
%feature_maps = dir(fullfile(feature_maps_dir, '*.mat'));


%% get feature map for a slide
for index = 1:length(files)
    filename = files(index).name;
    filename = extractBefore(filename, ".svs");
    filename

    file_feature_map = [];
    count = 0;
    sum = 0;
    max_file = -1000000;
    min_file = 1000000;
    for index1 = 1:length(feature_maps)
        file_feature_map_index1 = feature_maps(index1).name;
        file_feature_map_index1 = extractBefore(file_feature_map_index1, ".mat");
        file_feature_map_index1_split = split(file_feature_map_index1, "_");

        if strcmp(file_feature_map_index1_split{1}, filename)
            row = file_feature_map_index1_split(2);
            col = file_feature_map_index1_split(3);
            row = cellfun(@str2num, row);
            col = cellfun(@str2num, col);
            matrix = load(feature_maps_dir + file_feature_map_index1 + ".mat");
            count = count + 1;
            sum = sum + mean(matrix.matrix, 'all', 'omitnan');
            min_val = min(matrix.matrix, [], 'all', 'omitnan');
            max_val = max(matrix.matrix, [], 'all', 'omitnan');
            min_file = min(min_file, min_val);
            max_file = max(max_file, max_val);
            size_matrix = size(matrix.matrix);
            row = ((row / 3000) * 20) + 1;
            col = ((col / 3000) * 20) + 1;
            file_feature_map(row:row+size_matrix(1)-1, col:col+size_matrix(2)-1) = matrix.matrix;
        end
    end

    % update zero values in feature map to nan
    mean_file = sum/count;
    range_file = max_file - min_file;
    feature_matrix = [mean_file, min_file, max_file, range_file];
    csvwrite(collagen_masks_dir + filename + '.csv', feature_matrix);
    file_feature_map(file_feature_map == 0) = NaN;
    %save(collagen_masks_dir + filename + '.mat', "file_feature_map");

    % plot heatmap
    %figure
    %heat_map = heatmap(file_feature_map_avg, 'CellLabelColor','none');
    %heat_map.GridVisible = 'off';
    %heat_map.FontColor = 'none';
end