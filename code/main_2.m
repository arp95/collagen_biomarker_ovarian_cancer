%% Author: Arpit Aggarwal
clc
close all
clear 
addpath(genpath('pwd'))


% HPC Paths
files_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/collagen_feature_maps_200_final/";
feature_maps_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/collagen_feature_maps_450_2/";
files = dir(fullfile(files_dir, '*.mat'));
feature_maps = dir(fullfile(feature_maps_dir, '*.mat'));
collagen_masks_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/collagen_feature_maps_450_2_final/";

% hard-coded paths for masks and images
%files_dir = "../../ovarian_cancer_files/";
%feature_maps_dir = "../../ovarian_cancer_results/sample_collagen_masks/";
%collagen_masks_dir = "../../ovarian_cancer_results/final_collagen_feature_maps/";
%files = dir(fullfile(files_dir, '*.txt'));
%feature_maps = dir(fullfile(feature_maps_dir, '*.mat'));


%% get feature map for a slide
for index = 1:length(files)
    filename = files(index).name;
    filename = extractBefore(filename, ".mat");
    filename

    count = 0;
    sum1 = 0;
    max_file = -1000000;
    min_file = 1000000;
    mean_file = [];
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
            sum1 = sum1 + mean(matrix.matrix, 'all', 'omitnan');
            mean_file = [mean_file, mean(matrix.matrix, 'all', 'omitnan')];
            min_val = min(matrix.matrix, [], 'all', 'omitnan');
            max_val = max(matrix.matrix, [], 'all', 'omitnan');
            min_file = min(min_file, min_val);
            max_file = max(max_file, max_val);
        end
    end

    % update zero values in feature map to nan
    %mean_file = sum1/count;
    %range_file = max_file - min_file;


    % new code
    feature_1 = mean(mean_file);
    feature_4 = min_file;
    feature_5 = max_file;
    feature_6 = feature_5 - feature_4;
    feature_matrix = [feature_1, feature_4, feature_5, feature_6];
    if feature_5 > 0
        csvwrite(collagen_masks_dir + filename + '.csv', feature_matrix);
    end
    % plot heatmap
    %figure
    %heat_map = heatmap(file_feature_map_avg, 'CellLabelColor','none');
    %heat_map.GridVisible = 'off';
    %heat_map.FontColor = 'none';
end