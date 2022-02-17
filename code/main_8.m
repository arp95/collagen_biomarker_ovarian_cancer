%% Author: Arpit Aggarwal
clc
close all
clear
addpath(genpath('pwd'))


% hard-coded paths
files_dir_1 = "../../ovarian_cancer_results/collagen_final/collagen_feature_maps_200_final/";
files_1 = dir(fullfile(files_dir_1, '*.mat'));
files_dir_2 = "../../ovarian_cancer_results/collagen_final/collagen_feature_maps_250_final/";
files_2 = dir(fullfile(files_dir_2, '*.mat'));
files_dir_3 = "../../ovarian_cancer_results/collagen_final/collagen_feature_maps_300_final/";
files_3 = dir(fullfile(files_dir_3, '*.mat'));
files_dir_4 = "../../ovarian_cancer_results/collagen_final/collagen_feature_maps_350_final/";
files_4 = dir(fullfile(files_dir_4, '*.mat'));
files_dir_5 = "../../ovarian_cancer_results/collagen_final/collagen_feature_maps_400_final/";
files_5 = dir(fullfile(files_dir_5, '*.mat'));
files_dir_6 = "../../ovarian_cancer_results/collagen_final/collagen_feature_maps_450_final/";
files_6 = dir(fullfile(files_dir_6, '*.mat'));
files_dir_7 = "../../ovarian_cancer_results/collagen_final/collagen_feature_maps_500_final/";
files_7 = dir(fullfile(files_dir_7, '*.mat'));
files_dir_8 = "../../ovarian_cancer_results/collagen_final/collagen_feature_maps_550_final/";
files_8 = dir(fullfile(files_dir_8, '*.mat'));
files_dir_9 = "../../ovarian_cancer_results/collagen_final/collagen_feature_maps_600_final/";
files_9 = dir(fullfile(files_dir_9, '*.mat'));


%% get feature map for a slide
for index = 1:length(files_1)
    filename = files_1(index).name;
    filename = extractBefore(filename, ".mat");
    filename

    % load matrix
    matrix_1 = load(files_dir_1 + filename + ".mat");
    file_feature_map_1 = matrix_1.file_feature_map;
    mean_1 = mean(file_feature_map_1, 'all', 'omitnan');
    median_1 = median(file_feature_map_1, 'all', 'omitnan');
    min_1 = min(file_feature_map_1, [], 'all', 'omitnan');
    max_1 = max(file_feature_map_1, [], 'all', 'omitnan');
    range_1 = max_1 - min_1;
    std_1 = std(file_feature_map_1, 0, 'all', 'omitnan');
    skewness_1 = skewness(file_feature_map_1, 1, 'all');
    kurtosis_1 = kurtosis(file_feature_map_1, 1, 'all');

    matrix_2 = load(files_dir_2 + filename + ".mat");
    file_feature_map_2 = matrix_2.file_feature_map;
    mean_2 = mean(file_feature_map_2, 'all', 'omitnan');
    median_2 = median(file_feature_map_2, 'all', 'omitnan');
    min_2 = min(file_feature_map_2, [], 'all', 'omitnan');
    max_2 = max(file_feature_map_2, [], 'all', 'omitnan');
    range_2 = max_2 - min_2;
    std_2 = std(file_feature_map_2, 0, 'all', 'omitnan');
    skewness_2 = skewness(file_feature_map_2, 1, 'all');
    kurtosis_2 = kurtosis(file_feature_map_2, 1, 'all');

    matrix_3 = load(files_dir_3 + filename + ".mat");
    file_feature_map_3 = matrix_3.file_feature_map;
    mean_3 = mean(file_feature_map_3, 'all', 'omitnan');
    median_3 = median(file_feature_map_3, 'all', 'omitnan');
    min_3 = min(file_feature_map_3, [], 'all', 'omitnan');
    max_3 = max(file_feature_map_3, [], 'all', 'omitnan');
    range_3 = max_3 - min_3;
    std_3 = std(file_feature_map_3, 0, 'all', 'omitnan');
    skewness_3 = skewness(file_feature_map_3, 1, 'all');
    kurtosis_3 = kurtosis(file_feature_map_3, 1, 'all');

    matrix_4 = load(files_dir_4 + filename + ".mat");
    file_feature_map_4 = matrix_4.file_feature_map;
    mean_4 = mean(file_feature_map_4, 'all', 'omitnan');
    median_4 = median(file_feature_map_4, 'all', 'omitnan');
    min_4 = min(file_feature_map_4, [], 'all', 'omitnan');
    max_4 = max(file_feature_map_4, [], 'all', 'omitnan');
    range_4 = max_4 - min_4;
    std_4 = std(file_feature_map_4, 0, 'all', 'omitnan');
    skewness_4 = skewness(file_feature_map_4, 1, 'all');
    kurtosis_4 = kurtosis(file_feature_map_4, 1, 'all');

    matrix_5 = load(files_dir_5 + filename + ".mat");
    file_feature_map_5 = matrix_5.file_feature_map;
    mean_5 = mean(file_feature_map_5, 'all', 'omitnan');
    median_5 = median(file_feature_map_5, 'all', 'omitnan');
    min_5 = min(file_feature_map_5, [], 'all', 'omitnan');
    max_5 = max(file_feature_map_5, [], 'all', 'omitnan');
    range_5 = max_5 - min_5;
    std_5 = std(file_feature_map_5, 0, 'all', 'omitnan');
    skewness_5 = skewness(file_feature_map_5, 1, 'all');
    kurtosis_5 = kurtosis(file_feature_map_5, 1, 'all');

    matrix_6 = load(files_dir_6 + filename + ".mat");
    file_feature_map_6 = matrix_6.file_feature_map;
    mean_6 = mean(file_feature_map_6, 'all', 'omitnan');
    median_6 = median(file_feature_map_6, 'all', 'omitnan');
    min_6 = min(file_feature_map_6, [], 'all', 'omitnan');
    max_6 = max(file_feature_map_6, [], 'all', 'omitnan');
    range_6 = max_6 - min_6;
    std_6 = std(file_feature_map_6, 0, 'all', 'omitnan');
    skewness_6 = skewness(file_feature_map_6, 1, 'all');
    kurtosis_6 = kurtosis(file_feature_map_6, 1, 'all');

    matrix_7 = load(files_dir_7 + filename + ".mat");
    file_feature_map_7 = matrix_7.file_feature_map;
    mean_7 = mean(file_feature_map_7, 'all', 'omitnan');
    median_7 = median(file_feature_map_7, 'all', 'omitnan');
    min_7 = min(file_feature_map_7, [], 'all', 'omitnan');
    max_7 = max(file_feature_map_7, [], 'all', 'omitnan');
    range_7 = max_7 - min_7;
    std_7 = std(file_feature_map_7, 0, 'all', 'omitnan');
    skewness_7 = skewness(file_feature_map_7, 1, 'all');
    kurtosis_7 = kurtosis(file_feature_map_7, 1, 'all');

    matrix_8 = load(files_dir_8 + filename + ".mat");
    file_feature_map_8 = matrix_8.file_feature_map;
    mean_8 = mean(file_feature_map_8, 'all', 'omitnan');
    median_8 = median(file_feature_map_8, 'all', 'omitnan');
    min_8 = min(file_feature_map_8, [], 'all', 'omitnan');
    max_8 = max(file_feature_map_8, [], 'all', 'omitnan');
    range_8 = max_8 - min_8;
    std_8 = std(file_feature_map_8, 0, 'all', 'omitnan');
    skewness_8 = skewness(file_feature_map_8, 1, 'all');
    kurtosis_8 = kurtosis(file_feature_map_8, 1, 'all');

    matrix_9 = load(files_dir_9 + filename + ".mat");
    file_feature_map_9 = matrix_9.file_feature_map;
    mean_9 = mean(file_feature_map_9, 'all', 'omitnan');
    median_9 = median(file_feature_map_9, 'all', 'omitnan');
    min_9 = min(file_feature_map_9, [], 'all', 'omitnan');
    max_9 = max(file_feature_map_9, [], 'all', 'omitnan');
    range_9 = max_9 - min_9;
    std_9 = std(file_feature_map_9, 0, 'all', 'omitnan');
    skewness_9 = skewness(file_feature_map_9, 1, 'all');
    kurtosis_9 = kurtosis(file_feature_map_9, 1, 'all');

    % form feature matrix and save
    feature_matrix = [mean_1, mean_2, mean_3, mean_4, mean_5, mean_6, mean_7, mean_8, mean_9, std_1, std_2, std_3, std_4, std_5, std_6, std_7, std_8, std_9, median_1, median_2, median_3, median_4, median_5, median_6, median_7, median_8, median_9, min_1, min_2, min_3, min_4, min_5, min_6, min_7, min_8, min_9, max_1, max_2, max_3, max_4, max_5, max_6, max_7, max_8, max_9, skewness_1, skewness_2, skewness_3, skewness_4, skewness_5, skewness_6, skewness_7, skewness_8, skewness_9, range_1, range_2, range_3, range_4, range_5, range_6, range_7, range_8, range_9, kurtosis_1, kurtosis_2, kurtosis_3, kurtosis_4, kurtosis_5, kurtosis_6, kurtosis_7, kurtosis_8, kurtosis_9];
    csvwrite("../../ovarian_cancer_results/collagen_final/features/" + filename + '.csv', feature_matrix);
end