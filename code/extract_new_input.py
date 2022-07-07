# header files
import cv2
import os
import glob
import numpy as np


# paths
input_dir = glob.glob("/scratch/users/axa1399/dl_predict_outcome/test/input/high/*")
input_dir = input_dir[125000:150000]
epi_stroma_dir = "/scratch/users/axa1399/dl_predict_outcome/test/epi_stroma_masks/high/"
new_input_dir = "/scratch/users/axa1399/dl_predict_outcome/test/new_input/high/"

for index in range(0, len(input_dir)):
    input_filename = input_dir[index].split("/")[-1]
    epi_stroma_filename = epi_stroma_dir + input_filename
    
    # count pixels belonging to stroma
    epi_stroma_image = cv2.imread(epi_stroma_filename, 0)
    sum_zeros = sum(sum(epi_stroma_image == 0))
    if sum_zeros / 262144.0 > 0.75:
        input_image = cv2.imread(input_dir[index])
        cv2.imwrite(new_input_dir + input_filename, input_image)