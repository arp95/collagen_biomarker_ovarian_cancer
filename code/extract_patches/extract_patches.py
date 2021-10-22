"""
Script for extracting patches.
"""


# header files needed
import openslide
import glob
from PIL import Image
import numpy as np
import sys
import os

# get the options selected by user
patch_dim = 1000
data_path = "data/"
output_path = "results/patches"

# get input files
files = glob.glob(data_path + "*.svs")

# read file in 40x/20x and generate the corresponding patch
for file in files:
    ts = openslide.OpenSlide(file)
    size0, size1, size2, size3 = ts.level_dimensions
    size = size0
    output_mask_name = f'{output_path}/{file}'
    output_mask_name  = output_mask_name.replace(f'/data', "")

    # read and save patch
    for index2 in range(0, int(size[1]), patch_dim):
        for index1 in range(0, int(size[0]), patch_dim):
            index1 = 60000
            index2 = 6000
            patch = ts.read_region((index1, index2), 0, (patch_dim, patch_dim)).convert('RGB')
            path = output_mask_name.replace(f'.svs', "_" + str(index2) + "_" + str(index1) + ".png")
            Image.fromarray(np.array(patch).astype(np.uint8)).save(path)
            break
        break