import numpy as np
from PIL import Image
import cv2
import glob
import csv


patches = "/scratch/users/rnd27/tcga_ovarian_cancer/patches/"
epi_stroma_masks = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/epi_stroma_masks/"
output_patches = "/scratch/users/rnd27/tcga_ovarian_cancer/patches_epithelium/"

masks = glob.glob(patches + "*")
masks = masks[20000:25000]
for mask in masks:
    filename = mask.split("/")[-1]
    patch = cv2.imread(mask)
    epi_stroma_mask = cv2.imread(epi_stroma_masks + filename, 0)

    count = sum(sum(epi_stroma_mask > 0))
    if float(count)/9000000 > 0.9:
        cv2.imwrite(output_patches + filename, patch)