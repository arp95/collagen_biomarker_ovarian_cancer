import numpy as np
from PIL import Image
import cv2
import glob
import csv


macrophage_masks = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/upmc_ovarian_cancer/macrophage_nuclei_interim_masks_1/"
nuclei_masks = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/upmc_ovarian_cancer/nuclei_masks/"
output_masks = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/upmc_ovarian_cancer/macrophage_nuclei_masks_1/"

masks = glob.glob(macrophage_masks + "*")
masks = masks[15000:20000]
for mask in masks:
    filename = mask.split("/")[-1]
    nuclei_mask = cv2.imread(nuclei_masks + filename, 0)
    macrophage_mask = cv2.imread(mask, 0)
    
    if nuclei_mask is None:
        continue

    _, nuclei_mask_thresh = cv2.threshold(nuclei_mask, 1, 255, cv2.THRESH_BINARY)
    cnts = cv2.findContours(nuclei_mask_thresh.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    cnts = cnts[0]
    final_mask = np.zeros(nuclei_mask_thresh.shape)
    for cnt in cnts:
        (x, y, w, h) = cv2.boundingRect(cnt)
        area = cv2.contourArea(cnt)
        
        count1 = 0
        for index1 in range(max(y-5, 0), min(y+h+5, 3000)):
            for index2 in range(max(x-5, 0), min(x+w+5, 3000)):
                if macrophage_mask[index1, index2] > 0:
                    count1 += 1
        if count1 > 0.25*area:
            cv2.fillPoly(final_mask, pts=[cnt], color=(255, 255, 255))
    cv2.imwrite(output_masks+filename, final_mask)