import cv2
import glob
import csv
import numpy as np


files_dir = "/mnt/rstor/CSE_BME_AXM788/data/TCGA_Ovarian Cancer/TCGA_Ovarian_Diagnostic_Path/"
epi_stroma_masks_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/epi_stroma_masks/"
nuclei_masks_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/nuclei_masks/"
macrophage_masks_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/macrophage_nuclei_masks/"
output_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/macrophage_output/"

files = glob.glob(files_dir + "*")
files = files[30:40]
macrophage_masks = glob.glob(macrophage_masks_dir + "*")
for file in files:
    filename = file.split("/")[-1][:-4]
    count_macrophage = 0.0
    count_total = 0.0
    count1_macrophage = 0.0
    count1_total = 0.0
    for macrophage_mask in macrophage_masks:
        mask_filename = macrophage_mask.split("/")[-1]
        if filename in macrophage_mask:
            macrophage_image = cv2.imread(macrophage_mask, 0)
            nuclei_image = cv2.imread(nuclei_masks_dir + mask_filename, 0)
            epi_stroma_image = cv2.imread(epi_stroma_masks_dir + mask_filename, 0)
            
            # count macrophage
            _, macrophage_image_thresh = cv2.threshold(macrophage_image, 1, 255, cv2.THRESH_BINARY)
            cnts = cv2.findContours(macrophage_image_thresh.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            cnts = cnts[0]
            for cnt in cnts:
                (x, y, w, h) = cv2.boundingRect(cnt)
                area = cv2.contourArea(cnt)
                count1_macrophage += 1
                
                count1 = 0
                for index1 in range(max(y, 0), min(y+h, 3000)):
                    for index2 in range(max(x, 0), min(x+w, 3000)):
                        if epi_stroma_image[index1, index2] > 0:
                            count1 += 1
                if count1 > 0.95*area:
                    count_macrophage += 1
            
            # count nuclei
            _, nuclei_image_thresh = cv2.threshold(nuclei_image, 1, 255, cv2.THRESH_BINARY)
            cnts = cv2.findContours(nuclei_image_thresh.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            cnts = cnts[0]
            for cnt in cnts:
                (x, y, w, h) = cv2.boundingRect(cnt)
                area = cv2.contourArea(cnt)
                count1_total += 1
                
                count1 = 0
                for index1 in range(max(y, 0), min(y+h, 3000)):
                    for index2 in range(max(x, 0), min(x+w, 3000)):
                        if epi_stroma_image[index1, index2] > 0:
                            count1 += 1
                if count1 > 0.95*area:
                    count_total += 1
    if count_total > 0 and count1_total > 0:
        value = float(count_macrophage) / count_total
        value_1 = float(count1_macrophage) / count1_total
        with open(output_dir + filename + ".csv", 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile)
            spamwriter.writerow([value, value_1])