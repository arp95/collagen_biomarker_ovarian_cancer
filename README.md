# Collagen Biomarker for Ovarian Cancer

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)](LICENSE.md)
---


## Introduction to the Project 
Understanding the role of Collagen biomarker in Ovarian Cancer.


## Authors
Arpit Aggarwal


### Packages Required
The packages required for running this code are PyTorch, Numpy, Scipy, Sklearn and OpenCV.


### Procedure for running the code
For nuclei segmentation task, use the file (code/nuclei_segmentation.py) for getting the segmented nuclei as the output. Use following steps for running this file.
1. Make sure your current directory looks like code/, results/ etc. In the results folder, make three sub-folders, nuclei_csvs, nuclei_overlay and nuclei_masks.
2. Remember to install packages, PyTorch, Numpy, Scipy, Sklearn and OpenCV for running this file.
3. In the nuclei_segmentation.py file, there are 10 hyper parameters, namely:
	a. input_path - specify your input path where all patches are there.
	b. output_path remains as 'results/'
	c. batch_size remains as 8
	d. workers remains as 0 if using cpu or 8 if you are using a gpu
	e. nr_types are 5 if using Monusac checkpoint while they are 6 if using Pannuke checkpoint
	f. model_type remains as 5
	g. model_path - specify the model path of Pannuke checkpoint or Monusac checkpoint. Download the model files from this link: Pannuke checkpoint model (https://drive.google.com/file/d/1SbSArI3KOOWHxRlxnjchO7_MbWzB4lNR/view) and Monusac checkpoint model (https://drive.google.com/file/d/13qkxDqv7CUqxN-l5CpeFVmc24mDw6CeV/view)
	h. device as 'cpu' or 'cuda'
	i. image_shape corresponds to the shape of the patches. If patch dimension is (1000x1000) use 1000 or if patch dimension is (3000x3000) use 3000
4. Run the file using this command, "python code/nuclei_segmentation.py"
5. This generates outputs for all patches in folders results/nuclei_overlay, results/nuclei_masks and results/nuclei_csvs. Please check these folders for the next steps.