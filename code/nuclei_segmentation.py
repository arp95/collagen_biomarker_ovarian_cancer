'''
Modified by: Arpit Aggarwal
Hovernet Model for Nuclei Segmentation: https://github.com/vqdang/hover_net
'''


# header files
import torch
import os
import copy

# parameters to update
input_path = 'results/patches_example/'
output_path = 'results/'
batch_size = 1
workers = 0
nr_types = 6
model_type = 'fast'
model_path = 'model_files/hovernet_fast_pannuke_type_tf2pytorch.tar'
json_path = "model_files/type_info_pannuke.json"
device = 'cpu'
image_shape = 3000


# main code
method_args = {'method' : {'model_args' : {'nr_types': nr_types,'mode': model_type,}, 'model_path' : model_path,}, 'type_info_path': json_path, 'device': device}
run_args = {'batch_size': batch_size, 'nr_inference_workers': workers, 'nr_post_proc_workers': 0, 'input_shape': image_shape}

# depending on model type
if model_type == 'fast':
    run_args['patch_input_shape'] = 256
    run_args['patch_output_shape'] = 164
else:
    run_args['patch_input_shape'] = 270
    run_args['patch_output_shape'] = 80

sub_cmd = 'tile'
if sub_cmd == 'tile':
    run_args.update({'input_dir': input_path, 'output_dir': output_path, 'mem_usage': float(0.1), 'draw_dot': False, 'save_qupath': False, 'save_raw_map': True,})

if sub_cmd == 'tile':
    print("Started with nuclei segmentation task...")
    from nuclei_model import *
    infer = InferManager(**method_args)
    infer.process_file_list(run_args)
    print("Done with nuclei segmentation task...")