import cv2
import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision
from torch import Tensor
from torchvision.transforms.transforms import Resize
import random
import numpy as np
from random import shuffle
import glob
from PIL import Image


patches_dir = "/scratch/users/axa1399/dl_predict_outcome/test/new_input/low/"
output_dir = "/scratch/users/axa1399/dl_predict_outcome/test/results_resnet18/"
model_path = "/scratch/users/axa1399/dl_predict_outcome/dl_predict_outcome_best_model_resnet18.pth"


# evaluation code
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device = 'cpu'
model = torchvision.models.resnet18(pretrained=True)
model.fc = torch.nn.Sequential(
    torch.nn.Linear(512, 2, bias=True),
    torch.nn.LogSoftmax()
)
model.load_state_dict(torch.load(model_path, map_location=torch.device(device)))
model = model.to(device)
model.eval()


# transform for input patch
img_transform = torchvision.transforms.Compose([
        torchvision.transforms.CenterCrop((512, 512)),
        torchvision.transforms.ToTensor(),
        torchvision.transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])


patches = glob.glob(patches_dir + "*")
#patches = patches[75000:]
for patch in patches:
    filename = patch.split("/")[-1]
    image = Image.open(patch).convert('RGB')
    image = img_transform(image)
    image = np.array(image)
    image_to_save = cv2.imread(patch)
    
    current_patch = torch.from_numpy(np.array([image]))
    current_patch = current_patch.to(device, dtype=torch.float32)
    output_patch = model(current_patch)
    _, predicted = output_patch.max(1)
    pred = int(predicted[0])
    if pred == 0:
        cv2.imwrite(output_dir + 'low/' + filename, image_to_save)
    else:
        cv2.imwrite(output_dir + 'high/' + filename, image_to_save)