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


patches_dir = "/scratch/users/axa1399/upmc_ovarian_cancer/patches/"
output_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/upmc_ovarian_cancer/macrophage_nuclei_interim_masks_1/"
model_path = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/bestmodel_unet_2.pth"


# load model
class conv_block_nested(nn.Module):

    def __init__(self, in_ch, mid_ch, out_ch):
        super(conv_block_nested, self).__init__()
        self.activation = nn.ReLU(inplace=True)
        self.conv1 = nn.Conv2d(in_ch, mid_ch, kernel_size=3, padding=1, bias=True)
        self.bn1 = nn.BatchNorm2d(mid_ch)
        self.conv2 = nn.Conv2d(mid_ch, out_ch, kernel_size=3, padding=1, bias=True)
        self.bn2 = nn.BatchNorm2d(out_ch)

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.activation(x)

        x = self.conv2(x)
        x = self.bn2(x)
        output = self.activation(x)

        return output

class Nested_UNet(nn.Module):

    def __init__(self, in_ch=3, out_ch=1):
        super(Nested_UNet, self).__init__()

        n1 = 64
        filters = [n1, n1 * 2, n1 * 4, n1 * 8, n1 * 16]

        self.pool = nn.MaxPool2d(kernel_size=2, stride=2)
        self.Up = nn.Upsample(scale_factor=2, mode='bilinear', align_corners=True)

        self.conv0_0 = conv_block_nested(in_ch, filters[0], filters[0])
        self.conv1_0 = conv_block_nested(filters[0], filters[1], filters[1])
        self.conv2_0 = conv_block_nested(filters[1], filters[2], filters[2])
        self.conv3_0 = conv_block_nested(filters[2], filters[3], filters[3])
        self.conv4_0 = conv_block_nested(filters[3], filters[4], filters[4])

        self.conv0_1 = conv_block_nested(filters[0] + filters[1], filters[0], filters[0])
        self.conv1_1 = conv_block_nested(filters[1] + filters[2], filters[1], filters[1])
        self.conv2_1 = conv_block_nested(filters[2] + filters[3], filters[2], filters[2])
        self.conv3_1 = conv_block_nested(filters[3] + filters[4], filters[3], filters[3])

        self.conv0_2 = conv_block_nested(filters[0]*2 + filters[1], filters[0], filters[0])
        self.conv1_2 = conv_block_nested(filters[1]*2 + filters[2], filters[1], filters[1])
        self.conv2_2 = conv_block_nested(filters[2]*2 + filters[3], filters[2], filters[2])

        self.conv0_3 = conv_block_nested(filters[0]*3 + filters[1], filters[0], filters[0])
        self.conv1_3 = conv_block_nested(filters[1]*3 + filters[2], filters[1], filters[1])

        self.conv0_4 = conv_block_nested(filters[0]*4 + filters[1], filters[0], filters[0])

        self.final = nn.Conv2d(filters[0], out_ch, kernel_size=1)

    def forward(self, x):

        x0_0 = self.conv0_0(x)
        x1_0 = self.conv1_0(self.pool(x0_0))
        x0_1 = self.conv0_1(torch.cat([x0_0, self.Up(x1_0)], 1))

        x2_0 = self.conv2_0(self.pool(x1_0))
        x1_1 = self.conv1_1(torch.cat([x1_0, self.Up(x2_0)], 1))
        x0_2 = self.conv0_2(torch.cat([x0_0, x0_1, self.Up(x1_1)], 1))

        x3_0 = self.conv3_0(self.pool(x2_0))
        x2_1 = self.conv2_1(torch.cat([x2_0, self.Up(x3_0)], 1))
        x1_2 = self.conv1_2(torch.cat([x1_0, x1_1, self.Up(x2_1)], 1))
        x0_3 = self.conv0_3(torch.cat([x0_0, x0_1, x0_2, self.Up(x1_2)], 1))

        x4_0 = self.conv4_0(self.pool(x3_0))
        x3_1 = self.conv3_1(torch.cat([x3_0, self.Up(x4_0)], 1))
        x2_2 = self.conv2_2(torch.cat([x2_0, x2_1, self.Up(x3_1)], 1))
        x1_3 = self.conv1_3(torch.cat([x1_0, x1_1, x1_2, self.Up(x2_2)], 1))
        x0_4 = self.conv0_4(torch.cat([x0_0, x0_1, x0_2, x0_3, self.Up(x1_3)], 1))

        output = self.final(x0_4)
        return output


# evaluation code
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device = 'cpu'
model = Nested_UNet()
model.load_state_dict(torch.load(model_path, map_location=torch.device(device)))
model = model.to(device)
model.eval()


# transform for input patch
img_transform = torchvision.transforms.Compose([
        torchvision.transforms.Resize((1500, 1500)),
        torchvision.transforms.ToTensor(),
        torchvision.transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])


patches = glob.glob(patches_dir + "*")
patches = patches[23500:25000]
for patch in patches:
    filename = patch.split("/")[-1]
    image = Image.open(patch).convert('RGB')
    image = img_transform(image)
    image = np.array(image)
    output_mask = np.zeros((1500, 1500))
    
    for index1 in range(0, image.shape[1]-128, 128):
        for index2 in range(0, image.shape[2]-128, 128):
            current_patch = image[:, index1:index1+128, index2:index2+128]
            current_patch = torch.from_numpy(np.array([current_patch]))
            current_patch = current_patch.to(device, dtype=torch.float32)
            output_patch = torch.sigmoid(model(current_patch))
            output_patch = output_patch.squeeze()
            output_patch = output_patch.squeeze()
            output_patch = output_patch > 0.8
            output_mask[index1:index1+128, index2:index2+128] = output_patch
    output_mask = np.array(output_mask, dtype=np.uint8)*255
    output_mask = cv2.resize(output_mask, (3000, 3000))
    cv2.imwrite(output_dir+filename, output_mask)