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


patches_dir = "/scratch/users/rnd27/tcga_ovarian_cancer/patches/"
output_dir = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/macrophage_nuclei_interim_masks/"
model_path = "/mnt/rstor/CSE_BME_AXM788/home/axa1399/tcga_ovarian_cancer/bestmodel_unet.pth"


# load model
class UNet(nn.Module):
    def __init__(self, in_channels=3, n_classes=1, depth=4, wf=4, padding=True, batch_norm=True, up_mode='upconv' ,concat=True):
        """
        Implementation of
        U-Net: Convolutional Networks for Biomedical Image Segmentation
        (Ronneberger et al., 2015)
        https://arxiv.org/abs/1505.04597

        Using the default arguments will yield the exact version used
        in the original paper

        Args:
            in_channels (int): number of input channels
            n_classes (int): number of output channels
            depth (int): depth of the network
            wf (int): number of filters in the first layer is 2**wf
            padding (bool): if True, apply padding such that the input shape
                            is the same as the output.
                            This may introduce artifacts
            batch_norm (bool): Use BatchNorm after layers with an
                               activation function
            up_mode (str): one of 'upconv' or 'upsample'.
                           'upconv' will use transposed convolutions for
                           learned upsampling.
                           'upsample' will use bilinear upsampling.
        """
        super(UNet, self).__init__()
        assert up_mode in ('upconv', 'upsample')
        self.padding = padding
        self.depth = depth
        self.concat = concat
        prev_channels = in_channels
        self.down_path = nn.ModuleList()
        for i in range(depth):
            self.down_path.append(UNetConvBlock(prev_channels, 2**(wf+i),
                                                padding, batch_norm))
            prev_channels = 2**(wf+i)

        self.up_path = nn.ModuleList()
        for i in reversed(range(depth - 1)):
            self.up_path.append(UNetUpBlock(prev_channels, 2**(wf+i), up_mode,
                                            padding, batch_norm , concat))
            prev_channels = 2**(wf+i)
        self.last = nn.Conv2d(prev_channels, n_classes, kernel_size=1)

    def forward(self, x):
        blocks = []
        for i, down in enumerate(self.down_path):
            x = down(x)
            if i != len(self.down_path)-1:
                blocks.append(x)
                x = F.max_pool2d(x, 2)

        for i, up in enumerate(self.up_path):
            x = up(x, blocks[-i-1])
        return self.last(x)


class UNetConvBlock(nn.Module):
    def __init__(self, in_size, out_size, padding, batch_norm):
        super(UNetConvBlock, self).__init__()
        block = []

        block.append(nn.Conv2d(in_size, out_size, kernel_size=3,
                               padding=int(padding)))
        if batch_norm:
            block.append(nn.BatchNorm2d(out_size))
        block.append(nn.ReLU())

        block.append(nn.Conv2d(out_size, out_size, kernel_size=3,
                               padding=int(padding)))
        if batch_norm:
            block.append(nn.BatchNorm2d(out_size))
        block.append(nn.ReLU())
        self.block = nn.Sequential(*block)

    def forward(self, x):
        out = self.block(x)
        return out


class UNetUpBlock(nn.Module):
    def __init__(self, in_size, out_size, up_mode, padding, batch_norm , concat):
        super(UNetUpBlock, self).__init__()
        self.concat=concat
        if up_mode == 'upconv':
            self.up = nn.ConvTranspose2d(in_size, out_size, kernel_size=2,
                                         stride=2)
        elif up_mode == 'upsample':
            self.up = nn.Sequential(nn.Upsample(mode='bilinear', scale_factor=2),
                                    nn.Conv2d(in_size, out_size, kernel_size=1))

        self.conv_block = UNetConvBlock(in_size, out_size, padding, batch_norm)

    def center_crop(self, layer, target_size):
        _, _, layer_height, layer_width = layer.size()
        diff_y = (layer_height - target_size[0]) // 2
        diff_x = (layer_width - target_size[1]) // 2
        return layer[:, :, diff_y:(diff_y + target_size[0]), diff_x:(diff_x + target_size[1])]

    def forward(self, x, bridge):
        up = self.up(x)
        crop1 = self.center_crop(bridge, up.shape[2:])
        if(self.concat):
            out = torch.cat([up, crop1], 1)
        else:
            out = torch.cat([up, torch.rand_like(crop1) ], 1)
        out = self.conv_block(out)
        return out


# evaluation code
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device = 'cpu'
model = UNet()
model.load_state_dict(torch.load(model_path, map_location=torch.device(device)))
model = model.to(device)
model.eval()


# transform for input patch
img_transform = torchvision.transforms.Compose([
        torchvision.transforms.Resize((750, 750)),
        torchvision.transforms.ToTensor(),
        torchvision.transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])


patches = glob.glob(patches_dir + "*")
patches = patches[25000:30000]
for patch in patches:
    filename = patch.split("/")[-1]
    image = Image.open(patch).convert('RGB')
    image = img_transform(image)
    image = np.array(image)
    output_mask = np.zeros((750, 750))
    
    for index1 in range(0, image.shape[1]-128, 50):
        for index2 in range(0, image.shape[2]-128, 50):
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