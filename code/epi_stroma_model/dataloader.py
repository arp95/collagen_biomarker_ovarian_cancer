#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 11:04:44 2018
@author: zzl
"""
import torch.utils.data as data
from PIL import Image
import os


IMG_EXTENSIONS = ['.jpg', '.JPG', '.jpeg', '.JPEG', '.png', '.PNG', '.ppm', '.PPM', '.bmp', '.BMP', '.tiff']

# check if it is a valid image
def is_image_file(filename):
    return any(filename.endswith(extension) for extension in IMG_EXTENSIONS)

# get the images in the data root directory
def make_dataset(dir):
    images = []
    assert os.path.isdir(dir), '%s is not a valid directory' % dir

    for root, _, fnames in sorted(os.walk(dir)):
        for fname in fnames:
            if is_image_file(fname):
                path = os.path.join(root, fname)
                images.append(path)
    return images


# main class
class CustomDataLoader(data.Dataset):
      def __init__(self, root, transform=None, return_paths=True):
          imgs = make_dataset(root)
          if len(imgs) == 0:
              raise(RuntimeError("Found 0 images in: " + root + "\n"
                               "Supported image extensions are: " +
                               ",".join(IMG_EXTENSIONS)))
          self.root = root
          self.imgs = imgs
          self.transform = transform
          self.return_paths = return_paths

      def __getitem__(self, index):
          path = self.imgs[index]
          img = Image.open(path).convert('RGB')
          if self.transform is not None:
              img = self.transform(img)
          if self.return_paths:
              return img, path
          else:
              return img

      def __len__(self):
          return len(self.imgs)