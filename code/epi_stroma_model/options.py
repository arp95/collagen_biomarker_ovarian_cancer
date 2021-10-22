#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 16:02:15 2018
@author: zzl
Modified by Arpit for a particular use-case
"""
import argparse


class TestOptions():
    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.initialized = False

    def initialize(self):        
        self.parser.add_argument('--trained_model', type=str, default='./latest_net_G.pth', help='path to trained model')        
        self.parser.add_argument('--image_size', type=int, default=1000, help='resize images to this size')
        self.parser.add_argument('--dataroot', type=str, default='./data/') 
        self.parser.add_argument('--results_dir', type=str, default='./results/', help='saves results here.')
        self.initialized = True

    # Default use single gpu0
    def parse(self, save=True):        
        if not self.initialized:
            self.initialize()            
        self.opt = self.parser.parse_args()
        args = vars(self.opt)
        print('------------ Options -------------')
        for k, v in sorted(args.items()):
            print('%s: %s' % (str(k), str(v)))
        print('-------------- End ----------------')
        return self.opt