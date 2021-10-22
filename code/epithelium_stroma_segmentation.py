'''
Script provided by Haojia for running epi/stroma segmentation. Modified script for my use-case.
'''


# header files loaded
import os
os.environ["CUDA_VISIBLE_DEVICES"]="1"
from torch.autograd import Variable
import torchvision.transforms as transforms
import torch.utils.data as DATA
from PIL import Image
import numpy as np
import os
from skimage import io
import cv2
from epi_stroma_model import *
print("Header files loaded!")


# get the options selected by user
image_size = 1000
model_path = "model_files/latest_net_G.pth"
data_path = "results/patches"
output_path = "results/epithelium_stroma_masks/"
try:
    opt = TestOptions().parse()
except Exception as e:
    # add by haojia
    print('error accur here')
    exc_type, exc_value, exc_traceback = sys.exc_info()
    print("*** format_exc, first and last line:")
    formatted_lines = traceback.format_exc().splitlines()
    print(formatted_lines[0])
    print(formatted_lines[-1])
    print("*** format_exception:")
    print(repr(traceback.format_exception(exc_type, exc_value,
                                          exc_traceback)))
    print("*** extract_tb:")
    print(repr(traceback.extract_tb(exc_traceback)))
    print("*** format_tb:")
    print(repr(traceback.format_tb(exc_traceback)))
    print("*** tb_lineno:", exc_traceback.tb_lineno)


# convert PyTorch tensor to numpy
def tensor2im(image_tensor, imtype=np.uint8, normalize=True):
    image_numpy = image_tensor.cpu().float().numpy()
    if normalize:
        image_numpy = (np.transpose(image_numpy, (1, 2, 0)) + 1) / 2.0 * 255.0
    else:
        image_numpy = np.transpose(image_numpy, (1, 2, 0)) * 255.0      
    image_numpy = np.clip(image_numpy, 0, 255)
    if image_numpy.shape[2] == 1 or image_numpy.shape[2] > 3:        
        image_numpy = image_numpy[:,:,0]
    return image_numpy.astype(imtype)


# normalizing the data
data_transform = transforms.Compose([transforms.ToTensor(), transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5)),])    
DATASET = CustomDataLoader(root=data_path, transform=data_transform)
data_loader = DATA.DataLoader(dataset=DATASET, batch_size=1, shuffle=False, num_workers=0)
print('Finished loading dataset...')


# load model
netG = networks.define_G() 
networks.load_network(network=netG, save_path=model_path)
#netG.cuda()
print('Loaded model...')


# main loop for testing
for i, (data, img_path) in enumerate(data_loader):
    #input_concat = Variable(data, volatile=True).cuda()
    input_concat = Variable(data, volatile=True)
    fake_image = netG.forward(input_concat)
    
    image_pil = tensor2im(fake_image.data[0])  
    file_name = os.path.basename("".join(img_path))
    ret, thresh_eimg = cv2.threshold(image_pil[:,:,0], 1, 255, cv2.THRESH_BINARY)
    ret, thresh_simg = cv2.threshold(image_pil[:,:,1], 1, 255, cv2.THRESH_BINARY)
    eimg = Image.fromarray(thresh_eimg).resize((image_size, image_size), Image.BICUBIC)
    simg = Image.fromarray(thresh_simg).resize((image_size, image_size), Image.BICUBIC)
    image = np.array(eimg)
    image_inv = cv2.threshold(image, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]

    # filter using contour area and remove small noise
    cnts = cv2.findContours(image_inv, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    cnts = cnts[0] if len(cnts) == 2 else cnts[1]
    for c in cnts:
        area = cv2.contourArea(c)
        if area < 50:
            cv2.drawContours(image_inv, [c], -1, (0, 0, 0), -1)

    # filter using contour area and remove small noise
    output_mask = 255 - image_inv
    cnts = cv2.findContours(output_mask, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    cnts = cnts[0] if len(cnts) == 2 else cnts[1]
    for c in cnts:
        area = cv2.contourArea(c)
        if area < 50:
            cv2.drawContours(output_mask, [c], -1, (0, 0, 0), -1)

    cv2.imwrite(output_path + file_name[:-4] + '_epi_stroma_mask.png', output_mask)
    print('process image... %s' % img_path)