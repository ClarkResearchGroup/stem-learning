import numpy as np
from make_dataset import cut_data, process_image
from tifffile import imwrite
import matplotlib.pyplot as plt
import os

def stitch(size_x, size_y, sx, sy, images):
    '''
    Modified stitch function for better performance while processing larger image
    '''
    def stitch_sum(size_x, size_y, sx, sy, images):
        (num_cuts, lx, ly, num_class) = np.shape(images)
        Nx = len(np.arange(0, size_x - lx + 1, sx))
        Ny = len(np.arange(0, size_y - ly + 1, sy))

        L_X = ((num_cuts-1)  % Nx)*sx + lx
        L_Y = ((num_cuts-1) // Nx)*sy + ly
        final_img = np.zeros((L_X, L_Y, num_class))
        for idx, img in enumerate(images):
            nx = idx % Nx
            ny = idx //Nx
            final_img[nx*sx:nx*sx + lx,  ny*sy:ny*sy + ly, :] += img      
        return final_img
    stitch_result = stitch_sum(size_x, size_y, sx, sy, images)/stitch_sum(size_x, size_y, sx, sy, np.ones_like(images))
    return stitch_result
    
def get_s(size_x, lx):
    '''
    computes the stride so that the whole input image
    is processed
    '''
    #for sx in range(lx, 0, -1):
    #    if (size_x - lx) % sx == 0:
    #        return sx
    return lx - 16

def get_avg_pred(model, cut):
    (lx, ly, _) = cut.shape
    rot_cuts  = np.array([ np.rot90(cut, k) for k in range(4) ])
    flip_cuts = np.array([ np.fliplr(r_cut) for r_cut in rot_cuts ])

    pred_rots  = model(rot_cuts)
    pred_flips = model(flip_cuts)

    pred_rots  = np.reshape(pred_rots , [4, lx, ly, pred_rots.shape[-1]])
    pred_flips = np.reshape(pred_flips, [4, lx, ly, pred_flips.shape[-1]])

    pred_1 = [ np.rot90(pred_rots[k], -k)             for k in range(4) ]
    pred_2 = [ np.rot90(np.fliplr(pred_flips[k]), -k) for k in range(4) ]

    preds = pred_1 + pred_2

    return sum(preds)/8.0

def generate_image(model, input_file, fine_size, stride=None, avg=False, \
        plot=False, save_data=False, save_dir='./',fname='generated_input.tiff'):
    ''' evauates an input image given the model'''

    print("processing data")
    input_img = process_image(input_file)
    (size_x, size_y, _) = input_img.shape

    print("loading model")
    stride = stride if stride != None else get_s(size_x, fine_size)
    print("stride: {}".format(stride))

    input_cuts = cut_data(input_img, fine_size, stride)
    #input_cuts = np.reshape(input_cuts, [-1, fine_size, fine_size, 1])
    num_cuts = len(input_cuts)

    print("predicting data")
    predictions = np.array([get_avg_pred(model,cut) for cut in input_cuts] if avg else model(input_cuts))
    predictions = np.reshape(np.array(predictions), [num_cuts, fine_size, fine_size, -1])

    print("stitching data")
    a = stitch(size_x, size_y, stride, stride, predictions)[:,:,0]
    
    if plot:
        plt.figure(figsize=(10,10))
        plt.imshow(a, cmap='gray')
        plt.axis('off')
        plt.show()

    if save_data:
        print("saving data")
        imwrite(save_dir + fname, a.astype(np.float32))

    return a

def GAN_image_folder(model, input_dir, save_dir, fine_size, stride=None, avg=False):
    """given a model, and input_directory of images, creates a new folder that GANs each image"""
    os.makedirs(save_dir, exist_ok=True)
    for input_file in os.listdir(input_dir):
        generate_image(model, os.path.join(input_dir, input_file), fine_size, stride, \
            avg, plot=False, save_data=True, save_dir=save_dir,fname=input_file)