from __future__ import division
import numpy as np
from imageio import imread, imwrite

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


def flatten_channels(img):
    '''
    Takes an a 3D image array of length x width x num_channels and returns a 2D image
    of length x width. The channels are summed and the overall image is normed.
    '''
    if len(img.shape) == 3:
        c = min(3, img.shape[-1]) # picks out the rgb channels from rgba
        return np.dot(img[...,:c], np.ones(c))


def process_image(input_file, standardize=False):
    """
       Takes in an image file and outputs a grayscale image. This grayscale is
       computed by flattening the channels potentially standardizing the distribution
       of pixel values
    """

    input_img = imread(input_file).astype(np.float32)

    if len(input_img.shape) == 3:
        input_img = flatten_channels(input_img)

    if standardize:
        input_img = (input_img - np.mean(input_img))/np.std(input_img)

    img_shape = len(input_img.shape)
    assert  img_shape == 2, "improper image shape of " + str(img_shape)
    return input_img


def encode(label_img, nb_classes):
    """
    takes an image and replaces each pixel with an array of size nb_classes, where
    the array is a zero-filled except at the location of the value of the pixel. For
    example, if nb_classes = 3, then the pixel values of the images are either 0, 1,
    or 2, and then
    0 -> [0,0,1]
    1 -> [0,1,0]
    2 -> [1,0,0]
    """
    (nx, ny) = label_img.shape
    flattened = np.reshape(label_img, -1)
    encoded = np.eye(nb_classes)[flattened]
    return np.reshape(encoded, (nx, ny, -1))


def process_label(label_file_list, tol=1.e-5):
    '''
    creates a 2 dimensional array of integers ranging from 0 to len(label_files)
    label_files is an array of file names for different labels of the same input image
    '''
    lbl = imread(label_file_list[0]).astype(np.float64)
    label_shape = lbl.shape
    (nx, ny) = (label_shape[0], label_shape[1])
    num_classes = len(label_file_list) + 1
    label_data = np.zeros((nx, ny, num_classes))

    for i, label_file in enumerate(label_file_list):
        img   = process_image(label_file, standardize=False)
        img   = (img - np.min(img))/np.ptp(img)
        label_data[:,:,i+1] = ((img > tol).astype(int))

    label_data[:,:,0] = (np.sum(label_data[:,:,1:], axis=2) == 0).astype(int)
    return label_data


def cut_data(data, lx, ly, stride=(1, 1), standardize=False):
    """
    cuts up the data into pieces with dimension lx-by-ly
    data = 2-dimensional array with integer elements ranging from 0 to num_classes-1
    """
    data_shape = np.shape(data)
    (nx, ny) = (data_shape[0], data_shape[1])
    (sx, sy) = stride

    if lx > nx or ly > ny or sx > nx or sy > ny:
        print("Error: cut dimensions are bigger than the image")
        print(lx, ly)
        exit()

    cut_data = [data[i:i+lx, j:j+ly] for j in np.arange(0, ny - ly + 1, sy) \
                                     for i in np.arange(0, nx - lx + 1, sx)]

    if standardize:
        cut_data = [(x - np.mean(x))/np.std(x) for x in cut_data]

    return np.array(cut_data)



def diff_images(label_input_fn, input_fn, save_fn, tol=1e-5):
    '''
    Given an image, and a copy of it with labels overlayed on top, this function
    creates an image of the labels without the background, and returns the
    image
    '''
    label_input_img = Image.open(label_input_fn)
    input_img       = Image.open(input_fn)

    label_img = np.array((label_input_img - input_img)).astype(np.float64)
    label_img = flatten_channels(label_img)
    label_img = (label_img > tol).astype(int)
    imwrite(save_fn, label_img)
    return label_img

def sift_cuts(input_cuts, label_cuts, ones_percent):
    '''
    given an array of input images and label images, returns a new list of input and label images
    that have a proportion of ones of at least ones_percent in the label images.
    '''
    ones_percent /=100

    lbl_shape = label_cuts[0].shape
    tot_pixels = float(lbl_shape[0]*lbl_shape[1]*(lbl_shape[2]-1))
    new_input_cuts = []
    new_label_cuts = []
    for idx, lbl in enumerate(label_cuts):
        lbl_percent = np.sum(lbl[:,:,1:])/tot_pixels
        if lbl_percent >= ones_percent:
            new_input_cuts.append(input_cuts[idx])
            new_label_cuts.append(label_cuts[idx])

    return new_input_cuts, new_label_cuts


def get_classifier_label(label_img, ones_per_dot=75):
    """takes in a label image and returns 1 if there is a
       dopant
    """
    return np.sum(label_img) >= ones_per_dot

def save_data(data, file_name):
    with open(file_name, 'wb') as output:
        pickle.dump(data, output, pickle.HIGHEST_PROTOCOL)

def load_data(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)
