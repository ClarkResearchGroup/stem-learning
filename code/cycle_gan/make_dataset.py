import os
from tifffile import imread, imwrite
import numpy as np

def process_image(input_file):
    """
       Takes in an image file and outputs a grayscale image. This grayscale is
       computed by flattening the channels potentially standardizing the distribution
       of pixel values
    """
    input_img = imread(input_file).astype(np.float32)
    if len(input_img.shape)==2:
        lx, ly = input_img.shape
        input_img = input_img.reshape((lx,ly,1))
    return input_img


def cut_data(data, fine_size, stride):
    """
    cuts up the data into pieces with dimension lx-by-ly
    data = 2-dimensional array with integer elements ranging from 0 to num_classes-1
    """
    data_shape = np.shape(data)
    (nx, ny) = (data_shape[0], data_shape[1])
    (sx, sy) = stride, stride
    (lx, ly) = fine_size, fine_size

    if lx > nx or ly > ny or sx > nx or sy > ny:
        print("Error: cut dimensions are bigger than the image")
        print(lx, ly)
        exit()

    return np.array([ data[i:i+lx,j:j+ly,:] for j in np.arange(0, ny - ly + 1, sy) \
            for i in np.arange(0, nx - lx + 1, sx)])


def find_min_max_vals(image_dir, num_channels=1):
    """
    returns the min and max pixel values of all images in a directory
    """
    fn_list = [x for x in os.listdir(image_dir) if ".tif" in x]
    min_v, max_v = num_channels*[999999], num_channels*[-999999]
    for fn in fn_list:
        img = process_image(os.path.join(image_dir, fn))
        for i in range(num_channels):
            upper, lower = np.max(img[:,:,i]), np.min(img[:,:,i])
            min_v[i] = lower if lower < min_v[i] else min_v[i]
            max_v[i] = upper if upper > max_v[i] else max_v[i]

    for i in range(num_channels):
        if min_v[i] == max_v[i]:
            min_v[i], max_v[i] = 0, 1

    return min_v, max_v

def parse_and_save_image(fn, image_dir, min_v, max_v, save_dir="./save/", fine_size=256, stride=256, num_channels=1, gaussian=0):
    """
    takes a large image and creates subimages of it and stores it in a directory
    """
    input_file = image_dir + fn
    data = process_image(input_file)
    for i in range(num_channels):
        x = data[:,:,i]
        if gaussian == 0:
            data[:,:,i] = 2*( (x - min_v[i])/(max_v[i] - min_v[i]) - .5)
        else:
            lx, ly = data[:,:,i].shape
            data[:,:,i] = 2*( (x - np.min(x))/(np.max(x) - np.min(x)) - .5)
            data[:,:,i] += np.random.normal(0, gaussian, size=(lx, ly))
            data[:,:,i] = 2*( (x - np.min(x))/(np.max(x) - np.min(x)) - .5)

    arr = cut_data(data, fine_size, stride)
    os.makedirs(save_dir, exist_ok=True)
    for i, img in enumerate(arr):
        imwrite("{}{}_{}.tiff".format(save_dir, fn[:fn.find('.')], str(i).zfill(3)), img)
    return

def parse_and_save_dir(image_dir, save_dir="./save/", fine_size=256, stride=256, num_channels=1, gaussian=0):
    """
    takes images in a directory and cuts them into subimages
    """
    fn_list = [x for x in os.listdir(image_dir) if ".tif" in x]
    min_v, max_v = find_min_max_vals(image_dir, num_channels)
    for fn in fn_list:
        parse_and_save_image(fn, image_dir, min_v, max_v, save_dir, fine_size, stride, num_channels, gaussian)
    return

def load_train_data(fn, num_channels):
    arr = process_image(fn)
    if np.random.random() > 0.5:
        arr = np.fliplr(arr)
    arr = np.rot90(arr, np.random.randint(4))
    return arr[:,:,:num_channels]