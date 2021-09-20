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

    img_shape = len(input_img.shape)
    assert  img_shape == 2, "improper image shape of " + str(img_shape)
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

    return np.array([ data[i:i+lx, j:j+ly] for j in np.arange(0, ny - ly + 1, sy) \
            for i in np.arange(0, nx - lx + 1, sx)])


def find_min_max_vals(image_dir):
    """
    returns the min and max pixel values of all images in a directory
    """
    fn_list = [x for x in os.listdir(image_dir) if ".tif" in x]
    min_v, max_v = 999999, -999999
    for fn in fn_list:
        img = imread(os.path.join(image_dir, fn))
        upper, lower = np.max(img), np.min(img)
        min_v = lower if lower < min_v else min_v
        max_v = upper if upper > max_v else max_v
    return min_v, max_v


def parse_and_save_image(fn, image_dir, min_v, max_v, save_dir="./save/", fine_size=256, stride=256):
    """
    takes a large image and creates subimages of it and stores it in a directory
    """
    input_file = image_dir + fn
    data = process_image(input_file)
    arr = cut_data(data, fine_size, stride)

    arr = 2*( (arr - min_v)/(max_v - min_v) - .5)
    N = len(arr)

    os.makedirs(save_dir, exist_ok=True)
    for i, img in enumerate(arr):
        imwrite("{}{}_{}.tiff".format(save_dir, fn[:-4], str(i).zfill(3)), img)
    return

def parse_and_save_dir(image_dir, save_dir="./save/", fine_size=256, stride=256):
    """
    takes images in a directory and cuts them into subimages
    """
    fn_list = [x for x in os.listdir(image_dir) if ".tif" in x]
    min_v, max_v = find_min_max_vals(image_dir)
    for fn in fn_list:
        parse_and_save_image(fn, image_dir, min_v, max_v, save_dir, fine_size, stride)
    return

def load_train_data(fn):
    arr = np.array(imread(fn)).astype(np.float32)
    if np.random.random() > 0.5:
        arr = np.fliplr(arr)
    arr = np.rot90(arr, np.random.randint(4))
    (lx, ly) = arr.shape
    return arr.reshape([lx,ly,1])