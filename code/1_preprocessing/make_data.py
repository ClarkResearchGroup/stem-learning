from image_parse import *
import numpy as np
import matplotlib.pyplot as plt
from random import shuffle
from generate_training_set import make_augments
import os

def diff_labels(input_dir, label_list, data_dirs, ftype, tol=1e-5):
    '''
    This function is used when the labeled images have the input image in the background.
    It renames the labeled image and
    '''
    for data_dir in data_dirs:
        for label in label_list:
            label_fn       = input_dir + data_dir + "/label_" + label + ftype
            input_fn       = input_dir + data_dir + "/input" + ftype
            label_input_fn = input_dir + data_dir + "/label_input_" + label + ftype
            os.rename(label_fn, label_input_fn)
            diff_images(label_input_fn, input_fn, label_fn, tol)


def create_augments(input_dir, data_dirs, ftype):
    '''
    This function calls the bash script 'generate_training_set.sh'
    '''
    for data_dir in data_dirs:
        print("creating augments in " + data_dir)
        augdir = input_dir + data_dir
        make_augments(augdir, ftype)


def _save_data(info, parsed_dir, data):
    '''
    helper function for save_data. data is a list of
    (input, label) data tuples.
    '''
    (start, end, tr_bs, p_name) = info
    train_data = data[start:start + tr_bs]
    test_data = data[start + tr_bs:end]

    save_data(train_data, parsed_dir + "train/train_" + p_name + ".p")
    save_data(test_data,  parsed_dir + "test/test_"  + p_name + ".p")


def make_data(input_dir, label_list, data_dirs, l_shape, stride, ftype, parsed_dir_name='parsed',\
        tr_bs=100, ts_bs=10, ones_percent=0, tol=1e-5, show_plots=False, one_save=False):
    '''
    This function creates npy data ready for training and testing. The input directory
    must contain a folder called augments in which contains folders of transformations of the
    original image data. each folder in augments must contain an image called input.ftype, and
    one or more label images called label_X.ftype, where X is an element in label_list
    data_dirs is the set of raw images that are desired for training.
    The output is a folder called parsed which contains two folders, train and test. In these
    folders are npy data of each of the augmented images desired for training.
    '''

    (lx, ly) = l_shape
    (sx, sy) = stride

    # create a directory to store npys
    parsed_dir = input_dir + parsed_dir_name + "/"
    if not os.path.isdir(parsed_dir):
        os.mkdir(parsed_dir)
        os.mkdir(parsed_dir + "train/")
        os.mkdir(parsed_dir + "test/")



    #go through each image
    data = []
    i = 0
    for f in data_dirs:
        augments = input_dir + f + "/augments/"
        print(f)

        #go through each augment of the image
        for aug_dir in os.listdir(augments):
            print(aug_dir)
            full_aug_dir =  augments + aug_dir

            input_file = full_aug_dir + "/input" + ftype
            label_files = [ full_aug_dir + "/label_" + label + ftype for label in label_list]

            input_img = process_image(input_file, standardize=True)
            label_img = process_label(label_files, tol=tol)

            if show_plots:
                plt.imshow(input_img)
                for c in range(len(label_files) + 1):
                    plt.figure()
                    plt.imshow(label_img[:,:,c])
                plt.show()
                continue

            input_cuts = cut_data(input_img, lx, ly, stride=(sx, sy))
            label_cuts = cut_data(label_img, lx, ly, stride=(sx, sy))

            # only allow labels that have more than a certain proportion of ones
            (input_cuts, label_cuts) = sift_cuts(input_cuts, label_cuts, ones_percent)

            data += list(zip(input_cuts, label_cuts))


            if not one_save and len(data) >= tr_bs + ts_bs:
                # save the cut data into picked data files
                print("saving file {}".format(i))
                shuffle(data)
                info = (0, tr_bs + ts_bs, tr_bs, str(i).zfill(5))
                _save_data(info, parsed_dir, data)
                del data[:(tr_bs + ts_bs)]
                i += 1

        if one_save:
            shuffle(data)
            num_ex = len(data)
            print(str(num_ex) + " total examples")
            info = (0, num_ex, int(.9*num_ex), f)
            _save_data(info, parsed_dir, data)
            data = []



def check_data(parsed_fn, idx=-1, l_shape=(128,128)):
    '''
    this script looks at data ready to be fed into the CNN.
    parsed_fn = the path to the npy data
    idx       = optional index of npy data
    l_shape   = size of the training images
    '''

    data = load_data(parsed_fn)
    idx = np.random.randint(len(data)) if idx == -1 else idx
    [img, lbl] = data[idx]
    nb_classes = lbl.shape[-1]
    lbl = np.reshape(lbl, [l_shape[0], l_shape[1], nb_classes])

    plt.figure()
    plt.imshow(img)
    for c in range(nb_classes):
        plt.figure()
        plt.imshow(img)
        plt.imshow(lbl[:,:,c], alpha=0.5)
    plt.show()






