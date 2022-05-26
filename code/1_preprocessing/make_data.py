from image_parse import *
import numpy as np
import matplotlib.pyplot as plt
from random import shuffle
from generate_augments import make_augments
import os


def create_augments(input_dir, data_dirs, ftype):
    '''
    This function calls the bash script 'generate_training_set.sh'
    '''
    for data_dir in data_dirs:
        print("creating augments in " + data_dir)
        augdir = input_dir + data_dir
        make_augments(augdir, ftype)


def get_image_arr(input_file, label_files, tol, lx, ly, stride, ones_pcent):

    input_img = process_image(input_file, standardize=False)
    label_img = process_label(label_files, tol=tol)

    input_cuts = cut_data(input_img, lx, ly, stride=stride, standardize=True)
    label_cuts = cut_data(label_img, lx, ly, stride=stride, standardize=False)

    # only allow labels that have more than a certain proportion of ones
    if ones_pcent >0:
        (input_cuts, label_cuts) = sift_cuts(input_cuts, label_cuts, ones_pcent)

    return list(zip(input_cuts, label_cuts))


def make_data(input_dir, data_dirs, label_list, l_shape, stride, ftype, parsed_dir_name="parsed", prefix="train",\
        AUG=True, tol=1e-5, ones_pcent=0, one_save=True, fsize=1000, ret_only=False):

    (lx, ly) = l_shape

    parsed_dir = input_dir + parsed_dir_name + "/"
    if not os.path.isdir(parsed_dir):
        if not ret_only:
            os.mkdir(parsed_dir)

    data, i = [], 0

    for f in data_dirs:
        print("parsing directory ", f)
        #go through each augment of the image
        if AUG:
            augments = input_dir + f + "/augments/"
            for aug_dir in os.listdir(augments):
                full_aug_dir =  augments + aug_dir

                input_file  =   full_aug_dir + "/input" + ftype
                label_files = [ full_aug_dir + "/label_" + label + ftype for label in label_list]
                data += get_image_arr(input_file, label_files, tol, lx, ly, stride, ones_pcent)

                if not one_save and len(data) >= fsize:
                    shuffle(data)
                    while(len(data) >= fsize):
                        print("saving file {}".format(i))
                        save_data(data[:fsize], parsed_dir + prefix + "_" +  str(i).zfill(5) + ".p")
                        data = data[fsize:]
                        i += 1

        else:
            input_file  =   input_dir + f + "/input" + ftype
            label_files = [ input_dir + f + "/label_" + label + ftype for label in label_list]
            data += get_image_arr(input_file, label_files, tol, lx, ly, stride, ones_pcent)


    shuffle(data)
    if one_save and not ret_only:
        print("saving file {} with {} examples".format(i, len(data)))
        save_data(data, parsed_dir + prefix + "_" +  str(i).zfill(5) + ".p")
    return data


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
    plt.imshow(img, cmap='gray')
    for c in range(nb_classes):
        plt.figure()
        plt.imshow(img, cmap='gray')
        plt.imshow(lbl[:,:,c], alpha=0.5, cmap='gray')
    plt.show()






