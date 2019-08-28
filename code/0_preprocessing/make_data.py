from image_parse import *
import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from random import shuffle
import os, subprocess, pickle

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
        command = "./generate_training_set.sh " + augdir + " " + ftype[1:]
        subprocess.call(command, shell=True)


def _pickle_data(info, parsed_dir, data):
    '''
    helper function for pickle_data. data is a list of
    (input, label) data tuples.
    '''
    (start, end, tr_bs, p_name) = info
    train_data = data[start:start + tr_bs]
    test_data = data[start + tr_bs:end]

    save_data(train_data, parsed_dir + "train/train_" + p_name + ".p")
    save_data(test_data,  parsed_dir + "test/test_"  + p_name + ".p")


def make_data(input_dir, label_list, data_dirs, l_shape, stride, ftype, parsed_dir_name='parsed',\
        tr_bs=100, ts_bs=10, ones_percent=0, tol=1e-5, show_plots=False, one_pickle=False):
    '''
    This function creates pickled data ready for training and testing. The input directory
    must contain a folder called augments in which contains folders of transformations of the
    original image data. each folder in augments must contain an image called input.ftype, and
    one or more label images called label_X.ftype, where X is an element in label_list
    data_dirs is the set of raw images that are desired for training.
    The output is a folder called parsed which contains two folders, train and test. In these
    folders are pickled data of each of the augmented images desired for training.
    '''

    (lx, ly) = l_shape
    (sx, sy) = stride

    # create a directory to store pickles
    parsed_dir = input_dir + parsed_dir_name + "/"
    if not os.path.isdir(parsed_dir):
        os.mkdir(parsed_dir)
        os.mkdir(parsed_dir + "train/")
        os.mkdir(parsed_dir + "test/")



    #go through each image
    data = []
    for f in data_dirs:
        augments = input_dir + f + "/augments/"
        if f == 'parsed':
            continue
        print(f)

        #go through each augment of the image
        for aug_dir in os.listdir(augments):
            print(aug_dir)
            full_aug_dir =  augments + aug_dir
            aug_name = f + "_" + aug_dir

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

            if one_pickle:
                data += zip(input_cuts, label_cuts)

            else:
                data = zip(input_cuts, label_cuts)
                shuffle(data)


                # save the cut data into picked data files
                num_examples = len(input_cuts)
                ex_per_p = tr_bs + ts_bs
                num_pickles = int(num_examples/ex_per_p)
                last_bs = num_examples%ex_per_p
                include_last = (last_bs > .1*tr_bs)
                print(str(num_examples)+"examples and "+str(num_pickles+int(include_last))+" files")

                data_idx_list = [(i*ex_per_p, (i+1)*ex_per_p, tr_bs, aug_name + str(i).zfill(3)) \
                        for i in range(num_pickles)]
                if include_last:
                    data_idx_list.append((num_pickles*ex_per_p, num_pickles*ex_per_p + last_bs,\
                            last_bs - int(last_bs/10), aug_name + str(num_pickles).zfill(3)))

                [_pickle_data(info, parsed_dir, data) for info in  data_idx_list]

        if one_pickle:
            shuffle(data)
            num_ex = len(data)
            print(str(num_ex) + " total examples")
            info = (0, num_ex, int(.9*num_ex), f)
            _pickle_data(info, parsed_dir, data)
            data = []

def check_data(parsed_fn, idx=-1, l_shape=(128,128)):
    '''
    this script looks at data ready to be fed into the CNN.
    parsed_fn = the path to the pickled data
    idx       = optional index of pickled data
    l_shape   = size of the training images
    '''

    data = pickle.load(open(parsed_fn, 'rb'))
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






