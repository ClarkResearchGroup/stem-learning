from make_data import *
import os
parent_dir = "/mnt/d/stem-learning/"
fn = "sim_pristine_0.1_noised"
input_dir = parent_dir + "data/WSe/data_for_gan/generated/GANNED_DATA_FOLDER/data_folder_{}/".format(fn)
#input_dir = parent_dir + "data/WSe/manual_simulation_211222/"
all_dirs = os.listdir(input_dir)
train_data_dirs = all_dirs[:-10]
test_data_dirs = all_dirs[-10:]
#label_list = ["1Doped", "1vacancy", "2Doped", "2vacancy", "metal_Doped", "metal_vacancy"]
label_list = ["1vacancy"]

lbl = label_list[0]

parsed_dir_name='parsed_label_{}'.format(lbl)
ftype = '.tiff'
l_shape = (256,256)
stride = (256,256)
one_pickle=False
tr_bs = 1000
ts_bs = 100

ones_percent = .00
tol = 0.25
show_plots=False

create_augments(input_dir, train_data_dirs, ftype)


for lbl in label_list:
    ll = [lbl]
    parsed_dir_name='parsed_label_{}'.format(lbl)

    make_data(input_dir, train_data_dirs, ll, l_shape, stride, ftype, \
            parsed_dir_name=parsed_dir_name, prefix="train", AUG=True, tol=tol, \
            ones_pcent=ones_percent, one_save=one_pickle, fsize=tr_bs)

    make_data(input_dir, test_data_dirs, ll, l_shape, stride, ftype, \
              parsed_dir_name=parsed_dir_name, prefix="test", AUG=False, tol=tol,\
              ones_pcent=ones_percent, one_save=True, fsize=ts_bs)

import numpy as np
import matplotlib.pyplot as plt
from make_data import *
parsed_fn = input_dir + parsed_dir_name + "/train_00000.p"
check_data(parsed_fn, l_shape=l_shape)
