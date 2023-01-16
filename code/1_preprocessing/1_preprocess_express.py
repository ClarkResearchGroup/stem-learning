from make_data import *
import os
from sys import argv

parent_dir = "/u/skhan/stem-learning/"
data_dir   = parent_dir  + "data/WSe/data_for_gan/"
# pick model type
model_type = argv[2] #sim or exp

# locate checkpoint folder
idfn       = argv[1]
#identifier = "20220606_MODEL_unet_dist_gen_fft_10_SIM_pristine_gaussian_0.1_EXP_{}".format(idfn)
lam_cycle  = 10
gaussian   = 0.2
identifier = "20220708_MODEL_unet_ident_gen_fft_{}_SIM_pristine_gaussian_{}_EXP_{}".format(lam_cycle, gaussian, idfn)


input_dir  = data_dir + "generated/GANNED_DATA_FOLDER/{}/data_folder_{}/".format(model_type, identifier)

all_dirs = [x for x in os.listdir(input_dir) if "parsed_label" not in x]
train_data_dirs = all_dirs[:-10]
test_data_dirs = all_dirs[-10:]
lbl = "1vacancy"

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


ll = [lbl]
parsed_dir_name='parsed_label_{}'.format(lbl)

make_data(input_dir, train_data_dirs, ll, l_shape, stride, ftype, \
        parsed_dir_name=parsed_dir_name, prefix="train", AUG=True, tol=tol, \
        ones_pcent=ones_percent, one_save=one_pickle, fsize=tr_bs)

make_data(input_dir, test_data_dirs, ll, l_shape, stride, ftype, \
          parsed_dir_name=parsed_dir_name, prefix="test", AUG=False, tol=tol,\
          ones_pcent=ones_percent, one_save=True, fsize=ts_bs)

"""
import numpy as np
import matplotlib.pyplot as plt
from make_data import *
parsed_fn = input_dir + parsed_dir_name + "/train_00000.p"
check_data(parsed_fn, l_shape=l_shape)
"""


