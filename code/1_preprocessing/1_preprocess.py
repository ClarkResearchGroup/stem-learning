from make_data import *
import numpy as np

input_dir = "../../data/WSeTe/simulated/"
train_dirs = ["0", "1", "2"]
test_dirs  = ["3"]
label_list = ["Se"]
parsed_dir_name='parsed_label_Se'
ftype = '.tiff'

l_shape = (256,256)
stride = (64,64)

one_pickle=False
tr_fsize = 2000
ts_fsize = 200

tol = 0.05

ones_percent = 0.

create_augments(input_dir, train_dirs, ftype)

make_data(input_dir, train_dirs, label_list, l_shape, stride, ftype parsed_dir_name=parsed_dir_name, \
        prefix="train", AUG=True, tol=tol, ones_pcent=ones_percent, one_save=one_pickle, fsize=tr_fsize)
make_data(input_dir, test_dirs, label_list, l_shape, stride, ftype, parsed_dir_name=parsed_dir_name, \
        prefix="test", AUG=False, tol=tol, ones_pcent=ones_percent, one_save=one_pickle, fsize=ts_fsize)

parsed_fn = input_dir + parsed_dir_name + "/test_00000.p"
check_data(parsed_fn, l_shape=l_shape)

