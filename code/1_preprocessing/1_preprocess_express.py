from make_data import *
input_dir = "../../data/WSe/sim_50/"
train_data_dirs = ["image{}".format(i) for i in range(50) if i%5 != 0]
test_data_dirs  = ["image{}".format(i) for i in range(50) if i%5 == 0]
label_list = ["1Doped", "1vacancy", "2Doped", "2vacancy", "metal_Doped", "metal_vacancy"]

parsed_dir_name='parsed_label_{}'.format(label_list[0])
ftype = '.tiff'
l_shape = (256,256)
stride = (64,64)
one_pickle=False
tr_bs = 1000
ts_bs = 100

ones_percent = .00
tol = 0.5
show_plots=False


#create_augments(input_dir, train_data_dirs, ftype)


for lbl in label_list:
    ll = [lbl]
    parsed_dir_name='parsed_label_{}'.format(lbl)

    #make_data(input_dir, train_data_dirs, ll, l_shape, stride, ftype, \
    #        parsed_dir_name=parsed_dir_name, prefix="train", AUG=True, tol=tol, \
    #        ones_pcent=ones_percent, one_save=one_pickle, fsize=tr_bs)
    make_data(input_dir, test_data_dirs, ll, l_shape, stride, ftype, \
              parsed_dir_name=parsed_dir_name, prefix="test", AUG=False, tol=tol,\
              ones_pcent=ones_percent, one_save=True, fsize=ts_bs)


#import numpy as np
#from make_data import *
#parsed_fn = input_dir + parsed_dir_name + "/train_00000.p"
#check_data(parsed_fn, l_shape=l_shape)
