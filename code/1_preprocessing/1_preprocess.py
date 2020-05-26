from make_data import *

input_dir = "../../data/ZnPc/"
data_dirs = ["1149_0", "1149_1", "1149_2", "1149_3", "1153"]
label_list = ["ZnPc"]
parsed_dir_name='parsed_label_ZnPc'
ftype = '.tif'

l_shape = (256,256)
stride = (64,64)
one_pickle=False
tr_bs = 2000
ts_bs = 200
ones_percent = .00
tol = 0.05
show_plots=False

create_augments(input_dir, data_dirs, ftype)

make_data(input_dir, label_list, data_dirs, l_shape, stride, ftype,\
        parsed_dir_name=parsed_dir_name, tr_bs=tr_bs, ts_bs=ts_bs, ones_percent=ones_percent, \
        tol=tol, show_plots=show_plots, one_save=one_pickle)
