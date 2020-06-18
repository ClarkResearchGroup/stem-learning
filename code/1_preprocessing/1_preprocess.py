from make_data import *
from sys import argv

i = int(argv[1])

all_input_dirs = ["../../data/CuPcCl/original/", "../../data/CoMTPP/", "../../data/ZnPc/", "../../data/CuPcCl_CoMTPP/"]
all_data_dirs = [["0"], ["0"], ["1149_00"], ["combined"]]
all_label_list = [["CuPcCl"], ["CoMTPP"], ["ZnPc"], ["combined"]]


input_dir = all_input_dirs[i]
data_dirs = all_data_dirs[i]
label_list = all_label_list[i]
parsed_dir_name='parsed_label_{}'.format(label_list[0])
ftype = '.tif'

l_shape = (64,64)
stride = (32, 32)
one_pickle=True
tr_bs = 2000
ts_bs = 200
ones_percent = .00
tol = 0.05
show_plots=False

create_augments(input_dir, data_dirs, ftype)

make_data(input_dir, label_list, data_dirs, l_shape, stride, ftype,\
        parsed_dir_name=parsed_dir_name, tr_bs=tr_bs, ts_bs=ts_bs, ones_percent=ones_percent, \
        tol=tol, show_plots=show_plots, one_save=one_pickle)
