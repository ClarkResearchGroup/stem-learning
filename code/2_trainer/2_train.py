#!/usr/bin/env python
# coding: utf-8

from input_data import *
from models import *
from training_utilities import *
from sys import argv

i       = int(argv[1])
k_fac   = int(argv[2])
dropout = float(argv[3])

all_input_dirs = ["data/CuPcCl/original/", "data/CoMTPP/", "data/ZnPc/", "data/CuPcCl_CoMTPP/"]
all_label_list = [["CuPcCl"], ["CoMTPP"], ["ZnPc"], ["combined"]]
label = all_label_list[i][0]

parent_dir = '/home/skhan/stem-learning/'
data_dir   = '{}{}parsed_label_{}/'.format(parent_dir, all_input_dirs[i], label)
sess_name  = label
N          = 64
#k_fac      = 1
nb_classes = 2
#dropout    = 0.35





from os import makedirs
sess_dir = parent_dir + "results_dropout_{}_k_{}/".format(dropout, k_fac) + sess_name + "/"
makedirs(sess_dir, exist_ok=True)

model_weights_fn = sess_dir + "weights.h5"
model_fn         = sess_dir + "model.json"
diagnostics_fn   = sess_dir + "diagnostics.dat"

model = construct_model(N, k_fac, nb_classes, sess_dir, model_fn, model_weights_fn, dropout=dropout)
step = setup_diagnostics(diagnostics_fn)

train(step, data_dir, N, nb_classes, model, diagnostics_fn, model_weights_fn)

