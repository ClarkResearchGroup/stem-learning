from input_data import *
from models import *
from training_utilities import *

parent_dir = '/home/abid/Dropbox/Development/programs/stem-learning/'
data_dir   = parent_dir + 'data/WSeTe/simulated/parsed_label_Se/'
sess_name  = 'Se'
N          = 256
k_fac      = 16
nb_classes = 2
num_steps  = 10

from os import makedirs
sess_dir = parent_dir + "results/" + sess_name + "/"
makedirs(sess_dir, exist_ok=True)

model_weights_fn = sess_dir + "weights.h5"
model_fn         = sess_dir + "model.json"
diagnostics_fn   = sess_dir + "diagnostics.dat"

model = construct_model(N, k_fac, nb_classes, sess_dir, model_fn, model_weights_fn)
step = setup_diagnostics(diagnostics_fn)

train(step, data_dir, N, nb_classes, model, diagnostics_fn, model_weights_fn, num_steps=num_steps)

