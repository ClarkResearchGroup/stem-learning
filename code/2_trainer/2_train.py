from input_data import *
from models import *
from training_utilities import *
from sys import argv


parent_dir = "/u/skhan/stem-learning/"
data_path  = "data/WSe/data_for_gan/generated/GANNED_DATA_FOLDER/"
# pick model type
model_type = argv[2] #sim or exp

# locate checkpoint folder
idfn       = argv[1]
gaussian   = 0.2
lam_cycle  = 10
#identifier = "20220606_MODEL_unet_dist_gen_fft_10_SIM_pristine_gaussian_0.1_EXP_{}".format(idfn)
identifier = "20220708_MODEL_unet_ident_gen_fft_{}_SIM_pristine_gaussian_{}_EXP_{}".format(lam_cycle, gaussian, idfn)

sess_name  = '1vacancy'
data_dir  = "{}{}{}/data_folder_{}/parsed_label_{}/".format(parent_dir, data_path, model_type, identifier, sess_name)

N          = 256
k_fac      = 16
nb_classes = 2
num_steps  = 500



from os import makedirs
sess_dir = "{}results/{}/results_{}/{}/".format(parent_dir, model_type, identifier, sess_name)
makedirs(sess_dir, exist_ok=True)

model_weights_fn = sess_dir + "weights.h5"
model_fn         = sess_dir + "model.json"
diagnostics_fn   = sess_dir + "diagnostics.dat"



model = construct_model(N, k_fac, nb_classes, sess_dir, model_fn, model_weights_fn)
step = setup_diagnostics(diagnostics_fn)



train(step, data_dir, N, nb_classes, model, diagnostics_fn, model_weights_fn, num_steps=num_steps)

