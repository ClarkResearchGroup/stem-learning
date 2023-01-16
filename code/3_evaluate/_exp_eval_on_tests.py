from evaluate import *
import os
from sys import argv
sys.path.insert(0, '../cycle_gan')
from generate_image import generate_image
#from models import unet_gen

parent_dir  = "/u/skhan/stem-learning/"
idfn        = argv[1]
model_type  = argv[2]
lbl         = "1vacancy"

identifier  = "20220606_MODEL_unet_dist_gen_fft_10_SIM_pristine_gaussian_0.1_EXP_{}".format(idfn)
results_dir = "{}results/{}/results_{}/{}/".format(parent_dir, model_type, identifier, lbl)
save_dir    = "{}evals/".format(results_dir)

def evaluate_image(data_dir, results_dir, lbl, save_dir, save_fn_prefix):


    model_fn         = results_dir + "model.json"
    model_weights_fn = results_dir + "weights.h5"
    input_file       = data_dir + "input.tif"

    l_shape = (256, 256)
    stride  = (64,64)
    plot = False
    fname = "evaluated_labels_{}.tiff".format(save_fn_prefix)

    prediction = evaluate(model_fn, model_weights_fn, input_file, l_shape, stride,
                        avg=True, plot=plot, save_data=True, save_dir=save_dir, fname=fname)

    label_file_list = [data_dir + "label_{}.tif".format(lbl)]
    tol    = 0.5
    nconvs = 2
    r      = 9
    prefix = save_fn_prefix
    TP, FP, FN, TN, recall, precision, F1, bal_acc = calc_accuracy(prediction,
                                                            label_file_list, tol=tol, bdy=32,
                                                            nconvs=nconvs, r=r, TN=0,plot=plot,
                                                            save_data=True, save_dir=save_dir,
                                                            prefix=prefix, verbose=False)
    return TP, FP, FN


#def make_cycle_image(data_dir, save_dir):




if not os.path.isdir(save_dir):
        os.makedirs(save_dir, exist_ok=True)

label_csv = "{}label_row.csv".format(save_dir)
def write_to_file(vals):
    with open(label_csv, "a") as f:
        f.write(vals)

if os.path.exists(label_csv):
    os.remove(label_csv)

write_to_file("{}_{}".format(idfn, model_type))

for (d, pic) in [(211,1528), (211, 1922), (211, 2124), (107, 1750), (107, 1840), (107, 1847)]:
    data_dir = parent_dir + "data/WSe/test_experiment/{}_day/RR_{}/".format(d, pic)
    TP, FP, FN = evaluate_image(data_dir, results_dir, lbl, save_dir, "{}_{}".format(d, pic))
    write_to_file(",{},{},{}".format(TP,FP,FN))
