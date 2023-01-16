from evaluate import *
import os
from sys import argv
import sys
import tensorflow as tf

parent_dir  = "/u/skhan/stem-learning/"

sys.path.insert(0, '{}code/cycle_gan'.format(parent_dir))
from models import unet_generator

idfn        = argv[1]
model_type  = argv[2]
lbl         = "1vacancy"

gaussian   = 0.2
lam_cycle  = 10
identifier = "20220708_MODEL_unet_ident_gen_fft_{}_SIM_pristine_gaussian_{}_EXP_{}".format(lam_cycle, gaussian, idfn)
#identifier  = "20220606_MODEL_unet_dist_gen_fft_10_SIM_pristine_gaussian_0.1_EXP_{}".format(idfn)

results_dir = "{}results/{}/results_{}/{}/".format(parent_dir, model_type, identifier, lbl)
save_dir    = "{}evals/".format(results_dir)
gan_cpt_dir = "{}cycle_gan_results/checkpoints/checkpoint_{}".format(parent_dir, identifier)



if model_type == "sim":
    generator_exp = unet_generator(1, 1, "instancenorm")
    generator_sim = unet_generator(1, 1, "instancenorm")

    def cycle(x):
        #return generator_sim(generator_exp(x))
        return generator_sim(x)

    try:
        ckpt = tf.train.Checkpoint(generator_exp=generator_exp, generator_sim=generator_sim)
        ckpt_manager = tf.train.CheckpointManager(ckpt, gan_cpt_dir, max_to_keep=3)
        if ckpt_manager.latest_checkpoint:
            cpath = ckpt_manager.latest_checkpoint
            ckpt.restore(cpath)
    except Exception as e:
        print("loading checkpoint failed")
        print(e)


def cycle_image(cycle, input_file, save_dir, fname, num_channels=1):
    print("GANNING {}".format(input_file))
    generate_image(cycle, input_file, (256,256), stride=(256,256), avg=False,\
            plot=False, save_data=True, save_dir=save_dir, fname=fname)


def evaluate_image(model_type, data_dir, results_dir, lbl, save_dir, save_fn_prefix):

    if not os.path.isdir(save_dir):
        os.makedirs(save_dir, exist_ok=True)

    model_fn         = results_dir + "model.json"
    model_weights_fn = results_dir + "weights.h5"
    input_file = data_dir + "input.tif"
    if model_type == "sim":
        cycle_image(cycle, input_file, save_dir, "{}_cycled.tif".format(save_fn_prefix))
        input_file = "{}{}_cycled.tif".format(save_dir, save_fn_prefix)

    l_shape = (256, 256)
    stride  = (64,64)
    plot = False
    fname = "evaluated_labels_{}.tiff".format(save_fn_prefix)

    print("EVALUATING {}".format(input_file))
    prediction = evaluate(model_fn, model_weights_fn, input_file, l_shape, stride,
                        avg=True, plot=plot, save_data=True, save_dir=save_dir, fname=fname)

    label_file_list = [data_dir + "label_{}.tif".format(lbl)]
    tol    = 0.5
    nconvs = 2
    r      = 9
    prefix = save_fn_prefix
    TP, FP, FN, TN, recall, precision, F1, bal_acc = calc_accuracy(prediction, label_file_list, tol=tol, bdy=32,
                                                                nconvs=nconvs, r=r, TN=0, plot=plot,
                                                                save_data=True, save_dir=save_dir,
                                                                prefix=prefix, verbose=False)
    return TP, FP, FN


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
    print("DATA ({}, {})".format(d, pic))
    data_dir = parent_dir + "data/WSe/test_experiment/{}_day/RR_{}/".format(d, pic)
    TP, FP, FN = evaluate_image(model_type, data_dir, results_dir, lbl, save_dir, "{}_{}".format(d, pic))
    write_to_file(",{},{},{}".format(TP,FP,FN))
