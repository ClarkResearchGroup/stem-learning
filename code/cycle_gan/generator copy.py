import tensorflow as tf
from models import unet_generator
from generate_image import generate_image, GAN_image_folder
import os, shutil
from sys import argv


# locate checkpoint folder
idfn       = argv[1]
model_type = argv[2] #sim or exp
lam_cycle  = 10
exp_dir_fn = argv[3]
gaussian   = 0.2 if len(argv) < 4 else float(argv[4])

parent_dir = "/u/skhan/stem-learning/"
#data_dir   = "{}data/WSe/data_for_gan/".format(parent_dir)
data_dir   = "{}data/{}/{}_{}/".format(parent_dir, exp_dir_fn, exp_dir_fn, idfn)
# pick model type

identifier = "20220727_MODEL_unet_ident_gen_fft_{}_SIM_pristine_gaussian_{}_EXP_{}_{}".format(lam_cycle, gaussian, exp_dir_fn, idfn)

# image folder to GAN
#input_dir = "{}simulation/sim_pristine_{}_noised/".format(data_dir, gaussian)
input_dir               = "{}sim_1K_gauss_{}_256_slices/".format(data_dir, gaussian)
GANNED_dir              = "{}GANNED_SIM/{}/sim_abberation_{}/".format(data_dir, model_type, identifier)
# folder location to fofor FCN data
data_folder_with_labels = "{}simulation/data_folders_with_labels/".format(data_dir)
target_folder           = "{}generated/GANNED_DATA_FOLDER/{}/data_folder_{}".format(data_dir, model_type, identifier)

print(identifier)
print(input_dir)
print(GANNED_dir)


########################################################################################################################
num_channels = 1

def make_folder(data_folder_with_labels, input_folder, target_folder):
    shutil.copytree(data_folder_with_labels, target_folder, dirs_exist_ok=True)
    input_folder_list, target_folder_list = os.listdir(input_folder), os.listdir(target_folder)
    input_folder_list.sort(), target_folder_list.sort()
    for img_file, save_dir in zip(input_folder_list, target_folder_list):
        full_src = os.path.join(input_folder, img_file)
        full_dst = os.path.join(target_folder, save_dir, "input.tiff")
        print("copying from {} to {}".format(full_src, full_dst))
        shutil.copy(full_src, full_dst)

generator_exp = unet_generator(num_channels, 1, "instancenorm")
generator_sim = unet_generator(1, num_channels, "instancenorm")


def cycle(x):
    return generator_sim(generator_exp(x))

model = generator_exp if model_type == "exp" else cycle

try:
    checkpoint_path = "{}cycle_gan_results/checkpoints/checkpoint_{}".format(parent_dir, identifier)
    ckpt = tf.train.Checkpoint(generator_exp=generator_exp, generator_sim=generator_sim)
    ckpt_manager = tf.train.CheckpointManager(ckpt, checkpoint_path, max_to_keep=3)
    if ckpt_manager.latest_checkpoint:
        cpath = ckpt_manager.latest_checkpoint
        ckpt.restore(cpath)
except:
    print("loading checkpoint failed")
    exit()

GAN_image_folder(model, input_dir, GANNED_dir, 256, 256, False)
#make_folder(data_folder_with_labels, GANNED_dir, target_folder)











