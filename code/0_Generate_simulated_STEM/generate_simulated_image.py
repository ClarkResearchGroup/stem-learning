import os
from params import *
from post_process import *

def create_simulated_image(incostem_dir, save_path="./sim_image/", experimentalize=False, bkg_file=None, num_files=1, save_as_stack=False):
    cur_dir = os.getcwd()
    os.chdir(incostem_dir)
    ### generate param and xyz files ####
    generate_files(sample_param_dic, EM_param_dic, num_files)
    #####################################

    ### use incostem to create images ###
    bat_file = [f for f in os.listdir() if ".bat" in f][0]
    with open(bat_file, 'r') as bat:
        for cmd in bat:
            print(cmd)
            if os.system(cmd):
                print("ERROR EXECUTING COMMAND: ", cmd)


    [os.remove(f) for f in os.listdir() if ".xyz" in f or ".param" in f or ".bat" in f]
    #####################################

    ### do postprocessing
    Mypostprocess = post_process(image_path="./",file_num=num_files,defect_list = defect_list)
    Mypostprocess.read_image_and_label()
    if experimentalize:
        Mypostprocess.experimentalize(bkg_file)
    if save_as_stack:
        Mypostprocess.save_stack(save_path)
    else:
        Mypostprocess.save_as_image(save_path)

    if bkg_file:
        os.rename(bkg_file, "tmp")
    [os.remove(f) for f in os.listdir() if ".tif" in f]
    if bkg_file:
        os.rename("tmp", bkg_file)
    
    print("Done generating simulated images!")
    os.chdir(cur_dir)