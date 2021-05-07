import numpy as np
import os
from params import *
from post_process import *
import shutil

def create_simulated_image():
    ### generate param and xyz files ####
    generate_files(sample_param_dic, EM_param_dic, 1)
    #####################################

    ### use incostem to create images ###
    bat_file = [f for f in os.listdir() if ".bat" in f][0]
    with open(bat_file, 'r') as bat:
        for cmd in bat:
            print(cmd)
            os.system(cmd)

    [os.remove(f) for f in os.listdir() if ".xyz" in f or ".param" in f or ".bat" in f]
    #####################################

    ### do postprocessing

    def experimentalize(Mypostprocess):
        #Mypostprocess.add_horizental_sheer((0.05,0.025))
        #Mypostprocess.add_vertical_constrain((0.05,0.025))
        Mypostprocess.crop(1024,1024)


        Mypostprocess.image_stacks -= 10
        Mypostprocess.image_stacks *= .0175/117.5
        Mypostprocess.image_stacks += .086
        Mypostprocess.add_gaussian_noise(0,.0025)

        #bkg_stack = tifffile.imread('./bkg/Bkg_stack.tif')
        Mypostprocess.add_bkg(bkg_stack)

    Mypostprocess = post_process(image_path='./',file_num=1,defect_list = defect_list)
    Mypostprocess.read_image_and_label()
    #experimentalize(Mypostprocess)
    Mypostprocess.save_as_image('./sim_image/')

    [os.remove(f) for f in os.listdir() if ".tif" in f]







