import numpy as np
from os.path import isdir, isfile, join
from os import mkdir, listdir
from shutil import rmtree
from imageio import imread, imwrite
from PIL import Image

def _list_images(imgdir, ftype):
    return [f for f in listdir(imgdir) if isfile(join(imgdir, f)) and ftype in f]

def _list_dirs(augdir):
    return [f for f in listdir(augdir) if isdir(join(augdir,f))]


def make_augments(imgdir, ftype):
    augdir = join(imgdir,"augments")
    if isdir(augdir):
        rmtree(augdir)
    mkdir(augdir)


    print("inverting")
    mkdir(join(augdir,"orig"))
    mkdir(join(augdir,"flip"))
    for img in _list_images(imgdir, ftype):
        img_arr = imread(join(imgdir,img))
        img_flp = np.fliplr(img_arr)
        imwrite(join(augdir, "orig",img), img_arr)
        imwrite(join(augdir, "flip",img), img_flp)

    print("rotating")
    for d in _list_dirs(augdir):
        mkdir(join(augdir, "rot0_" + d))
        mkdir(join(augdir, "rot1_" + d))
        mkdir(join(augdir, "rot2_" + d))
        mkdir(join(augdir, "rot3_" + d))
        for img in _list_images(join(augdir, d), ftype):
            img_arr = imread(join(augdir, d, img))
            img_rt1 = np.rot90(img_arr, 1)
            img_rt2 = np.rot90(img_arr, 2)
            img_rt3 = np.rot90(img_arr, 3)

            imwrite(join(augdir, "rot0_" + d,img), img_arr)
            imwrite(join(augdir, "rot1_" + d,img), img_rt1)
            imwrite(join(augdir, "rot2_" + d,img), img_rt2)
            imwrite(join(augdir, "rot3_" + d,img), img_rt3)
        rmtree(join(augdir, d))

    print("magnifying")
    for d in _list_dirs(augdir):
        mkdir(join(augdir, "mag0_" + d))
        mkdir(join(augdir, "mag1_" + d))
        mkdir(join(augdir, "mag2_" + d))
        for img in _list_images(join(augdir, d), ftype):
            img_arr = imread(join(augdir, d, img))
            imwrite(join(augdir, "mag0_" + d,img), img_arr)
            imwrite(join(augdir, "mag1_" + d,img), img_arr)
            imwrite(join(augdir, "mag2_" + d,img), img_arr)

        f  = join(augdir, d, "input" + ftype)
        img_arr = Image.open(f)
        l, w  = img_arr.size
        l2 = int(l/2)
        l4 = int(l/4)

        l2_img = (img_arr.resize((l2, l2), Image.LANCZOS)).resize((l,l), Image.LANCZOS)
        l4_img = (img_arr.resize((l4, l4), Image.LANCZOS)).resize((l,l), Image.LANCZOS)

        l2_img.save(join(augdir, "mag1_" + d, "input" + ftype))
        l4_img.save(join(augdir, "mag2_" + d, "input" + ftype))

        rmtree(join(augdir, d))






if __name__ == "__main__":
    imgdir = "/home/abid/Dropbox/Development/programs/stem-learning/data/WSeTe/simulated/0"
    ftype = ".tiff"
    make_augments(imgdir, ftype)

