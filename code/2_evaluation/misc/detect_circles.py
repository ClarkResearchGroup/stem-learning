import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import imsave
import cv2
import pickle
from math import ceil, sin, cos
import sys
sys.path.insert(0, '../preprocessing')
from image_parse import *


def extract_defects(input_img, label_img, threshold=0.85, num_convs=5, radius=5, l=50):
    ''' takes in an input image and label image in the form of 2D numpy arrays and returns a
        list of images with defects at the centers. The size of each output image is l.
        The threshold value should be a positive number less than 1. the lower the threshold, the
        more it accepts questionable defects (i.e. labels that are not as intensely labeled)
        This depends on num_convs, which determines how many convolutions are done on the labeled image.
        radius is the radius of the circle that indicates a defect.
    '''

    (size_x, size_y) = label_img.shape

    def convolve(num_convs):
        print "convolving"
        kernel = np.ones((3,3))
        old_label_img = (label_img - np.min(label_img))/(np.max(label_img) - np.min(label_img))
        for k in range(num_convs):
            print k
            conv_label_img = np.array([[np.sum(np.multiply(kernel, old_label_img[i-1:i+2, j-1:j+2])) for j in range(1, size_y-1)] for i in range(1, size_x-1)])
            conv_label_img = np.array([ np.append([0], list(conv_label_img[i])) for i in range(size_x-2)])
            conv_label_img = np.array([ np.append(list(conv_label_img[i]), [0]) for i in range(size_x-2)])
            conv_label_img = np.append([np.zeros(size_y)], conv_label_img, axis=0)
            conv_label_img = np.append(conv_label_img, [np.zeros(size_y)], axis=0)

            conv_label_img = (conv_label_img - np.min(conv_label_img))/(np.max(conv_label_img) - np.min(conv_label_img))
            conv_label_img = np.array([[conv_label_img[i,j] if conv_label_img[i,j] > .6 else 0 for j in range(size_y)] for i in range(size_x)])

            old_label_img = np.copy(conv_label_img)
        return conv_label_img

    def threshold_points(conv_label_img, threshold, radius):
        print "thresholding"
        points_list = [ (i, j, conv_label_img[i,j]) for i in range(size_x) for j in range(size_y) ]
        center_list = []
        for (i, j, c) in points_list:
            if c > threshold:
                center_list.append((i, j, c))

        print "averaging"
        new_center_list = []
        while len(center_list) > 0:
            print len(center_list)
            avg_list = []
            (i, j, c) = center_list.pop(0)
            avg_list.append((i, j, c))

            k = len(center_list)
            while k > 0:
                (ik, jk, ck) = center_list.pop(0)
                rsq = (i - ik)*(i - ik) + (j - jk)*(j - jk)
                if rsq < radius*radius:
                    avg_list.append((ik, jk, ck))
                else:
                    center_list.append((ik, jk, ck))
                k -= 1

            avg_list = zip(*avg_list)
            i_new = int(round(np.sum(np.multiply(avg_list[0], avg_list[2])/np.sum(avg_list[2]))))
            j_new = int(round(np.sum(np.multiply(avg_list[1], avg_list[2])/np.sum(avg_list[2]))))
            c_new = np.mean(avg_list[2])
            new_center_list.append((i_new, j_new, c_new))

        print len(new_center_list), "centers found"
        return new_center_list


    def cut_defects(input_img, center_list, l):
        print "cutting defects"
        (x, y, c) = zip(*center_list)
        image_list = []
        num_defects_list = []
        for q, (i, j, k) in enumerate(center_list):
            x_start = i - l/2
            x_end   = i + l/2 + l%2
            y_start = j - l/2
            y_end   = j + l/2 + l%2
            if 0 <= x_start < x_end < size_x and 0 <= y_start < y_end < size_y:
                num_centers = sum(map(lambda r: x_start < r[0] < x_end and y_start < r[1] < y_end, center_list))
                im = input_img[x_start:x_end, y_start:y_end]
                image_list.append(im)
                num_defects_list.append(num_centers)
        return image_list,num_defects_list




    conv_label_img = convolve(num_convs)
    center_list = threshold_points(conv_label_img, threshold, radius)
    image_list, num_defects_list = cut_defects(input_img, center_list, l)
    new_label_img = (conv_label_img > threshold).astype(int)
    return image_list, center_list, conv_label_img, num_defects_list, new_label_img



def rotateImage(image, image_center, angle):
  rot_mat = cv2.getRotationMatrix2D(image_center, angle, 1.0)
  return cv2.warpAffine(image, rot_mat, image.shape[1::-1], flags=cv2.INTER_CUBIC)


def translateImage(image, t):
    translation_matrix = np.float32([ [1,0,t[0]], [0,1,t[1]] ])
    return cv2.warpAffine(image, translation_matrix, image.shape[1::-1])


def cropImage(image, center, size):
    x_start = center[0] - size[0]/2
    x_end   = center[0] + size[0]/2 + size[0]%2
    y_start = center[1] - size[1]/2
    y_end   = center[1] + size[1]/2 + size[0]%2
    return image[x_start:x_end, y_start:y_end]


def align_images(img_main, img_list, s=5, theta=5, size=(30,30)):
    ''' given a list of images and an image to compare with, this function
        translates and rotates the images to maximize the L2 norm between the
        images.
        s gives how many translations to provide.
        theta gives how many rotations to try
        size is the final size of the image to output
    '''
    (sx, sy) = img_main.shape
    (x_c, y_c) = (sx/2, sy/2)
    cropped_img_main = cropImage(img_main, (x_c,y_c), size)
    aligned_list = [cropped_img_main]

    for i, img in enumerate(img_list):
        print i
        best_img = None
        L2 = np.inf
        for i in range(-s/2, s/2 + 1):
            for j in range(-s/2, s/2 + 1):
                for th in range(-theta, theta + 1):
                    new_img = rotateImage(img, (x_c, y_c), th)
                    new_img = translateImage(new_img, (i, j))
                    new_img = cropImage(new_img, (x_c + i, y_c + j), size)
                    new_L2 = np.linalg.norm(new_img - cropped_img_main)
                    if new_L2 < L2:
                        best_img = new_img
                        L2 = new_L2
        aligned_list.append(best_img)
    return aligned_list


if __name__ == "__main__":

        final_image_size = (64, 64)
        theta_deg = 5


        l_size = int(ceil(final_image_size[0]*(sin(theta_deg*np.pi/180) + cos(theta_deg*np.pi/180))))


        # Locations of the images
	data_dir = "../../data/"
	dataset_dirs = ["01_8Mx_1K_0.019876nm_RGB_PPT/", "02_10Mx_1K_0.015806nm_RGB_PPT/", "03_8Mx_1K_0.019876nm_RGB_PPT/"]
	print "processing data"
	input_file = data_dir + "raw_data/" + dataset_dirs[0] + "input.png"
	label_file = data_dir + "raw_data/" + dataset_dirs[0] + "label2.png"


        # turns images to np arrays
	(input_img, label_img) = process_data(input_file, label_file, tol=0.5, diff=False)

        # IGNORE #######################################
	#pickle.dump(label_img, open( "label.p", "wb" ) )
	#import stitcher
	#label_img = stitcher.a
	#pickle.dump(label_img, open( "stitched_label.p", "wb" ) )
	#label_img = pickle.load( open( "stitched_label.p", "rb" ) )
        # IGNORE #######################################

        # NOTE the number of centers in an image returned here are within the lxl image and not the final cropped image from aligning
	image_list, center_list, conv_label_img, num_defects_list, new_label_img = extract_defects(input_img, label_img, threshold=0.85, num_convs=5, radius=5, l=l_size)

        # saves/loads image_list and num_center_list
        # you can save/load other variables here
	pickle.dump((image_list, num_center_list), open( "defect_data.p", "wb" ) )
	(image_list,num_center_list) = pickle.load( open( "defect_data.p", "rb" ) )

	aligned_images = align_images(image_list[0], image_list[1:], s=5, theta=theta_deg, size=final_image_size)


	for i, (img,num_center) in enumerate(zip(aligned_images, num_center_list)):
		imsave("aligned_defects/" + str(num_center) + "_defect_" + "_" + str(i).zfill(3) + ".png", img)

	(x, y, z) = zip(*center_list)
	plt.imshow(label_img)
	plt.scatter(y, x)
	plt.figure()
	plt.imshow(label_img)
	plt.figure()
	plt.imshow(conv_label_img)
	plt.show()






