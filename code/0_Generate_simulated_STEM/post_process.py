# Import necessary package
import numpy as np
import os
from skimage import transform
import tifffile
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (16, 9)

defect_list = ['1Doped','2Doped','1vacancy','2vacancy', 'metal_vacancy', 'metal_Doped']

class post_process():

    #Two gloabal helper functions#
    global _get_normal_distribution
    global _generate_random_bkg


    def __init__(self, image_path, file_num, defect_list):

        #Input parameters of this class:
        #string: image_path, the path of you images and labels
        #int: file_num
        #list: defect_list, list of str, includes the defect labels you need

        self.image_path = image_path
        self.file_num = file_num
        self.defect_list = defect_list


    def _get_normal_distribution(param):
        '''
        Get positive/negative normal distribution from a tuple
        Imput: tuple: (mean, std)
        return: float: normal_distribution_number
        '''
        num = np.random.normal(param[0],param[1])
        pos = np.random.rand()

        if pos>0.5:
            return num
        else:
            return -num


    def read_image_and_label(self):

        '''
        Put image and its defect into one image stack for following process
        self.image_stack: np.array shape = (num_image, x,y,num_defects+1)
        '''

        all_files = os.listdir(self.image_path)
        image_list = [image_file for image_file in all_files if image_file[:5] == 'Image']
        img = plt.imread(self.image_path+image_list[0])
        #get image size
        x,y = np.shape(img)
        image_stacks = np.zeros((self.file_num,x,y,len(self.defect_list)+1))
        for i in range(self.file_num):
            image_name = image_list[i]
            temp_image = plt.imread(self.image_path+image_name)
            image_stacks[i,:,:,0] = temp_image
            for j in range(len(self.defect_list)):
                temp_defect_name = self.defect_list[j]+'_'+image_name[5:]
                try:
                    temp_defect_image = plt.imread(self.image_path+temp_defect_name)
                except:
                    temp_defect_image = np.zeros((x,y))
                image_stacks[i,:,:,j+1] = temp_defect_image
        self.image_stacks = image_stacks.copy()

    def add_horizental_shear(self, shear_rate):

        """
        Add horizental shear to simulate horizental sample draft during taking the image in STEM
        shear_rate = (shear_mean,shear_std) in guassian distribution
        single image has a linear shear, while in different images, shear rates are different
        """
        num_image, x, y, num_defects = np.shape(self.image_stacks)
        for i in range(num_image):
            temp_stack = self.image_stacks[i,:,:,:]
            shear = _get_normal_distribution(shear_rate)
            afine_tf = transform.AffineTransform(shear = shear)
            modified = transform.warp(temp_stack, inverse_map=afine_tf)
            self.image_stacks[i,:,:,:] = modified




    def add_vertical_contraction(self, contraction_rate):
        '''
        Add vertical contraction to simulate vertical sample draft during taking the image in STEM
        contraction_rate = (contraction_mean,contraction_std) in Gaussian distribution
        single image has a constant contraction rate, while in different images, contraction rate are different
        '''

        num_image, x, y, num_defects = np.shape(self.image_stacks)
        for i in range(num_image):
            temp_stack = self.image_stacks[i,:,:,:]
            contraction = _get_normal_distribution(contraction_rate)+1
            afine_tf = transform.AffineTransform(scale = (1,contraction))
            modified = transform.warp(temp_stack, inverse_map=afine_tf)
            self.image_stacks[i,:,:,:] = modified



    def rotate(self, degree):

        '''
        rotate the images into the same degree of real images, if there are multiple orientations in real images,
        you can genreate different degrees in training set
        '''

        num_images,x,y,num_defects = np.shape(self.image_stacks)
        for i in range(num_images):
            for j in range(num_defects):
                 self.image_stacks[i,:,:,j] = transform.rotate(self.image_stacks[i,:,:,j],degree)
        pass

    def crop(self, target_x, target_y, rnd_shift):

        '''
        The raw images is a bit larger because the shear and rotation will cause some blank in the image
        Crop the center image to shape(target_x, target_y) and avoid the blank
        shift the cropping position a little bit
        '''

        num_image, x, y, num_defects = np.shape(self.image_stacks)
        angle = 2*np.pi*np.random.rand()
        shift_vec = rnd_shift*np.array([np.cos(angle),np.sin(angle)])
        cropped_x = abs(x-target_x)//2 + int(shift_vec[0])
        cropped_y = abs(y-target_y)//2 + int(shift_vec[1])
        self.image_stacks = self.image_stacks[:,cropped_x:cropped_x+target_x,cropped_y:cropped_y+target_y,:]
        pass

    def change_brightness_and_contrast(self, target_mean, target_std):

        '''
        Change the B/C of simulate images to the real image
        You can get the target_mean and target_std from a real image

        '''

        num_images,x,y,num_defects = np.shape(self.image_stacks)
        for i in range(num_images):
            temp_image = self.image_stacks[i,:,:,0]
            temp_image += target_mean - np.mean(temp_image)
            temp_image *= target_std/np.std(temp_image)

            self.image_stacks[i,:,:,0] = temp_image


    def add_gaussian_noise(self, mean, std):

        '''
        add gaussian noise in the image
        '''

        for i in range(self.file_num):
            self.image_stacks[i,:,:,0] += np.random.normal(mean,std,self.image_stacks[i,:,:,0].shape)

    def _generate_random_bkg(bkg_stack):

        num, x, y = np.shape(bkg_stack)
        i = np.random.randint(num)
        bkg_image = bkg_stack[i,:,:]
        bkg_image = np.rot90(bkg_image, np.random.randint(4))
        return bkg_image



    def add_bkg(self, bkg_stack, bkg_strength):

        '''
        Add random, ununiform background on the images,
        you can extract bkg through high pass filter of STEM image
        '''

        _,x,y,_ = np.shape(self.image_stacks)
        for i in range(self.file_num):
            bkg_image = _generate_random_bkg(bkg_stack)
            self.image_stacks[i,:,:,0] += bkg_strength * bkg_image[:x,:y]
            self.image_stacks[i,:,:,1:] += bkg_strength * np.mean(bkg_image[:x,:y]) #Keep the image and lables with the same intensity range
        pass
    
    def add_jitter(self, jitter_std):
        
        '''
        Add Gaussian distributed jittering to add streak-like feature to the STEM images
        '''
        _,x,y,_ = np.shape(self.image_stacks)
        for i in range(self.file_num):
            jitter_arr = np.int8(np.random.normal(0, jitter_std, size=(x)))
            for line_num in range(x):
                self.image_stacks[i,line_num,:,:] = np.roll(self.image_stacks[i,line_num,:,:], jitter_arr[line_num], axis = 0)
        pass
    
    def save_as_image(self, save_path):
        '''
        Save images stacks as image
        '''
        if os.path.exists(save_path) == False:
            os.makedirs(save_path)
        num_images,x,y,num_defects = np.shape(self.image_stacks)

        for i in range(num_images):
            temp_path = save_path+'image{}/'.format(i)
            if os.path.exists(temp_path) == False:
                os.makedirs(temp_path)
            tifffile.imsave(temp_path+'input.tiff',self.image_stacks[i,:,:,0].astype(np.float32))
            for j in range(len(self.defect_list)):
                tifffile.imsave(temp_path+'label_'+defect_list[j]+'.tiff',
                                self.image_stacks[i,:,:,j+1].astype(np.float32))
    
    def save_stack(self, save_path):
        if os.path.exists(save_path) == False:
            os.makedirs(save_path)
        num_images = np.shape(self.image_stacks)[0]
        
        for i in range(num_images):
            tifffile.imsave(save_path+'sim_stack{}.tiff'.format(i),np.rollaxis(self.image_stacks[i,:,:,:], 2, 0).astype(np.float32))
        tifffile.imsave(save_path+'sim_images_stack.tiff', self.image_stacks[:,:,:,0].astype(np.float32))

    def save_as_npy_files(self, save_path):

        '''
        Save images files as .npy files, it is easy to load into deep learning code
        '''

        if os.path.exists(save_path) == False:
            os.makedirs(save_path)

        np.save(save_path+'images_files.npy',self.image_stacks[:,:,:,0])
        for i in range(len(self.defect_list)):
            np.save(save_path+defect_list[i]+'_files.npy',self.image_stacks[:,:,:,i+1])
    
    def experimentalize(self, bkg_file=None):
        self.add_horizental_shear((0.01,0.01))
        self.add_vertical_contraction((0.01,0.01))
        self.rotate(-90)
        self.add_jitter(1)
        self.crop(1024, 1024, rnd_shift=10)
        
        self.image_stacks -= self.image_stacks.min()
        self.image_stacks *= 0.0174/115 #Get these numbers from noise-free images
        self.image_stacks += 0.0034
        self.add_gaussian_noise(0,.0030)
        if bkg_file:
            bkg_stack = tifffile.imread(bkg_file)
            self.add_bkg(bkg_stack, 1)
        self.image_stacks *= 65535 #Get these numbers from noise-free images



