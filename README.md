# stem-learning
Package to use an FCN to learn defects from STEM imaging
## Generating Training Set
This part will show you why we need generate simualted STEM images as training set of deep learning model and how to use these code generate and process these images.
### Purpose
In deep learning, whether people can get a good model always depends on the amounts and quality of their training set. While in nano-scale, especially for electron microscope images, there are no specific images sets for training DL models. If people  want to use deep learning to recognize features, like point defects in electron microscope images, how to get enough high quality, labeled training image is the key question. In previous work, researchers in EM field use several methods to get training set, including manually labeling and using fourier filter to label features. But these methods are time consuming or unaccurate. Here we come up a method to generate large amount simulated STEM images with high quality labels to train DL models.
Our method is based on computem written by Earl J. Kirkland, which is used for Transmission Electron Microscope Image Simulation. With computem people can get single high quality EM images, our purpose is to batchly generate EM images and its defect maps for training.
### How to use the code.
The code used for generating training sets hase two parts: the first part is used for generate .xyz and .param files, which can be input into computem to get simulate images and defect labels, the second part is used to add noises on simulated images to make them more "realistic".
