# stem-learning
Package to use an FCN to learn defects from STEM imaging
## Generating Training Set
This part will show you why we need generate simualted STEM images as training set of deep learning model and how to use these code generate and process these images.
### Purpose
In deep learning, whether people can get a good model always depends on the amounts and quality of their training set. While in nano-scale, especially for electron microscope images, there are no specific images sets for training DL models. If people  want to use deep learning to recognize features, like point defects in electron microscope images, how to get enough high quality, labeled training image is the key question. In previous work, resaercher in EM field use several methods to get training set, including manually label and using fourier filter to label features. But these methods are time consuming or unaccurate. Here we come up a method to generate large amount simulated STEM images with high quality labels to train DL models.
