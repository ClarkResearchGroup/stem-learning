import numpy as np
from image_parse import *
from evaluator import predict, model_load
from sys import argv
import theano
from keras import backend as K


#arg = "08-23-2018+1"
arg = argv[1]
exp_idx = int(argv[2]) - 1

parent_dir = "../../"

exp_dirs = ["01_8Mx_1K_0.019876nm_RGB_PPT/", "02_10Mx_1K_0.015806nm_RGB_PPT/", "03_8Mx_1K_0.019876nm_RGB_PPT/", "02_10Mx_1K_0.015806nm_RESIZED/"]
print "loading model"
model_fn = parent_dir + "model/model" + arg + ".h5"
model_weights_fn = parent_dir + "model/model" + arg + "_weigths.h5"
model = model_load(model_fn, model_weights_fn)

print "processing data"
data_dir = parent_dir + "data/"
dataset_dir = "raw_data/" + exp_dirs[exp_idx]
#dataset_dir = "simulated_data/"
input_file = data_dir + dataset_dir + "input.png"
label_file = data_dir + dataset_dir + "label.png"

(input_img, label_img) = process_data(input_file, label_file)#, tol=0.5, diff=False)
(size_x, size_y) = label_img.shape

print "cutting data"
(sx, sy) = (32, 32)
(lx, ly) = (128, 128)
input_cuts = np.array(cut_data(input_img, (lx, ly), (sx, sy), rotflip=False))
input_cuts = np.reshape(input_cuts, [-1, lx, ly, 1])


def get_dims(val):
        fac_list = [(i, val / i) for i in range(1, int(val**0.5)+1) if val % i == 0]
        num_facs = len(fac_list)
        return fac_list[-1]


print "plotting"
import matplotlib.pyplot as plt
img_to_visualize = 0
plt.figure()
plt.title("Image used: #%d (digit=%d)" % (img_to_visualize, 1))
plt.imshow(np.reshape(input_cuts[img_to_visualize], (128,128)))



for layer in range(1, len(model.layers)-2):

    convout1_f = K.function([model.layers[0].input], [model.layers[layer].output])
    convolutions = convout1_f([input_cuts[:5]])
    convolutions = convolutions[0]
    convolutions = np.transpose(convolutions, (0,3,1,2))

    (n1, n2) = get_dims(convolutions.shape[1])

    fig, axes = plt.subplots(n1, n2, figsize=(15,15), squeeze=False)
    fig.suptitle("layer " + str(layer))

    print "The second dimension tells us how many convolutions do we have: %s (%d convolutions)" % (str(convolutions.shape),convolutions.shape[1])

    for idx, convolution in enumerate(convolutions[img_to_visualize]):
        i = idx/n2
        j = idx - i*n2
        axes[i,j].imshow(convolution)
    plt.savefig("image%s_layer_%s" % (str(img_to_visualize).zfill(3), str(layer).zfill(2)))


'''
print "predicting data"
predictions = predict(model, input_cuts)
predictions = [np.reshape(predictions[i,:,1], [lx,ly]) for i in range(len(predictions))]

print "stitching data"
a = stitch(size_x, size_y, sx, sy, predictions)
a = ( a > .5).astype(int)
b = np.abs(a - label_img)

print "plotting data"
if avg:
    avg_path = "avg/"
else:
    avg_path = "non_avg/"

path = "/home/aakhan3/" + arg + "/" + exp_dirs[exp_idx] + avg_path
import matplotlib.pyplot as plt
plt.imshow(a)
plt.savefig(path + "prediction.png")
plt.figure()
plt.imshow(label_img)
plt.savefig(path + "label.png")
plt.figure()
plt.imshow(b, vmin = 0, vmax = 1)
plt.savefig(path + "diff.png")
plt.figure()
'''
