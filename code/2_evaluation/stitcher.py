from keras.models import model_from_json
from image_parse import *
import numpy as np
import sys
from scipy.misc import imsave
sys.path.insert(0, '../preprocessing')


def model_load(model_fn, model_weights_fn):
    with open(model_fn, 'r') as f:
        model = model_from_json(f.read())
    model.load_weights(model_weights_fn)
    return model


def predict(model, images):
    return model.predict_on_batch(images)


def stitch(size_x, size_y, sx, sy, images):
    '''
    Takes in a set of images, and stitches them into one image
    '''
    num_images = len(images)
    lx = len(images[0])
    ly = len(images[0][0])
    Nx = len(np.arange(0, size_x - lx + 1, sx))

    final_img = [[[] for j in range(size_y)] for i in range(size_x)]
    for idx, img in enumerate(images):
        nx = idx % Nx
        ny = idx / Nx
        [final_img[nx*sx + x][ny*sy + y].append(img[x, y]) for x in range(lx) for y in range(ly)]

    ret = [[np.median(np.array(final_img[i][j]), axis=0) for j in range(size_y)]
           for i in range(size_x)]

    return np.nan_to_num(np.array(ret))


def get_avg_pred(model, cut):
    (lx, ly, nb_classes) = cut.shape
    rot_cuts = np.array([np.rot90(cut, k) for k in range(4)])
    flip_cuts = np.array([np.fliplr(r_cut) for r_cut in rot_cuts])

    pred_rots = predict(model, rot_cuts)
    pred_flips = predict(model, flip_cuts)

    pred_rots = np.reshape(pred_rots, [4, lx, ly, pred_rots.shape[-1]])
    pred_flips = np.reshape(pred_flips, [4, lx, ly, pred_flips.shape[-1]])

    pred_1 = [np.rot90(pred_rots[k], -k) for k in range(4)]
    pred_2 = [np.rot90(np.fliplr(pred_flips[k]), -k) for k in range(4)]

    preds = pred_1 + pred_2

    return sum(preds)/8.0


def make_prediction(model_fn, model_weights_fn, input_file, label_file_list, Tol, avg, save_dir,
                    thresh=-1, prefix="", plot=False, save_data=False):

    print "processing data"
    input_img = process_image(input_file)
    label_img = process_label(label_file_list, tol=Tol)
    (size_x, size_y, nb_classes) = label_img.shape

    print "loading model"
    model = model_load(model_fn, model_weights_fn)

    print "cutting data"
    (sx, sy) = (32, 32)
    (lx, ly) = (128, 128)
    input_cuts = cut_data(input_img, lx, ly, (sx, sy))
    input_cuts = np.reshape(input_cuts, [-1, lx, ly, 1])
    num_cuts = len(input_cuts)

    print "predicting data"
    predictions = [get_avg_pred(model, cut) for cut in input_cuts] if avg else predict(model, input_cuts)
    predictions = np.reshape(np.array(predictions), [num_cuts, lx, ly, nb_classes])

    print "stitching data"
    a = stitch(size_x, size_y, sx, sy, predictions)
    a = np.argmax(a, axis=2)
    label_img = np.argmax(label_img, axis=2)
    b = a - label_img

    if save_data:
        print "saving data"
        imsave(save_dir + prefix + "prediction.png", a)
        imsave(save_dir + prefix + "label.png", label_img)
        imsave(save_dir + prefix + "diff.png", b)

    if plot:
        if plot:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.imshow(a)
            plt.colorbar()

            plt.figure()
            plt.imshow(label_img)
            plt.colorbar()

            plt.figure()
            plt.imshow(b)
            plt.colorbar()

            plt.show()
