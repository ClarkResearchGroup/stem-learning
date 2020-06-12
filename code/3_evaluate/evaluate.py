import numpy as np
import sys
from matplotlib.pyplot import imsave
sys.path.insert(0, '../1_preprocessing')
sys.path.insert(0, '../2_trainer')
from image_parse import *
from accuracy import *
from keras.models import model_from_json
import matplotlib.pyplot as plt



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
        ny = idx // Nx
        [ final_img[nx*sx + x][ny*sy + y].append(list(img[x,y])) for x in range(lx) for y in range(ly) ]
    for i in range(len(final_img)):
        for j in range(len(final_img[0])):
            if len(final_img[i][j]) == 0:
                final_img[i][j].append([1, 0])

    ret = [[np.mean(np.array(final_img[i][j]), axis=0) for j in range(size_y)] \
            for i in range(size_x)]

    return np.nan_to_num(np.array(ret))


def get_avg_pred(model, cut):
    (lx, ly, nb_classes) = cut.shape
    rot_cuts  = np.array([ np.rot90(cut, k) for k in range(4) ])
    flip_cuts = np.array([ np.fliplr(r_cut) for r_cut in rot_cuts ])

    pred_rots  = predict(model, rot_cuts)
    pred_flips = predict(model, flip_cuts)

    pred_rots  = np.reshape(pred_rots , [4, lx, ly, pred_rots.shape[-1]])
    pred_flips = np.reshape(pred_flips, [4, lx, ly, pred_flips.shape[-1]])

    pred_1 = [ np.rot90(pred_rots[k], -k)             for k in range(4) ]
    pred_2 = [ np.rot90(np.fliplr(pred_flips[k]), -k) for k in range(4) ]

    preds = pred_1 + pred_2

    return sum(preds)/8.0

def get_s(size_x, lx):
    '''
    computes the stride so that the whole input image
    is processed
    '''
    #for sx in range(lx, 0, -1):
    #    if (size_x - lx) % sx == 0:
    #        return sx
    return lx - 16


def evaluate(model_fn, model_weights_fn, input_file, l_shape, stride=None, avg=False, \
        plot=False, save_data=False, save_dir='./', prefix=""):
    ''' evauates an input image given the model'''

    print("processing data")
    input_img = process_image(input_file, standardize=True)
    (size_x, size_y) = input_img.shape

    print("loading model")
    model = model_load(model_fn, model_weights_fn)
    (lx, ly) = l_shape
    (sx, sy) = stride if stride != None else (get_s(size_x, lx), get_s(size_y, ly))
    print("strides: ({}, {})".format(sx, sy))

    input_cuts = cut_data(input_img, lx, ly, (sx, sy))
    input_cuts = np.reshape(input_cuts, [-1, lx, ly, 1])
    num_cuts = len(input_cuts)

    print("predicting data")
    predictions = np.array([get_avg_pred(model,cut) for cut in input_cuts] if avg else predict(model,input_cuts))
    predictions = np.reshape(np.array(predictions), [num_cuts, lx, ly, -1])

    print("stitching data")
    a = stitch(size_x, size_y, sx, sy, predictions)
    a = np.argmax(a, axis=2)

    if save_data:
        print("saving data")
        imsave(save_dir + prefix + "prediction.png", a)

    if plot:
        plt.figure()
        plt.imshow(a)
        plt.colorbar()
        plt.show()
    return a

def calc_accuracy(evals_img, label_file_list, tol=0, nconvs=1, r=2, TN=0,
        plot=False, save_data=False, save_dir="./", prefix="", verbose=True):

    label_img = np.argmax(process_label(label_file_list, tol=tol), axis=2)

    conv_label_img = convolve(nconvs, label_img)
    conv_evals_img = convolve(nconvs, evals_img)

    conv_label_cen = get_center_list(conv_label_img, r)
    conv_evals_cen = get_center_list(conv_evals_img, r)

    match_list, label_list, evals_list = detect_diff(conv_label_cen, conv_evals_cen)


    fig = plt.figure()
    m_xy      = [[],[]] if len(match_list)     == 0 else list(zip(*match_list))
    l_xy      = [[],[]] if len(label_list)     == 0 else list(zip(*label_list))
    e_xy      = [[],[]] if len(evals_list)     == 0 else list(zip(*evals_list))

    plt.scatter(list(l_xy[0]), list(l_xy[1]), c='r', label='FN')
    plt.scatter(list(e_xy[0]), list(e_xy[1]), c='k', label='FP')
    plt.scatter(list(m_xy[0]), list(m_xy[1]), c='b', label='TP')
    plt.legend(loc='best')
    plt.xlim(0, len(evals_img))
    plt.ylim(0, len(evals_img[0]))
    plt.gca().invert_yaxis()

    if plot:
        plt.show()
    if save_data:
        fig.savefig("{}{}accuracy_plot.png".format(save_dir, prefix))

    TP = len(match_list)
    FP = len(evals_list)
    FN = len(label_list)

    TNR = -1 if (TN + FP) == 0 else TN/(TN + FP)
    TPR = -1 if (TP + FN) == 0 else TP/(TP + FN)

    recall    = -1 if (TP + FN)            == 0 else TP/(TP + FN)
    precision = -1 if (TP + FP)            == 0 else TP/(TP + FP)
    F1        = -1 if (recall + precision) <= 0 else 2*recall*precision/(recall + precision)
    bal_acc   = -1 if (TNR + TPR)          <= 0 else 0.5*(TNR + TPR)

    if verbose:
        print("TP: {}".format(TP))
        print("FP: {}".format(FP))
        print("FN: {}".format(FN))
        print("TN: {}".format(TN))
        print("")
        print("recall:    {}".format(recall))
        print("precision: {}".format(precision))
        print("F1:        {}".format(F1))
        print("bal_acc:   {}".format(bal_acc))

    return TP, FP, FN, TN, recall, precision, F1, bal_acc


def get_diagnostic_data(defect_dir, verbose=False):
    d = open(defect_dir + "diagnostics.dat", 'r')
    lines = [line.split() for line in d]
    vals = lines[0]
    data = []
    for line in lines:
        try:
            if len(line) == len(vals):
                x = np.array(line).astype(float)
                x[x<0] = 0
                data.append(x)
        except:
            continue
    if verbose:
        print("values: {}".format(vals))
        print("number of data points: {}".format(len(data)))
    return np.array(data)


def get_diff(evals_img, label_file_list, tol=.1, plot=False, save_data=False,
        save_dir="./", prefix=""):

    label_img = process_label(label_file_list, tol=tol)
    label_img = np.argmax(label_img, axis=2)
    diff = evals_img - label_img

    if plot:
        plt.figure()
        plt.imshow(diff)
        plt.colorbar()

    if save_data:
        print("saving data")
        imsave(save_dir + prefix + "diff.png", diff)

    return diff



