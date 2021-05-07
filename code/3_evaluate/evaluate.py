import numpy as np
import sys
sys.path.insert(0, '../1_preprocessing')
sys.path.insert(0, '../2_trainer')
from image_parse import *
from accuracy import *
from tensorflow.keras.models import model_from_json
import matplotlib.pyplot as plt
from tifffile import imsave



def model_load(model_fn, model_weights_fn):
    with open(model_fn, 'r') as f:
            model = model_from_json(f.read())
    model.load_weights(model_weights_fn)
    return model


def predict(model, images):
    return model.predict_on_batch(images)


def stitch(size_x, size_y, sx, sy, images):
    '''
    Modified stitch function for better performance while processing larger image
    '''
    def stitch_sum(size_x, size_y, sx, sy, images):
        (num_cuts, lx, ly, num_class) = np.shape(images)
        Nx = len(np.arange(0, size_x - lx + 1, sx))
        Ny = len(np.arange(0, size_y - ly + 1, sy))

        L_X = ((num_cuts-1)  % Nx)*sx + lx
        L_Y = ((num_cuts-1) // Nx)*sy + ly
        final_img = np.zeros((L_X, L_Y, num_class))
        for idx, img in enumerate(images):
            nx = idx % Nx
            ny = idx //Nx
            final_img[nx*sx:nx*sx + lx,  ny*sy:ny*sy + ly, :] += img      
        return final_img
    stitch_result = stitch_sum(size_x, size_y, sx, sy, images)/stitch_sum(size_x, size_y, sx, sy, np.ones_like(images))
    return stitch_result


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
        plot=False, save_data=False, save_dir='./', fname="evaluated.tiff"):
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

    if plot:
        plt.figure()
        plt.imshow(a)
        plt.colorbar()
        plt.show()

    if save_data:
        print("saving data")
        imsave(save_dir + fname, a)

    return a

def calc_accuracy(evals_img, label_file_list, tol=0, nconvs=1, r=2, TN=0,
        plot=False, save_data=False, save_dir="./", prefix="", verbose=True):

    label_img = np.argmax(process_label(label_file_list, tol=tol), axis=2)

    conv_label_img = convolve(nconvs, label_img)
    conv_evals_img = convolve(nconvs, evals_img)

    conv_label_cen = get_center_list(conv_label_img, r)
    conv_evals_cen = get_center_list(conv_evals_img, r)

    match_list, label_list, evals_list = detect_diff(conv_label_cen, conv_evals_cen)


    fig = plt.figure(figsize=(10,10))
    m_xy      = [[],[]] if len(match_list)     == 0 else list(zip(*match_list))
    l_xy      = [[],[]] if len(label_list)     == 0 else list(zip(*label_list))
    e_xy      = [[],[]] if len(evals_list)     == 0 else list(zip(*evals_list))

    if plot:
        plt.scatter(list(l_xy[1]), list(l_xy[0]), c='r', label='FN')
        plt.scatter(list(e_xy[1]), list(e_xy[0]), c='k', label='FP')
        plt.scatter(list(m_xy[1]), list(m_xy[0]), c='b', label='TP')
        plt.legend(loc='best')
        plt.xlim(0, len(evals_img[0]))
        plt.ylim(0, len(evals_img))
        plt.gca().invert_yaxis()
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


def get_diagnostic_data(defect_dir_list, verbose=False):
    data_list = []
    for defect_dir in defect_dir_list:
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
        data_list.append(np.array(data))
    return data_list


def plot_diagnostics(data_list, label_list, diag="loss", log=True, invert=False, N=1, save=False, prefix=''):

    i1, i2 = 1, 3
    if   diag == "accuracy":
        i1, i2 = 2, 4
    elif diag == "recall":
        i1, i2 = 9, -1
    elif diag == "precision":
        i1, i2 = 10, -1

    def f(x):
        g = 1 - x if invert else x
        g = np.log10(g) if log else g
        return np.convolve(g, np.ones((N,))/N, mode='valid')

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    if i2 >=0:
        t, = ax.plot(f(data_list[0][:,i1]), 'k',   lw=2)
        v, = ax.plot(f(data_list[0][:,i2]), 'k--', lw=2)


    for i in range(len(data_list)):
        ax.plot(f(data_list[i][:,i1]), 'C{}'.format(i),  label=label_list[i], lw=2)
    if i2 >=0:
        for i in range(len(data_list)):
            ax.plot(f(data_list[i][:,i2]), 'C{}--'.format(i), lw=2)

    ax.tick_params(labelsize=15)
    plt.ylabel(diag, size=16)
    plt.xlabel("Number of epochs", size=16)
    leg1 = ax.legend(loc='best')
    if i2 >=0:
        leg2 = ax.legend([t,v],['Training','Validation'], loc='lower left')
    ax.add_artist(leg1)
    plt.tight_layout()
    if save:
        fig.savefig("{}loss.png".format(prefix), dpi=500)
    plt.show()


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



