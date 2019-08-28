import numpy as np
from scipy.signal import convolve2d


def convolve(num_convs, label_img, thresh=0.6):
    kernel = np.ones((3,3))/9.
    conv_label_img = (label_img - np.min(label_img))/np.ptp(label_img)
    conv_label_img = np.array(conv_label_img >= thresh).astype(int)
    for k in range(num_convs):
        conv_label_img = convolve2d(conv_label_img, kernel, mode='same')
        conv_label_img = np.array(conv_label_img >= thresh).astype(int)
    return conv_label_img


def get_center_list(conv_label_img, radius):
    size_x, size_y = conv_label_img.shape
    center_list = [[i, j] for i in range(size_x) for j in range(size_y)\
                  if conv_label_img[i,j]==1]

    new_center_list = []
    while len(center_list) > 0:
        avg_list = []
        i, j = center_list.pop(0)
        avg_list.append([i, j])

        for k in range(len(center_list)):
            ik, jk = center_list.pop(0)
            rsq = (i - ik)*(i - ik) + (j - jk)*(j - jk)
            if rsq < radius*radius:
                avg_list.append([ik, jk])
            else:
                center_list.append([ik, jk])

        avg_list = zip(*avg_list)
        i_new = int(round(np.mean(avg_list[0])))
        j_new = int(round(np.mean(avg_list[1])))
        new_center_list.append((i_new, j_new))

    #print len(new_center_list), "centers found"
    return new_center_list


def detect_diff(label_center, evals_center, radius=7.5):
    match_list = []
    label_list = list(label_center)
    evals_list = list(evals_center)

    for j in range(len(label_list)):
        (lx, ly) = label_list.pop(0)
        match_found = False
        for k in range(len(evals_list)):
            (ex, ey) = evals_list[k]
            rsq = (lx - ex)*(lx - ex) + (ly - ey)*(ly - ey)
            if rsq <= radius*radius:
                evals_list.pop(k)
                match_coord = ((lx + ex)/2, (ly + ey)/2)
                match_list.append(match_coord)
                match_found = True
                break
        if not match_found:
            label_list.append((lx, ly))
    return match_list, label_list, evals_list

def calculate_accuracy(label_file_list, evals_file_list, num_convs):

    label_img = process_label(label_file_list)[:,:,1]
    evals_img = process_label(evals_file_list)[:,:,1]

    conv_label_img = convolve(num_convs, label_img)
    conv_evals_img = convolve(num_convs, evals_img)

    conv_label_cen = get_center_list(conv_label_img, 7.5)
    conv_evals_cen = get_center_list(conv_evals_img, 7.5)

    match_list, label_list, evals_list = detect_diff(conv_label_cen, conv_evals_cen)

    TP = len(match_list)
    FP = len(evals_list)
    FN = len(label_list)
    TN = 304*16 - TP - FP - FN

    return TP ,FP, FN, TN
