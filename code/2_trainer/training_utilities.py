import os
import numpy as np
from models import *
from input_data import *
from accuracy import *
import matplotlib.pyplot as plt

def construct_model(N, k_fac, nb_classes, sess_dir, model_fn, model_weights_fn, dropout=0.1):
    input_img = Input(shape=(N,N,1))
    model = model_resunet(input_img, N, k_fac, nb_classes, dropout=dropout)
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics = ['accuracy'])
    model_json = model.to_json()

    with open(model_fn, 'w') as f:
            f.write(model_json)

    if os.path.isfile(model_weights_fn):
        print("loading weights")
        model.load_weights(model_weights_fn)
    return model


def setup_diagnostics(diagnostics_fn):
    step = 0
    if os.path.isfile(diagnostics_fn):
        with open(diagnostics_fn, 'r') as f:
            file_lines = f.readlines()
            if len(file_lines) > 1:
                step = int(file_lines[-1].split()[0])
                print("continuing session with step" + str(step))
    else:
        print("creating new session")
        with open(diagnostics_fn, "w") as f:
             f.write('{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\n'.format(\
                'step','loss','acc','val_loss','val_acc', 'TP', 'FP', 'FN', 'TN', 'recall', 'precision', 'F1', 'bal_acc'))
    return step

def calc_accuracy(model, x_test, y_true, N, nb_classes, plots):
    predictions = np.reshape(np.array(model.predict_on_batch(x_test)), [-1, N, N, nb_classes])
    y_evals = np.argmax(predictions, axis=3)
    y_trues = np.argmax(np.reshape(np.array(y_true), [-1, N, N, nb_classes]), axis=3)

    r = np.sum(y_evals)/(len(y_evals)*len(y_evals[0])*len(y_evals[0][0]))
    print(r)
    if r > 0.05:
        return -1, -1, -1, -1, -1, -1, -1, -1

    TP, FP, FN = 0, 0, 0
    for evals_img, label_img in zip(y_evals, y_trues):
        conv_label_img = convolve(1, label_img)
        conv_evals_img = convolve(1, evals_img)

        conv_label_cen = get_center_list(conv_label_img, 7.5)
        conv_evals_cen = get_center_list(conv_evals_img, 7.5)

        match_list, label_list, evals_list = detect_diff(conv_label_cen, conv_evals_cen)

        TP += len(match_list)
        FP += len(evals_list)
        FN += len(label_list)
    TN = 304*len(y_evals) - TP - FP - FN

    if plots:
        print("plotting")
        label_cen = [[],[]] if len(conv_label_cen) == 0 else list(zip(*conv_label_cen))
        evals_cen = [[],[]] if len(conv_evals_cen) == 0 else list(zip(*conv_evals_cen))
        m_xy      = [[],[]] if len(match_list)     == 0 else list(zip(*match_list))
        l_xy      = [[],[]] if len(label_list)     == 0 else list(zip(*label_list))
        e_xy      = [[],[]] if len(evals_list)     == 0 else list(zip(*evals_list))

        plt.subplot(131)
        plt.scatter(label_cen[0], label_cen[1], label='label')
        plt.legend(loc='best')
        plt.subplot(132)
        plt.scatter(evals_cen[0], evals_cen[1], label='evals')
        plt.legend(loc='best')
        plt.subplot(133)
        plt.scatter(list(m_xy[0]), list(m_xy[1]), label='match')
        plt.scatter(list(l_xy[0]), list(l_xy[1]), label='label')
        plt.scatter(list(e_xy[0]), list(e_xy[1]), label='evals')
        plt.legend(loc='best')

        plt.figure()
        plt.imshow(evals_img)
        plt.title("evaluation")

        plt.figure()
        plt.imshow(label_img)
        plt.title("label")

        plt.figure()
        plt.imshow(evals_img - label_img)
        plt.title("diff")
        plt.show()

    TNR = -1 if (TN + FP) == 0 else TN/(TN + FP)
    TPR = -1 if (TP + FN) == 0 else TP/(TP + FN)

    recall    = -1 if (TP + FN)            == 0 else TP/(TP + FN)
    precision = -1 if (TP + FP)            == 0 else TP/(TP + FP)
    F1        = -1 if (recall + precision) <= 0 else 2*recall*precision/(recall + precision)
    bal_acc   = -1 if (TNR + TPR)          <= 0 else 0.5*(TNR + TPR)

    return TP, FP, FN, TN, recall, precision, F1, bal_acc

def train_step(model, train_stem, test_stem, model_weights_fn, plots, epochs=1, batch_size=32):
    '''
    trains on the model for one step
    '''
    N = train_stem.N
    nb_classes = train_stem.nb_classes
    [tr_bs, ts_bs] = [train_stem.data._num_examples, test_stem.data._num_examples]
    def get_xy(batch_lbl):
        batch = train_stem.data.next_batch(tr_bs) if batch_lbl == 'train' else \
                test_stem.data.next_batch(ts_bs)
        return (np.reshape(batch[0], [-1,N,N,1]),  np.reshape(batch[1], [-1,N*N,nb_classes]))

    (x_train, y_train) = get_xy('train')
    (x_test , y_test ) = get_xy('test')

    history = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,\
            validation_data=(x_test, y_test), verbose=2)

    print("\tcalculating accuracy")
    TP, FP, FN, TN, recall, precision, F1, bal_acc = calc_accuracy(model, x_test, y_test, N, \
            nb_classes, plots)
    print("\tdone")

    model.save_weights(model_weights_fn)

    return history, TP, FP, FN, TN, recall, precision, F1, bal_acc


def train(step, data_dir, N, nb_classes, model, diagnostics_fn, model_weights_fn, num_steps=-1, plots=False):
    '''
    trains continuously until force stopped
    '''
    train_list = [f for f in os.listdir(data_dir) if "train" in f]
    test_f     = data_dir + "test_00000.p"
    test_stem  = read_data_set(test_f, N, nb_classes)

    num_files = len(train_list)
    steps_left = num_steps
    while steps_left > 0 or num_steps==-1:
        #pick a random augmented dataset from parsed_dir
        i = np.random.randint(num_files)
        train_f = train_list[i]
        print("training step: " + str(step) + "\ttraining file: " + train_f)

        # grab data for training
        print("\tgrabbing data")
        train_stem = grab_data(data_dir, train_f, N, nb_classes)
        print("\tdone")

        # train
        history, TP, FP, FN, TN, recall, precision, F1, bal_acc = train_step(model, train_stem,
                test_stem, model_weights_fn, plots)
        print("TP = {}, FP = {}, FN = {}, TN = {}".format(TP, FP, FN, TN))
        print("recall = {}, precision = {}, F1 = {}, bal_acc = {}".format(recall, precision, F1,\
                bal_acc))
        # record training step results
        loss     = history.history['loss'][0]
        val_loss = history.history['val_loss'][0]
        acc      = history.history['accuracy'][0]
        val_acc  = history.history['val_accuracy'][0]
        with open(diagnostics_fn, "a") as f:
             f.write('{:15d}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\n'.format(\
                step ,loss,acc,val_loss,val_acc,TP, FP, FN, TN, recall, precision, F1, bal_acc))
        step += 1
        steps_left -= 1
