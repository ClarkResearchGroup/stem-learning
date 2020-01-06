import os
import numpy as np
from models import *
from input_data import *
from accuracy import *

def construct_model(N, k_fac, nb_classes, sess_dir, model_fn, model_weights_fn):
    input_img = Input(shape=(N,N,1))
    model = model_resunet(input_img, N, k_fac, nb_classes)
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
             f.write('{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\n'.format(\
                'step','loss','acc','val_loss','val_acc', 'recall', 'precision', 'F1', 'bal_acc'))
    return step

def calc_accuracy(model, x_test, y_true, N, nb_classes):
    predictions = np.reshape(np.array(model.predict_on_batch(x_test)), [-1, N, N, nb_classes])
    y_evals = np.argmax(predictions, axis=3)
    y_trues = np.argmax(np.reshape(np.array(y_true), [-1, N, N, nb_classes]), axis=3)
    TP, FP, FN, TN = 0, 0, 0, 0
    for evals_img, label_img in zip(y_evals, y_trues):
        conv_label_img = convolve(2, label_img)
        conv_evals_img = convolve(2, evals_img)

        conv_label_cen = get_center_list(conv_label_img, 7.5)
        conv_evals_cen = get_center_list(conv_evals_img, 7.5)

        match_list, label_list, evals_list = detect_diff(conv_label_cen, conv_evals_cen)

        TP += len(match_list)
        FP += len(evals_list)
        FN += len(label_list)
        TN += 304 - TP - FP - FN

    TNR = TN/(TN + FP)
    TPR = TP/(TP + FN)

    recall    = TP/(TP + FN)
    precision = TP/(TP + FP)
    F1        = 2*recall*precision/(recall + precision)
    bal_acc   = 0.5*(TNR + TPR)
    return recall, precision, F1, bal_acc

def train_step(model, stem, model_weights_fn, epochs=1, batch_size=32):
    '''
    trains on the model for one step
    '''
    N = stem.N
    nb_classes = stem.nb_classes
    [tr_bs, ts_bs] = [stem.train._num_examples, stem.test._num_examples]
    def get_xy(batch_lbl):
        batch = stem.train.next_batch(tr_bs) if batch_lbl == 'train' else stem.test.next_batch(ts_bs)
        return (np.reshape(batch[0], [-1,N,N,1]),  np.reshape(batch[1], [-1,N*N,nb_classes]))

    (x_train, y_train) = get_xy('train')
    (x_test , y_test ) = get_xy('test')


    history = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,\
            validation_data=(x_test, y_test), verbose=1)

    print("\tcalculating accuracy")
    recall, precision, F1, bal_acc = calc_accuracy(model, x_test, y_test, N, nb_classes)
    print("\tdone")

    model.save_weights(model_weights_fn)

    return history, recall, precision, F1, bal_acc


def train(step, data_dir, N, nb_classes, model, diagnostics_fn, model_weights_fn):
    '''
    trains continuously until force stopped
    '''
    file_list = os.listdir(data_dir + "train/")
    num_files = len(file_list)
    while True:
        #pick a random augmented dataset from parsed_dir
        i = np.random.randint(num_files)
        train_f = file_list[i]
        print("training step: " + str(step) + "\ttraining file: " + train_f)

        # grab data for training
        print("\tgrabbing data")
        stem = grab_data(data_dir, train_f, N, nb_classes)
        print("\tdone")

        # train
        history, recall, precision, F1, bal_acc = train_step(model, stem, model_weights_fn)
        print("recall = {}, precision = {}, F1 = {}, bal_acc = {}".format(recall, \
                        precision, F1, bal_acc))
        # record training step results
        loss     = history.history['loss'][0]
        val_loss = history.history['val_loss'][0]
        acc      = history.history['acc'][0]
        val_acc  = history.history['val_acc'][0]
        with open(diagnostics_fn, "a") as f:
             f.write('{:15d}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\n'.format(\
                step ,loss,acc,val_loss,val_acc,recall, precision, F1, bal_acc))
        step += 1
