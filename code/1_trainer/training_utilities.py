import os
import numpy as np
from models import *
from input_data import *


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
            f.write('{:>15}\t{:>15}\t{:>15}\t{:>15}\t{:>15}\n'.format(\
                    'step','loss','acc','val_loss','val_acc'))
    return step


def train_step(model, stem, epochs=1, batch_size=32):
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
    model.save_weights(model_weights_fn)

    return history


def train(step, data_dir, N, nb_classes, model, diagnostics_fn):
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
        stem = grab_data(data_dir, train_f, N, nb_classes)

        # train
        history = train_step(model, stem)

        # record training step results
        loss     = history.history['loss'][0]
        val_loss = history.history['val_loss'][0]
        acc      = history.history['acc'][0]
        val_acc  = history.history['val_acc'][0]
        with open(diagnostics_fn, "a") as f:
            f.write('{:15d}\t{:.15e}\t{:.15e}\t{:.15e}\t{:.15e}\n'.format(\
                    step ,loss,acc,val_loss,val_acc))
        step += 1
