from keras.models import model_from_json
import pickle
import numpy as np
import matplotlib.pyplot as plt


N = 128
def load_data(data_file):
    data = pickle.load(open(data_file, 'rb'))

    images = np.asarray([data[i][0] for i in range(len(data))], dtype=np.float32)
    labels = np.asarray([data[i][1] for i in range(len(data))], dtype=np.float32)

    images = np.reshape(images, [-1,N,N,1])
    labels = np.reshape(labels, [-1,N*N,2])

    return (images, labels)


def model_load(model_fn, model_weights_fn):
    with open(model_fn, 'r') as f:
            model = model_from_json(f.read())
    model.load_weights(model_weights_fn)
    return model


def predict(model, images):
    return model.predict_on_batch(images)

def plot_data(images, labels, predis, idx):
    img = np.reshape(images[idx], [N,N])
    lbl = np.reshape(labels[idx,:,1], [N,N])
    prd = np.reshape(predis[idx,:,1], [N,N])

    plt.imshow(img, cmap='gray', vmin=0, vmax=1)
    plt.title("input")
    plt.figure()
    plt.imshow(lbl, cmap='gray', vmin=0, vmax=1)
    plt.title("label")
    plt.figure()
    plt.imshow(prd, cmap='gray', vmin=0, vmax=1)
    plt.title("prediction")
    plt.show()


if __name__ == "__main__":
    arg = "10-24-2018+1"
    parent_dir = "/media/bkcgroup2/abid/stem-learning_experimental/"
    data_dir = parent_dir + "code/preprocessing/generator/parsed_data/test/"
    model_fn = parent_dir + "model/model" + arg + ".h5"
    model_weights_fn = parent_dir + "model/model" + arg + "_weigths.h5"
    test_data = "test_flipped_shrink_rotate_043.p"

    print "loading dataset"
    (images, labels) = load_data(data_dir + test_data)
    print "loading model"
    model = model_load(model_fn, model_weights_fn)
    print "predicting outputs"
    predis = predict(model, images)
    print "plotting data"
    idx = 0
    plot_data(images, labels, predis, idx)






