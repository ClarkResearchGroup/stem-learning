#!/usr/bin/env python

"""Functions for downloading and reading MNIST data."""
import numpy as np
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


class DataSet(object):

    def __init__(self, images, labels):
        assert np.shape(images)[0] == np.shape(labels)[0], (
                "images.shape: %s labels.shape: %s" % (images.shape, labels.shape))
        self._num_examples = images.shape[0]
        self._images = images
        self._labels = labels
        self._epochs_completed = 0
        self._index_in_epoch = 0

    @property
    def images(self):
        return self._images

    @property
    def labels(self):
        return self._labels

    @property
    def num_examples(self):
        return self._num_examples


    @property
    def epochs_completed(self):
        return self._epochs_completed


    def next_batch(self, batch_size, shuff=True):
        """Return the next `batch_size` examples from this data set."""
        start = self._index_in_epoch
        self._index_in_epoch += batch_size
        if self._index_in_epoch > self._num_examples:
            # Finished epoch
            self._epochs_completed += 1
            # Shuffle the data
            if shuff:
                perm = np.arange(self._num_examples)
                np.random.shuffle(perm)
                self._images = self._images[perm]
                self._labels = self._labels[perm]
            # Start next epoch
            start = 0
            self._index_in_epoch = batch_size
            assert batch_size <= self._num_examples
        end = self._index_in_epoch
        return self._images[start:end], self._labels[start:end]


def read_data_set(data_file, N, nb_classes):
    class DataSets(object):
        pass
    data_sets = DataSets()
    data = pickle.load(open(data_file, 'rb'))
    data_images=np.asarray([data[i][0] for i in range(len(data))], dtype=np.float32)
    data_labels=np.asarray([data[i][1] for i in range(len(data))], dtype=np.float32)

    data_sets.data = DataSet(data_images, data_labels)
    data_sets.N = N
    data_sets.nb_classes = nb_classes

    return data_sets


stem_dict = {}
def grab_data(data_dir, train_f, N, nb_classes):
    '''
    creates a dictionary of DataSet objects
    '''
    if not train_f in stem_dict.keys():
        train_path = data_dir + train_f
        stem = read_data_set(train_path, N, nb_classes)
        try:
            stem_dict[train_f] = stem
        except:
            return stem
    return stem_dict[train_f]
