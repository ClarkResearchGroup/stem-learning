# Deep learning models for defect identification in STEM images

## Quick Overview
This "stem-learning" repo contains 2 parts:
1. Generating and post-processing of simulated STEM images
2. Training fully convolutional network (FCN) models that can identify atomic defects like substitutions or vacancies

These models are used to identify point defects like Te substitutions or Se vacancies in WSe<sub>2-2x</sub>Te<sub>2x</sub> and have made sub-pm precision measurement of strain fields of single-atom defects possible. More discussion of the results can be found on this NanoLetter paper <https://pubs.acs.org/doi/abs/10.1021/acs.nanolett.0c00269>.

 
## Getting Started

The code can be acquired by cloning this repository to your computer, using the green "Clone or download" button, or by typing into the command line

```
git clone https://github.com/ClarkResearchGroup/stem-learning.git
```

See `/code/training_tutorial.ipynb` for a walkthrough of the preprocessing, training, and evaluation of the FCNs on your data.
There are also jupyter notebooks in each folder describing the functionality of each section.

It's recommended to use Google Colab for the model training/evaluation since it's not taking too much time. If you want to train/evaluate on your local machine, here's our build.

### Local build
Python 3, Keras 2.2.4, TensorFlow 1.3

If you have any questions, feel free to contact 

Abid Khan:     aakhan3@illinois.edu

Chia-Hao Lee: chiahao3@illinois.edu
