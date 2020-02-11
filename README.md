# Deep learning models for defect identification in STEM images

## Quick Overview
This "stem-learning" repo contains 2 parts:
1. Generating and post-processing of simulated STEM images
2. Training FCN models that can identify atomic defects like substitutions or vacancies

These models are used to identify point defects like Te substitutions or Se vacancies in WSe<sub>2-2x</sub>Te<sub>2x</sub> and have made sub-pm precision measurement of strain fields of single-atom defects possible. More discussion of the results can be found on arXiv <https://arxiv.org/abs/2001.08233>.

 
## Getting Started

The code can be acquired by cloning this repository to your computer, using the green "Clone or download" button, or by typing into the command line

```
git clone https://github.com/ClarkResearchGroup/stem-learning.git
```

It's recommended to use Google Colab for the model training/evaluation since it's not taking too much time. If you want to train/evaluate on your local machine, here's our build.

### Local build
Python 3, Keras 2.2.4, TensorFlow 1.3
