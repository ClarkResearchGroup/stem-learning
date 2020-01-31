from stitcher import *

'''
#arg = "08-23-2018+1"
arg     = sys.argv[1]
exp_idx = int(sys.argv[2]) - 1
labeled = int(sys.argv[4])

parent_dir = "../../"
exp_dirs   = ["01_8Mx_1K_0.019876nm_RGB_PPT/"  , "02_10Mx_1K_0.015806nm_RGB_PPT/", \
              "03_8Mx_1K_0.019876nm_RGB_PPT/"  , "04_17b_1K_8Mx_0.019876nm/", \
              "05_19b_1K_8Mx_0.019876nm/"      , "06_20b_1K_8Mx_0.019876nm/",\
              "07_23b_1K_8Mx_0.019876nm/"      , "08_25b_1K_8Mx_0.019876nm/",\
              "02_10Mx_1K_0.015806nm_RESIZED/"\
              ]

data_dir     = parent_dir + "data/"
dataset_dir  = "simulated_data/" if exp_idx < 0 else "raw_data/" + exp_dirs[exp_idx]
avg_path     = "avg/" if avg else "non_avg/"
datasave_dir = "simulation/" if exp_idx < 0 else exp_dirs[exp_idx]


model_fn         = parent_dir + "model/model" + arg + ".h5"
model_weights_fn = parent_dir + "model/model" + arg + "_weigths.h5"
input_file       = data_dir + dataset_dir + "input.png"
label_file      = data_dir + dataset_dir + "label.png" if labeled else input_file
Tol              = 1e-5 if exp_idx else 0.5
Diff             = exp_idx
avg              = int(sys.argv[3])
path             = "/home/aakhan3/" + arg + "/" + datasave_dir + avg_path

'''

parent_dir  = "/home/abid/Dropbox/Development/programs/stem-learning/"
results_dir = parent_dir + "results_v2/2Te/"
model_fn = results_dir + "model.json"
model_weights_fn = results_dir + "weights.h5"
data_dir = parent_dir + "data/WSeTe/simulated/0/"
input_file =  data_dir + "input.tiff"
label_file_list = [data_dir + "label_2Te.tiff"]#, data_dir + "label_TeSe.tif"]
print(input_file)
#label_file_list = [data_dir + "label_W.tif"]
Tol = 0.05
avg = 1
path = "./"
thresh = -1
make_prediction(model_fn, model_weights_fn, input_file, label_file_list, Tol, avg, path,\
        thresh=thresh,plot=True)




