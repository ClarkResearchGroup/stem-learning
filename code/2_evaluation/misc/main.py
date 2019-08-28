from stitcher import make_prediction

arg = "08-23-2018+1"
#arg = "08-23-2018+1"
parent_dir = "../../"
model_fn         = parent_dir + "model/model" + arg + ".h5"
model_weights_fn = parent_dir + "model/model" + arg + "_weigths.h5"
Tol              = 0.5
Diff             = False
avg              = True
thresh           = 0.5

#for sampling
test_dir = "20181114_Sampling_tolerance_test/"
file_list = ["00_01_8Mx_0.019876nm_1024", "01_01_8Mx_0.039752nm_512",   \
        "02_01_8Mx_0.04969nm_409", "03_01_8Mx_0.059628nm_341", "04_01_8Mx_0.079504nm_256"]



#raw_data_dir = ["01_8Mx_1K_0.019876nm_RGB_PPT",  "02_10Mx_1K_0.015806nm_RESIZED",\
#        "03_8Mx_1K_0.019876nm_RGB_PPT"]
raw_data_dir = ["01_8Mx_1K_0.019876nm_RGB_PPT_512"]#, "02_10Mx_1K_0.015806nm_RGB_PPT_RESIZED_downsampled", "03_8Mx_1K_0.019876nm_RGB_PPT_downsampled"]
#"02_10Mx_1K_0.015806nm_RGB_PPT", "02_10Mx_1K_0.015806nm_RGB_PPT_RESIZED_downsampled",  "03_8Mx_1K_0.019876nm_RGB_PPT_downsampled", "01_8Mx_1K_0.019876nm_RGB_PPT_downsampled"]
for data_dir in raw_data_dir:
    print data_dir

    input_file = parent_dir + "data/raw_data/" + data_dir + "/input.png"
    label_file = parent_dir + "data/raw_data/" + data_dir + "/label2.png"
    save_dir   = parent_dir + "code/postprocessing/"
    prefix     = data_dir + "_"

    make_prediction(model_fn, model_weights_fn, input_file, label_file, Tol, Diff, avg, save_dir,\
            thresh, prefix, plot=False)

#for noise
#test_dir = "8bit_demo_file_Resize/"

#for file_fn in file_list:
'''
for i in range(1, 11):
    if i != 7:
        continue
    file_fn = "MoWTe2_incostem_60_36_1_" + str(i).zfill(3)
    print file_fn

    input_file = parent_dir + "data/simulated_data/input_image/"+ file_fn + ".tif"
    label_file = parent_dir + "data/simulated_data/label_image/W_"+ file_fn + ".tif"
    save_dir   = parent_dir + "code/postprocessing/"
    prefix     = file_fn + "_"

    make_prediction(model_fn, model_weights_fn, input_file, label_file, Tol, Diff, avg, save_dir,\
            thresh, prefix, plot=False)
'''
