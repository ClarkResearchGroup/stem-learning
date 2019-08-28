from detect_circles import extract_defects
from scipy.misc import imread
import matplotlib.pyplot as plt
from image_parse import process_data
from center_diff import detect_diff, get_center
import pickle

arg = "08_23_18"
parent_dir = "../../"
raw_data_dir = parent_dir + "data/raw_data/01_8Mx_1K_0.019876nm_RGB_PPT/"

'''
(input_img, label_img) = process_data(raw_data_dir + "input.png", raw_data_dir + "label2.png", tol=0.5, diff=False)
label_eval = imread("pngs_08_23_18/01_8Mx_1K_0.019876nm_RGB_PPT_prediction.png")

image_list, center_list, conv_label_img, num_defects_list, new_label_img = extract_defects(input_img, label_img, threshold=0.1, num_convs=3, radius=10, l=50)
image_list2, center_list2, conv_label_img2, num_defects_list2, new_label_img2 = extract_defects(input_img, label_eval, threshold=0.1, num_convs=3, radius=10, l=50)

pickle.dump(center_list, open("center_list.p", 'wb'))
pickle.dump(center_list2, open("center_list2.p", 'wb'))
'''
#evals_center = pickle.load(open("center_list2.p", 'rb'))
#label_center = pickle.load(open("center_list.p", 'rb'))

raw_data_dirs = ["01_8Mx_1K_0.019876nm_RGB_PPT",  "02_10Mx_1K_0.015806nm_RESIZED",\
        "03_8Mx_1K_0.019876nm_RGB_PPT"]

data_dir = raw_data_dirs[1]

label_file = raw_data_dir + "label2.png"
evals_file = "pngs_" + arg +"/" + data_dir + "_prediction.png"

save_fn = arg + "_" + data_dir + ".png"

label_center = get_center(label_file)
evals_center = get_center(evals_file)

match_list, evals_cen, label_cen = detect_diff(label_center, evals_center, 7.5)

(x, y) = zip(*match_list)if len(match_list) > 0 else ([],[])
(x2, y2) = zip(*evals_cen)if len(evals_cen) > 0 else ([],[])
(x3, y3) = zip(*label_cen) if len(label_cen) > 0 else ([],[])

fig = plt.figure(frameon=False)
fig.set_size_inches(10,10)
ax = plt.Axes(fig, [0.,0,1,1])
ax.set_axis_off()
fig.add_axes(ax)
ax.scatter(x,y, label='matched')
ax.scatter(x2,y2, color='r',label='false-positives')
ax.scatter(x3,y3, color='k',label='missed')
ax.set_ylim([1024,0])
fig.savefig(save_fn)

#plt.legend(loc='best')

'''
(x, y, z) = zip(*center_list)
(x2, y2, z2) = zip(*center_list2)
#plt.imshow(label_img)
plt.scatter(y, x)
#plt.imshow(label_img2)
plt.scatter(y2, x2)

plt.imshow(label_img)
plt.figure()
plt.imshow(conv_label_img)
plt.figure()
plt.imshow(new_label_img)
'''
plt.show()


