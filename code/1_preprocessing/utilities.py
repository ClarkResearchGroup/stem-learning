import shutil, os

def make_folder(data_folder_with_labels, input_folder, target_folder):
    shutil.copytree(data_folder_with_labels, target_folder, dirs_exist_ok=True)
    for img_file, save_dir in zip(os.listdir(input_folder), os.listdir(target_folder)):
        full_src = os.path.join(input_folder, img_file)
        full_dst = os.path.join(target_folder, save_dir, "input.tiff")
        shutil.copy(full_src, full_dst)