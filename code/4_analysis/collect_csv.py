import os

parent_dir = "/u/skhan/stem-learning/"
lbl = "1vacancy"

save_dir = parent_dir

prefix = "20220705_MODEL_unet_ident_gen_fft"

results_exp_dir = os.listdir("{}results/exp/".format(parent_dir))
results_sim_dir = os.listdir("{}results/sim/".format(parent_dir))

results_exp_files = ["{}results/exp/{}/{}/evals/label_row.csv".format(parent_dir, idfn, lbl) \
                        for idfn in results_exp_dir if prefix in idfn]
results_sim_files = ["{}results/sim/{}/{}/evals/label_row.csv".format(parent_dir, idfn, lbl) \
                        for idfn in results_exp_dir if prefix in idfn]


results_files = results_exp_files + results_sim_files


label_exp_csv = "{}all_evals.csv".format(save_dir)

if os.path.exists(label_exp_csv):
    os.remove(label_exp_csv)

with open(label_exp_csv, "a") as f:
    for fn in results_files:
        try:
            with open(fn, "r") as r:
                f.write(r.read())
                f.write("\n")
        except Exception as e:
            print(e)
