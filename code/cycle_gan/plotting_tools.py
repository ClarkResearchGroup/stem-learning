import numpy as np
import matplotlib.pyplot as plt

def get_discriminator_acc(discriminator, real_img, fake_img, threshold=0.5):
    
    logit_real = discriminator(real_img, training=False).numpy()
    logit_fake = discriminator(fake_img, training=False).numpy()

    p_real = 1/(1+np.exp(-logit_real))
    p_fake = 1/(1+np.exp(-logit_fake))

    vmin = -np.max([np.max(np.abs(p_real)), np.max(np.abs(p_fake))])
    vmax = -vmin

    plt.figure()
    plt.subplot(1,2,1)
    plt.title("real")
    plt.imshow(np.mean(p_real, axis=0)[:,:,0], cmap='bwr', vmin=0, vmax=1)
    plt.colorbar()

    plt.subplot(1,2,2)
    plt.title("fake")
    plt.imshow(np.mean(p_fake, axis=0)[:,:,0], cmap='bwr', vmin=0, vmax=1)
    plt.colorbar()
    plt.show()
    
    disc_real = p_real > threshold
    disc_fake = p_fake > threshold

    print(np.size(disc_fake), np.size(disc_real))

    TP = np.sum(disc_real)
    FP = np.size(disc_real) - TP
    FN = np.sum(disc_fake)
    TN = np.size(disc_fake) - FN
    print("TP: {}".format(TP))
    print("FP: {}".format(FP))
    print("FN: {}".format(FN))
    print("TN: {}".format(TN))
    print("recall:    {}".format(TP/(TP+FN)))
    print("precision: {}".format(TP/(TP+FP)))
    print("accuracy:  {}".format((TP + TN)/(TP + TN + FP + FN)))

    acc = float(TP + TN)/float(TP + TN + FP + FN)
    return acc

def generate_images(model_x, model_y, test_input):
  prediction = model_x(test_input)
  same       = model_y(prediction)
    
  plt.figure(figsize=(12, 18))

  display_list = [test_input[0], prediction[0], same[0]]
  title = ['Input Image', 'Predicted Image', 'Cycle Image']

  for i in range(3):
    plt.subplot(1, 3, i+1)
    plt.title(title[i])
    plt.imshow(display_list[i][:,:,0], cmap='gray')
    plt.axis('off')
  plt.show()

def generate_losses(gen_exp_losses,   gen_sim_losses,\
                    cycle_exp_losses, cycle_sim_losses,\
                    ident_exp_losses, ident_sim_losses,\
                    disc_sim_losses,  disc_exp_losses, epoch, \
                    last=None):
    ymin, ymax = 1e-2, 1e2
    x = list(range(epoch))
    N = len(x) if last is None or last > len(x) else last
    fig = plt.figure(figsize=(12,6))
   
    plt.subplot(1,2,1)
    plt.semilogy(x[-N:], disc_sim_losses[-N:], label='disc sim')
    plt.semilogy(x[-N:], disc_exp_losses[-N:], label='disc exp')
    plt.xlabel("epoch")
    plt.ylabel("loss")
    plt.legend(loc='best')

    plt.subplot(1,2,2)
    plt.semilogy(x[-N:], gen_exp_losses[-N:], label='gen_exp')
    plt.semilogy(x[-N:], gen_sim_losses[-N:], label='gen_sim')
    plt.semilogy(x[-N:], cycle_exp_losses[-N:], label='cycle_exp')
    plt.semilogy(x[-N:], cycle_sim_losses[-N:], label='cycle_sim')
    plt.semilogy(x[-N:], ident_exp_losses[-N:], label='ident_exp')
    plt.semilogy(x[-N:], ident_sim_losses[-N:], label='ident_sim')
    plt.xlabel("epoch")
    plt.ylabel("loss")
    plt.legend(loc="best")

    plt.show()

def generate_accuracies(disc_g_acc, disc_f_acc, epoch, last=None):
    x = list(range(epoch))
    N = len(x) if last is None else last
    plt.figure(figsize=(12,6))
    plt.plot(x[-N:], disc_g_acc[-N:], label="disc sim")
    plt.plot(x[-N:], disc_f_acc[-N:], label="disc_exp")
    plt.xlabel("epoch")
    plt.ylabel("accuracy")
    plt.ylim(0,1)
    plt.legend(loc='best')
    plt.show()
