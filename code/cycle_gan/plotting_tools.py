import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf

def get_discriminator_acc(discriminator, real_img, fake_img, threshold=0.5, with_logits=False):
    
    p_real = discriminator(real_img, training=False).numpy()
    p_fake = discriminator(fake_img, training=False).numpy()

    if with_logits:
      p_real = 1/(1+np.exp(-p_real))
      p_fake = 1/(1+np.exp(-p_fake))

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

def get_logf2(img_arr):
  fine_size = img_arr.shape[1]
  fft_list = tf.signal.fftshift(tf.signal.fft2d(tf.cast(tf.reshape(img_arr, [-1, fine_size, fine_size]), tf.complex64)))
  re_list, im_list = tf.math.real(fft_list), tf.math.imag(fft_list)
  f2_list = tf.math.multiply(re_list, re_list) + tf.math.multiply(im_list, im_list)
  return tf.reshape(tf.math.log(tf.clip_by_value(f2_list, 1e-36, 1e36))/tf.math.log(10.), [-1, fine_size, fine_size, 1])

def generate_images(model_x, model_y, test_input):
  prediction = model_x(test_input, training=True)
  same       = model_y(prediction, training=True)
    
  plt.figure(figsize=(12, 18))

  display_list = [test_input[0], prediction[0], same[0]]
  display_list_fft = [get_logf2(x)[0] for x in display_list]
  title = ['Input Image', 'Predicted Image', 'Cycle Image']

  for i in range(3):
    plt.subplot(2, 3, i+1)
    plt.title(title[i])
    plt.imshow(display_list[i][:,:,0], cmap='gray')
    plt.axis('off')

    plt.subplot(2, 3, i+4)
    plt.imshow(display_list_fft[i][:,:,0], cmap='gray')
    plt.axis('off')
  plt.show()

def generate_losses(gen_exp_losses,   gen_sim_losses,  \
                    cycle_exp_losses, cycle_sim_losses,\
                    ident_exp_losses, ident_sim_losses,\
                    disc_sim_losses,  disc_exp_losses, \
                    disc_sim_fft_losses,  disc_exp_fft_losses, \
                    epoch, last=None):
    ymin, ymax = 1e-2, 1e2
    x = list(range(epoch))
    N = len(x) if last is None or last > len(x) else last
    fig = plt.figure(figsize=(12,6))
   
    plt.subplot(1,2,1)
    plt.semilogy(x[-N:], disc_sim_losses[-N:], label='disc sim')
    plt.semilogy(x[-N:], disc_exp_losses[-N:], label='disc exp')
    plt.semilogy(x[-N:], disc_sim_fft_losses[-N:], label='disc sim_fft')
    plt.semilogy(x[-N:], disc_exp_fft_losses[-N:], label='disc exp_fft')
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
