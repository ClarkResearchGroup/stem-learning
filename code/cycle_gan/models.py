import tensorflow as tf

class InstanceNormalization(tf.keras.layers.Layer):
  """Instance Normalization Layer (https://arxiv.org/abs/1607.08022)."""

  def __init__(self, epsilon=1e-5):
    super(InstanceNormalization, self).__init__()
    self.epsilon = epsilon

  def build(self, input_shape):
    self.scale = self.add_weight(
        name='scale',
        shape=input_shape[-1:],
        initializer=tf.random_normal_initializer(1., 0.02),
        trainable=True)

    self.offset = self.add_weight(
        name='offset',
        shape=input_shape[-1:],
        initializer='zeros',
        trainable=True)

  def call(self, x):
    mean, variance = tf.nn.moments(x, axes=[1, 2], keepdims=True)
    inv = tf.math.rsqrt(variance + self.epsilon)
    normalized = (x - mean) * inv
    return self.scale * normalized + self.offset


'''
  2D Reflection Padding
  Attributes:
    - padding: (padding_width, padding_height) tuple
'''
class ReflectionPadding2D(tf.keras.layers.Layer):
    def __init__(self, padding=(1, 1), **kwargs):
        self.padding = tuple(padding)
        super(ReflectionPadding2D, self).__init__(**kwargs)

    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[1] + 2 * self.padding[0], input_shape[2] + 2 * self.padding[1], input_shape[3])

    def call(self, input_tensor, mask=None):
        padding_width, padding_height = self.padding
        return tf.pad(input_tensor, [[0,0], [padding_height, padding_height], [padding_width, padding_width], [0,0] ], 'REFLECT')


def conv2d(dim, ks, s, padding="valid", activation=None):
  initializer = tf.keras.initializers.TruncatedNormal(stddev=0.02)
  return tf.keras.layers.Conv2D(dim, ks, s, padding=padding, activation=activation, kernel_initializer=initializer)

def conv2dtranspose(dim, ks, s, padding="valid"):
  initializer = tf.keras.initializers.TruncatedNormal(stddev=0.02)
  return tf.keras.layers.Conv2DTranspose(dim, ks, s, padding=padding, kernel_initializer=initializer)

def generator_resnet(gf_dim, ic=1, oc=1):

  def residule_block(x, dim, ks=3, s=1):
    p = int((ks - 1) / 2)
    y = ReflectionPadding2D(padding=(p,p))(x)
    y = InstanceNormalization()(conv2d(dim, ks, s, padding='valid')(y))
    y = ReflectionPadding2D(padding=(p,p))(tf.keras.layers.ReLU()(y))
    y = InstanceNormalization()(conv2d(dim, ks, s, padding='valid')(y))
    return tf.keras.layers.Add()([x, y])

  # Justin Johnson's model from https://github.com/jcjohnson/fast-neural-style/
  # The network with 9 blocks consists of: c7s1-32, d64, d128, R128, R128, R128,
  # R128, R128, R128, R128, R128, R128, u64, u32, c7s1-3
  image = tf.keras.layers.Input(shape=[None, None, ic])       

  c0 = ReflectionPadding2D(padding=(3,3))(image)
  c1 = tf.keras.layers.ReLU()(InstanceNormalization()(conv2d(gf_dim  , 7, 1, padding='valid')(c0)))
  c2 = tf.keras.layers.ReLU()(InstanceNormalization()(conv2d(gf_dim*2, 3, 2, padding='same')(c1)))
  c3 = tf.keras.layers.ReLU()(InstanceNormalization()(conv2d(gf_dim*4, 3, 2, padding='same')(c2)))


  # define G network with 9 resnet blocks
  r1 = residule_block(c3, gf_dim*4)
  r2 = residule_block(r1, gf_dim*4)
  r3 = residule_block(r2, gf_dim*4)
  r4 = residule_block(r3, gf_dim*4)
  r5 = residule_block(r4, gf_dim*4)
  r6 = residule_block(r5, gf_dim*4)
  r7 = residule_block(r6, gf_dim*4)
  r8 = residule_block(r7, gf_dim*4)
  r9 = residule_block(r8, gf_dim*4)

  d1 = conv2dtranspose(gf_dim*2, 3, 2, padding='same')(r9)
  d1 = tf.keras.layers.ReLU()(InstanceNormalization()(d1))
  d2 = conv2dtranspose(gf_dim, 3, 2, padding='same')(d1)
  d2 = tf.keras.layers.ReLU()(InstanceNormalization()(d2))
  d2 = ReflectionPadding2D(padding=(3,3))(d2)
  pred = conv2d(oc, 7, 1, padding='valid', activation='tanh')(d2)

  return tf.keras.Model(inputs=image, outputs=pred)


def downsample(filters, size, norm_type='batchnorm', apply_norm=True):
  """Downsamples an input.

  Conv2D => Batchnorm => LeakyRelu

  Args:
    filters: number of filters
    size: filter size
    norm_type: Normalization type; either 'batchnorm' or 'instancenorm'.
    apply_norm: If True, adds the batchnorm layer

  Returns:
    Downsample Sequential Model
  """
  initializer = tf.random_normal_initializer(0., 0.02)

  result = tf.keras.Sequential()
  result.add(
      tf.keras.layers.Conv2D(filters, size, strides=2, padding='same',
                             kernel_initializer=initializer, use_bias=False))

  if apply_norm:
    if norm_type.lower() == 'batchnorm':
      result.add(tf.keras.layers.BatchNormalization())
    elif norm_type.lower() == 'instancenorm':
      result.add(InstanceNormalization())

  result.add(tf.keras.layers.LeakyReLU())

  return result


def upsample(filters, size, norm_type='batchnorm', apply_dropout=False):
  """Upsamples an input.

  Conv2DTranspose => Batchnorm => Dropout => Relu

  Args:
    filters: number of filters
    size: filter size
    norm_type: Normalization type; either 'batchnorm' or 'instancenorm'.
    apply_dropout: If True, adds the dropout layer

  Returns:
    Upsample Sequential Model
  """

  initializer = tf.random_normal_initializer(0., 0.02)

  result = tf.keras.Sequential()
  result.add(
      tf.keras.layers.Conv2DTranspose(filters, size, strides=2,
                                      padding='same',
                                      kernel_initializer=initializer,
                                      use_bias=False))

  if norm_type.lower() == 'batchnorm':
    result.add(tf.keras.layers.BatchNormalization())
  elif norm_type.lower() == 'instancenorm':
    result.add(InstanceNormalization())

  if apply_dropout:
    result.add(tf.keras.layers.Dropout(0.5))

  result.add(tf.keras.layers.ReLU())

  return result


def unet_generator(input_channels, output_channels, norm_type='batchnorm'):
  """Modified u-net generator model (https://arxiv.org/abs/1611.07004).

  Args:
    output_channels: Output channels
    norm_type: Type of normalization. Either 'batchnorm' or 'instancenorm'.

  Returns:
    Generator model
  """

  down_stack = [
      downsample(64, 4, norm_type, apply_norm=False),  # (bs, 128, 128, 64)
      downsample(128, 4, norm_type),  # (bs, 64, 64, 128)
      downsample(256, 4, norm_type),  # (bs, 32, 32, 256)
      downsample(512, 4, norm_type),  # (bs, 16, 16, 512)
      downsample(512, 4, norm_type),  # (bs, 8, 8, 512)
      downsample(512, 4, norm_type),  # (bs, 4, 4, 512)
      downsample(512, 4, norm_type),  # (bs, 2, 2, 512)
      downsample(512, 4, norm_type),  # (bs, 1, 1, 512)
  ]

  up_stack = [
      upsample(512, 4, norm_type, apply_dropout=True),  # (bs, 2, 2, 1024)
      upsample(512, 4, norm_type, apply_dropout=True),  # (bs, 4, 4, 1024)
      upsample(512, 4, norm_type, apply_dropout=True),  # (bs, 8, 8, 1024)
      upsample(512, 4, norm_type),  # (bs, 16, 16, 1024)
      upsample(256, 4, norm_type),  # (bs, 32, 32, 512)
      upsample(128, 4, norm_type),  # (bs, 64, 64, 256)
      upsample(64, 4, norm_type),  # (bs, 128, 128, 128)
  ]

  initializer = tf.random_normal_initializer(0., 0.02)
  last = tf.keras.layers.Conv2DTranspose(
      output_channels, 4, strides=2,
      padding='same', kernel_initializer=initializer,
      activation='tanh')  # (bs, 256, 256, 3)

  concat = tf.keras.layers.Concatenate()

  inputs = tf.keras.layers.Input(shape=[None, None, input_channels])
  x = inputs

  # Downsampling through the model
  skips = []
  for down in down_stack:
    x = down(x)
    skips.append(x)

  skips = reversed(skips[:-1])

  # Upsampling and establishing the skip connections
  for up, skip in zip(up_stack, skips):
    x = up(x)
    x = concat([x, skip])

  x = last(x)

  return tf.keras.Model(inputs=inputs, outputs=x)


def discriminator(norm_type='batchnorm', ic=1):
  """PatchGan discriminator model (https://arxiv.org/abs/1611.07004).

  Args:
    norm_type: Type of normalization. Either 'batchnorm' or 'instancenorm'.
    target: Bool, indicating whether target image is an input or not.

  Returns:
    Discriminator model
  """

  initializer = tf.random_normal_initializer(0., 0.02)

  inp = tf.keras.layers.Input(shape=[None, None, ic], name='input_image')
  x = inp

  down1 = downsample(64, 4, norm_type, False)(x)  # (bs, 128, 128, 64)
  down2 = downsample(128, 4, norm_type)(down1)  # (bs, 64, 64, 128)
  down3 = downsample(256, 4, norm_type)(down2)  # (bs, 32, 32, 256)

  zero_pad1 = tf.keras.layers.ZeroPadding2D()(down3)  # (bs, 34, 34, 256)
  conv = tf.keras.layers.Conv2D(
      512, 4, strides=1, kernel_initializer=initializer,
      use_bias=False)(zero_pad1)  # (bs, 31, 31, 512)

  if norm_type.lower() == 'batchnorm':
    norm1 = tf.keras.layers.BatchNormalization()(conv)
  elif norm_type.lower() == 'instancenorm':
    norm1 = InstanceNormalization()(conv)

  leaky_relu = tf.keras.layers.LeakyReLU()(norm1)

  zero_pad2 = tf.keras.layers.ZeroPadding2D()(leaky_relu)  # (bs, 33, 33, 512)

  last = tf.keras.layers.Conv2D(
      1, 4, strides=1,
      kernel_initializer=initializer, activation='sigmoid')(zero_pad2)  # (bs, 30, 30, 1)

  return tf.keras.Model(inputs=inp, outputs=last)