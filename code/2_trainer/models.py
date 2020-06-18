from tensorflow.keras.layers import Input, MaxPooling2D, UpSampling2D, Add
from tensorflow.keras.layers import Dropout, Activation, Reshape, concatenate
from tensorflow.keras.layers import Conv2D, Conv2DTranspose
from tensorflow.keras.models import Model
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.initializers import glorot_uniform

def model_lattice(input_img, N, k_fac, nb_classes):
    '''
    creates a convolution net with (assuming k_fac = 4)
    input (128 x 128 x 1) -> convolution (128 x 128 x  8) -> maxpooling ( 64 x  64 x  8) ->
                             convolution ( 64 x  64 x 16) -> maxpooling ( 32 x  32 x 16) ->
                             convolution ( 32 x  32 x 32) -> maxpooling ( 16 x  16 x 32) ->
                             convolution ( 16 x  16 x 32) -> upsampling ( 32 x  32 x 32) ->
                             convolution ( 32 x  32 x 16) -> upsampling ( 64 x  64 x 16) ->
                             convolution ( 64 x  64 x  8) -> upsampling (128 x 128 x  8) ->
                convolution (128 x 128 x  3) -> convolution (128 x 128 x  3) ->
    output (128 x 128 x nb_classes)
    '''
    x = Conv2D(2*k_fac, (3, 3), activation='relu', padding='same')(input_img)
    x = MaxPooling2D((2, 2), padding='same')(x)
    x = Conv2D(4*k_fac, (3, 3), activation='relu', padding='same')(x)
    x = MaxPooling2D((2, 2), padding='same')(x)
    x = Conv2D(8*k_fac, (3, 3), activation='relu', padding='same')(x)
    x = MaxPooling2D((2, 2), padding='same')(x)

    x = Conv2D(8*k_fac, (3, 3), activation='relu', padding='same')(x)
    x = UpSampling2D((2, 2))(x)
    x = Conv2D(4*k_fac, (3, 3), activation='relu', padding='same')(x)
    x = UpSampling2D((2, 2))(x)
    x = Conv2D(2*k_fac, (3, 3), activation='relu', padding='same')(x)
    x = UpSampling2D((2, 2))(x)

    x = Conv2D(nb_classes, (3, 3), activation = 'linear', padding='same')(x)
    x = Conv2D(nb_classes, (1, 1), activation = 'linear', padding='same')(x)
    x = Reshape((N*N, nb_classes))(x)
    output = Activation('softmax')(x)

    return Model(input_img, output)


def model_unet(input_img, N, n_filters = 16, dropout = 0.1, batchnorm = True, nb_classes = 2):

    def conv2d_block(input_tensor, n_filters, kernel_size = 3, batchnorm = True):
        '''Function to add 2 convolutional layers with the parameters passed to it'''
        # first layer
        x = Conv2D(filters = n_filters, kernel_size = (kernel_size, kernel_size),\
                  kernel_initializer = 'he_normal', padding = 'same')(input_tensor)
        if batchnorm:
            x = BatchNormalization()(x)
        x = Activation('relu')(x)

        # second layer
        x = Conv2D(filters = n_filters, kernel_size = (kernel_size, kernel_size),\
                  kernel_initializer = 'he_normal', padding = 'same')(x)
        if batchnorm:
            x = BatchNormalization()(x)
        x = Activation('relu')(x)

        return x

    # Contracting Path
    c1 = conv2d_block(input_img, n_filters * 1, kernel_size = 3, batchnorm = batchnorm)
    p1 = MaxPooling2D((2, 2))(c1)
    p1 = Dropout(dropout)(p1)

    c2 = conv2d_block(p1, n_filters * 2, kernel_size = 3, batchnorm = batchnorm)
    p2 = MaxPooling2D((2, 2))(c2)
    p2 = Dropout(dropout)(p2)

    c3 = conv2d_block(p2, n_filters * 4, kernel_size = 3, batchnorm = batchnorm)
    p3 = MaxPooling2D((2, 2))(c3)
    p3 = Dropout(dropout)(p3)

    c4 = conv2d_block(p3, n_filters * 8, kernel_size = 3, batchnorm = batchnorm)
    p4 = MaxPooling2D((2, 2))(c4)
    p4 = Dropout(dropout)(p4)

    c5 = conv2d_block(p4, n_filters = n_filters * 16, kernel_size = 3, batchnorm = batchnorm)

    # Expansive Path
    u6 = Conv2DTranspose(n_filters * 8, (3, 3), strides = (2, 2), padding = 'same')(c5)
    u6 = concatenate([u6, c4])
    u6 = Dropout(dropout)(u6)
    c6 = conv2d_block(u6, n_filters * 8, kernel_size = 3, batchnorm = batchnorm)

    u7 = Conv2DTranspose(n_filters * 4, (3, 3), strides = (2, 2), padding = 'same')(c6)
    u7 = concatenate([u7, c3])
    u7 = Dropout(dropout)(u7)
    c7 = conv2d_block(u7, n_filters * 4, kernel_size = 3, batchnorm = batchnorm)

    u8 = Conv2DTranspose(n_filters * 2, (3, 3), strides = (2, 2), padding = 'same')(c7)
    u8 = concatenate([u8, c2])
    u8 = Dropout(dropout)(u8)
    c8 = conv2d_block(u8, n_filters * 2, kernel_size = 3, batchnorm = batchnorm)

    u9 = Conv2DTranspose(n_filters * 1, (3, 3), strides = (2, 2), padding = 'same')(c8)
    u9 = concatenate([u9, c1])
    u9 = Dropout(dropout)(u9)
    c9 = conv2d_block(u9, n_filters * 1, kernel_size = 3, batchnorm = batchnorm)

    x = Conv2D(nb_classes, (1, 1), activation = 'linear', padding='same')(c9)
    x = Reshape((N*N, nb_classes))(x)
    output = Activation('softmax')(x)

    return Model(input_img, output)



def res_helper(Input, filters, stride=1, basic=True, f=3):
    ks = (3,3) if basic else (1,1)
    ##### MAIN PATH #####
    # First component of main path
    X = Conv2D(filters=filters[0], kernel_size=ks, strides=(stride, stride), padding='same',\
            kernel_initializer=glorot_uniform(seed=0))(Input)
    X = BatchNormalization(axis=3)(X)
    X = Activation('relu')(X)

    # Second component of main path
    if not basic:
        X = Conv2D(filters=filters[1], kernel_size=(f, f), strides=(1, 1), padding='same',\
                kernel_initializer=glorot_uniform(seed=0))(X)
        X = BatchNormalization(axis=3)(X)
        X = Activation('relu')(X)

    # Third component of main path
    X = Conv2D(filters=filters[-1], kernel_size=ks, strides=(1, 1), padding='same',\
            kernel_initializer=glorot_uniform(seed=0))(X)
    X = BatchNormalization(axis=3)(X)

    X = Activation('relu')(X)
    return X


def res_conv_block(X, filters, stage, block, basic=True, f=3, conv=False, dropout=0.1, stride=1):
    X_shortcut = X
    res = res_helper(X, filters, stride, basic, f)
    if not conv:
        res = Add()([res, X_shortcut])
        Y = MaxPooling2D((2, 2))(res)
    else:
        Y = MaxPooling2D((2, 2))(res) if stride == 1 else res
        X_shortcut = Conv2D(filters=filters[-1], kernel_size=(1, 1), strides=(2, 2), \
                padding='valid', kernel_initializer=glorot_uniform(seed=0))(X_shortcut)
        X_shortcut = BatchNormalization(axis=3)(X_shortcut)
        Y = Add()([Y, X_shortcut])

    Y = Dropout(dropout)(Y)
    return Y, res


def res_deconv_block(X, filters, stage, block, basic=True, f=3, dropout=0.1):
    Y_shortcut = Conv2DTranspose(filters[0], (3, 3), strides = (2, 2), padding = 'same')(X)
    Y = res_helper(Y_shortcut, filters, stride=1, basic=basic, f=f)
    Y = Add()([Y, Y_shortcut])
    return Y

def res_upsamp_block(c, r, n_filters, k, n_fac, l1, l2, basic, f, dropout, stride):
    u = res_deconv_block(c, [n_filters*n_fac]*k, l1, l2, basic, f, dropout)
    u = concatenate([u, r])
    u = Dropout(dropout)(u)
    return res_helper(u, [n_filters*n_fac]*k, stride=stride, basic=basic, f=f)



def model_resunet(input_img, N, n_filters, nb_classes, basic=True, dropout=0.1, stride=1, f=3, conv=True):
    """Creates a Deep Learning model for defect identification"""
    # Contracting Path
    k = 2 if basic else 3
    c1, r1 = res_conv_block(input_img, [n_filters*1]*k, '1', 'a', basic, f, conv, dropout, stride)
    c2, r2 = res_conv_block(c1,        [n_filters*2]*k, '1', 'b', basic, f, conv, dropout, stride)
    c3, r3 = res_conv_block(c2,        [n_filters*4]*k, '1', 'c', basic, f, conv, dropout, stride)
    c4, r4 = res_conv_block(c3,        [n_filters*8]*k, '1', 'd', basic, f, conv, dropout, stride)

    c5 = res_helper(c4, [n_filters*16]*k, stride=1, basic=basic, f=f)

    # Expansive Path
    c6 = res_upsamp_block(c5, r4, n_filters, k, 8, '2', 'a', basic, f, dropout, stride)
    c7 = res_upsamp_block(c6, r3, n_filters, k, 4, '2', 'b', basic, f, dropout, stride)
    c8 = res_upsamp_block(c7, r2, n_filters, k, 2, '2', 'c', basic, f, dropout, stride)
    c9 = res_upsamp_block(c8, r1, n_filters, k, 1, '2', 'd', basic, f, dropout, stride)

    x = Conv2D(nb_classes, (1, 1), activation = 'linear', padding='same')(c9)
    x = Reshape((N*N, nb_classes))(x)
    output = Activation('softmax')(x)

    return Model(input_img, output)
