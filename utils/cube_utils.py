import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
from photutils.segmentation import make_2dgaussian_kernel


def freq_smooth(cube, bin_size=0, smooth=0):
    '''
    bin or smooth the data cube along the frequency axis.
    '''
    if bin_size:    # bin
        shape = cube.shape
        cube_new = np.zeros((shape[0]//bin_size, shape[1], shape[2]))
        for freq in range(cube.shape[0]//bin_size):
            cube_new[freq] = np.sum(cube[freq*bin_size : (freq+1)*bin_size], axis=0)
    elif smooth:    # smooth
        cube_new = np.zeros_like(cube)
        kernel = Gaussian1DKernel(smooth)
        for i in range(cube.shape[1]):
            for j in range(cube.shape[2]):
                cube_new[:,i,j] = convolve(cube[:,i,j], kernel)
    else: cube_new = cube  # do nothing
    return cube_new


def extend(array, num=1, smooth=False):
    '''
    enlarge and smooth an image.
    '''
    shape = np.shape(array)
    dim = len(shape)
    if dim == 3: extend_array = np.zeros((shape[0], num*shape[1], num*shape[2]))
    elif dim==2: extend_array = np.zeros((num*shape[0], num*shape[1]))
    else: raise ValueError("Dimension must be 2 or 3.")
    for i in range(shape[-2]):
        for j in range(shape[-1]):
            if dim == 3: 
                for zz in range(shape[0]):
                    extend_array[zz, i*num:(i+1)*num, j*num:(j+1)*num] = array[zz,i,j]
            elif dim==2: 
                extend_array[i*num:(i+1)*num, j*num:(j+1)*num] = array[i][j]
    # 3d smoothing is time-consuming
    if smooth:
        kernel2d = make_2dgaussian_kernel(fwhm=num, size=4*num+1)
        if dim == 2: 
            extend_array = convolve(extend_array, kernel2d)
        elif dim == 3: 
            for freq in range(len(extend_array)):
                extend_array[freq] = convolve(extend_array[freq], kernel2d)
    return extend_array