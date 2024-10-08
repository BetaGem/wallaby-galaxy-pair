# A simplified 3D version of the 'find_peak' function in photutils.detection
# https://photutils.readthedocs.io/en/stable/api/photutils.detection.find_peaks.html

import warnings
import numpy as np
from astropy.table import QTable

def find_peaks(data, threshold, box_size=3, footprint=None, mask=None,
               border_width=None, npeaks=np.inf):
    
    from scipy.ndimage import maximum_filter
    
    data = np.asanyarray(data)
    
    # remove NaN values to avoid runtime warnings
    nan_mask = np.isnan(data)
    if np.any(nan_mask):
        data = np.copy(data)  # ndarray
        data[nan_mask] = np.nanmin(data)

    if footprint is not None:
        data_max = maximum_filter(data, footprint=footprint, mode='nearest')
    else:
        data_max = maximum_filter(data, size=box_size, mode='nearest')
    
    peak_goodmask = (data == data_max)  # good pixels are True

    if mask is not None:
        mask = np.asanyarray(mask)
        if data.shape != mask.shape:
            raise ValueError('data and mask must have the same shape')
        peak_goodmask = np.logical_and(peak_goodmask, ~mask)

    if border_width is not None:
        for i in range(peak_goodmask.ndim):
            peak_goodmask = peak_goodmask.swapaxes(0, i)
            peak_goodmask[:border_width] = False
            peak_goodmask[-border_width:] = False
            peak_goodmask = peak_goodmask.swapaxes(0, i)

    peak_goodmask = np.logical_and(peak_goodmask, (data > threshold))
    
    z_peaks, y_peaks, x_peaks = peak_goodmask.nonzero()
    peak_values = data[z_peaks, y_peaks, x_peaks]

    nxpeaks = len(x_peaks)
    if nxpeaks > npeaks:
        idx = np.argsort(peak_values)[::-1][:npeaks]
        x_peaks = x_peaks[idx]
        y_peaks = y_peaks[idx]
        z_peaks = z_peaks[idx]
        peak_values = peak_values[idx]

    # construct the output table
    meta = {'version': 0.0}
    colnames = ['z_peak', 'y_peak', 'x_peak', 'peak_value']
    coldata = [z_peaks, y_peaks, x_peaks, peak_values]
    table = QTable(coldata, names=colnames, meta=meta)

    return table