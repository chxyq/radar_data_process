 #*- coding: UTF-8 -*-
 #encoding:utf-8
 
import numpy as np
from numpy import fft
import parameter

n_samples = parameter.n_samples
n_chirps = parameter.n_chirps
numTx = parameter.numTx
numRx = parameter.numRx
calibrationInterp = parameter.calibrationInterp

def hann_local(len):
    win = (np.arange(1 , len + 1)) / float(len + 1)
    win = 0.5 - 0.5 * np.cos(2 * np.pi * win)
    return win

def Doppler_FFT(rangeFFT_in):
    data_to_Doppler = np.zeros((n_samples * calibrationInterp, numTx * numRx, n_chirps), dtype=complex)
    Doppler_out = np.zeros((n_samples * calibrationInterp, numTx * numRx, n_chirps), dtype=complex)
    window = hann_local(n_samples)
    window2D = np.zeros((n_samples, n_chirps), dtype=float)
    for i in range(n_chirps):
        window2D[:, i] = window
    for i in range(numTx * numRx):
        data_to_Doppler[: , i , :] = rangeFFT_in[:, i, :]       
        Doppler_out[:, i , :] = np.fft.fftshift(np.fft.fft(data_to_Doppler[: , i , :], axis=1), axes = 1)
    return Doppler_out
