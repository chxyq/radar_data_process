 #*- coding: UTF-8 -*-
 #encoding:utf-8
import numpy as np
from numpy import fft
import parameter
import data_parse
import cmath
import matplotlib.pyplot as plt

n_samples = parameter.n_samples
n_chirps = parameter.n_chirps
numTx = parameter.numTx
numRx = parameter.numRx
calibrationInterp = parameter.calibrationInterp
rangeResolution = parameter.rangeResolution
corner_reflecter_range = parameter.corner_reflecter_range
search_index = parameter.search_index

def hann_local(len):
    win = (np.arange(1 , len + 1)) / float(len + 1)
    win = 0.5 - 0.5 * np.cos(2 * np.pi * win)
    return win

def genCalibrationMatric(path):
    full_data, datatype = data_parse.read_16bit_data(path)
    data_to_range = np.zeros((n_samples, numTx * numRx), dtype=complex)
    rangeFFT_out = np.zeros((n_samples * calibrationInterp, numTx * numRx), dtype=complex)
    rangeMat = np.zeros((numTx * numRx), dtype=np.int16)
    peakValMat = np.zeros((numTx * numRx), dtype=complex)  
    window = hann_local(n_samples)
    for i in range(numTx * numRx):
        data_to_range[: , i] = full_data[:, i, :].mean(axis = 1)
        start_point = int(corner_reflecter_range / rangeResolution - search_index) * calibrationInterp + 10
        end_point = int(corner_reflecter_range / rangeResolution + search_index) * calibrationInterp
        rangeFFT_out[:, i] = np.fft.fft(data_to_range[ : , i] * window, n_samples * calibrationInterp, axis=0)
        rangeMat[i] = np.argmax(abs(rangeFFT_out[start_point : end_point, i])) + start_point
        peakValMat[i] = rangeFFT_out[rangeMat[i], i]
    peakV = peakValMat[0] / peakValMat
    rangeM = (rangeMat[0] - rangeMat) * 2. * np.pi / n_samples / float(calibrationInterp)
    np.savetxt("./freq.txt", rangeM)
    np.savetxt("./amp.txt", abs(peakV))
    phase = np.zeros((numTx * numRx), dtype=float)
    for i in range(numTx * numRx):
        phase[i] = cmath.phase(peakV[i])
    np.savetxt("./phase.txt", phase)
    return rangeMat, peakValMat, rangeFFT_out

def Calculate_FreqCalMatric(rangeMat, idx):
    global f_calibration
    f_calibration = np.zeros((n_samples, numTx * numRx), dtype=float)
    n = np.arange(0, n_samples)
    for i in range(numTx * numRx):
        f_calibration[: , i] = (rangeMat[i] - rangeMat[idx]) * 2 * np.pi  * (1) / n_samples / calibrationInterp * n
    #return f_calibration

def Calculate_Phase_Angle_Calmatric(full_data, f_calibration, idx):
    global phase_angle_calibration
    phase_angle_calibration = np.zeros((n_samples, numTx * numRx), dtype=complex)
    data_to_range = np.zeros((n_samples, numTx * numRx), dtype=complex)
    rangeFFT_out = np.zeros((n_samples * calibrationInterp, numTx * numRx), dtype=complex)
    peakValMat = np.zeros((numTx * numRx), dtype=complex)
    window = hann_local(n_samples)
    rms = np.sqrt(np.sum(window ** 2) / n_samples)
    window /= rms  
    rangeMat = np.zeros((numTx * numRx), dtype=np.int16)
    for i in range(n_chirps):
        full_data[:, :, i] *= np.exp(-1j *f_calibration)
    for i in range(numTx * numRx):
        data_to_range[: , i] = full_data[:, i, :].mean(axis = 1)        
        rangeFFT_out[:, i] = np.fft.fft(data_to_range[: , i] * window, n_samples * calibrationInterp, axis=0)
        rangeMat[i] = np.argmax(abs(rangeFFT_out[:100 , i])) + 0
        peakValMat[i] = rangeFFT_out[rangeMat[i], i] 
    for i in range(n_samples):
        phase_angle_calibration[i, :] = peakValMat[idx] / peakValMat
    #return phase_angle_calibration  

def apply_calibration(full_data, f_calibration, phase_angle_calibration):
    for i in range(n_chirps):
        full_data[:, :, i] *= np.exp(-1j *f_calibration) * phase_angle_calibration
    return full_data

#rangeMat, peakValMat, rangeFFT_out = genCalibrationMatric("/home/tusimple/Downloads/5723_frame.bin")
