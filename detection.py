 #*- coding: UTF-8 -*-
 #encoding:utf-8
 
import numpy as np
import parameter

calibrationInterp = parameter.calibrationInterp
n_samples = parameter.n_samples * calibrationInterp
n_chirps = parameter.n_chirps
RangeGuardUnitLen = parameter.RangeGuardUnitLen
RefUnitLen = parameter.RefUnitLen
SpeedGuardUnitLen = parameter.SpeedGuardUnitLen
cfra_threshold = parameter.cfra_threshold

def Cfar(Doppler_out):
    target = []
    cfar_out = 20 * np.log10(abs(Doppler_out).sum(axis = 1))
    noise_gate = np.zeros((n_samples, n_chirps), dtype=float)
    for i in range(n_samples):
        for j in range(n_chirps):
            num_of_point = 4. * RefUnitLen
            if i <= RangeGuardUnitLen:
                rangeRefer1 = np.zeros(1, dtype=float)
                num_of_point -= RefUnitLen
                rangeRefer2 = cfar_out[i + RangeGuardUnitLen + 1: i + RangeGuardUnitLen + RefUnitLen + 1, j]
            elif i <= RangeGuardUnitLen + RefUnitLen:
                rangeRefer1 = cfar_out[: i - RangeGuardUnitLen, j]
                num_of_point -= (RefUnitLen - len(rangeRefer1))
                rangeRefer2 = cfar_out[i + RangeGuardUnitLen + 1: i + RangeGuardUnitLen + RefUnitLen + 1, j]
            elif i >= n_samples - RangeGuardUnitLen - 1:
                rangeRefer1 = cfar_out[i - RangeGuardUnitLen - RefUnitLen : i - RangeGuardUnitLen, j]
                rangeRefer2 = np.zeros(1, dtype=float)
                num_of_point -= RefUnitLen
            elif i >= n_samples - RangeGuardUnitLen - RefUnitLen - 1:
                rangeRefer1 = cfar_out[i - RangeGuardUnitLen - RefUnitLen : i - RangeGuardUnitLen, j]
                rangeRefer2 = cfar_out[i + RangeGuardUnitLen + 1:, j]
                num_of_point -= (RefUnitLen - len(rangeRefer2))
            else:
                rangeRefer1 = cfar_out[i - RangeGuardUnitLen - RefUnitLen : i - RangeGuardUnitLen, j]
                rangeRefer2 = cfar_out[i + RangeGuardUnitLen + 1: i + RangeGuardUnitLen + RefUnitLen + 1, j]
            
            if j <= SpeedGuardUnitLen:
                speedRefer1 = np.zeros(1, dtype=float)
                num_of_point -= RefUnitLen
                speedRefer2 = cfar_out[i, j + SpeedGuardUnitLen + 1 : j + SpeedGuardUnitLen + RefUnitLen + 1]
            elif j <= SpeedGuardUnitLen + RefUnitLen:
                speedRefer1 = cfar_out[i, : j - SpeedGuardUnitLen]
                num_of_point -= (RefUnitLen - len(speedRefer1))
                speedRefer2 = cfar_out[i, j + SpeedGuardUnitLen + 1 : j + SpeedGuardUnitLen + RefUnitLen + 1]
            elif j >= n_chirps - SpeedGuardUnitLen - 1:
                speedRefer1 = cfar_out[i, j - SpeedGuardUnitLen - RefUnitLen : j - SpeedGuardUnitLen]
                speedRefer2 = np.zeros(1, dtype=float)
                num_of_point -= RefUnitLen
            elif j >= n_chirps - SpeedGuardUnitLen - RefUnitLen - 1:
                speedRefer1 = cfar_out[i, j - SpeedGuardUnitLen - RefUnitLen : j - SpeedGuardUnitLen]
                speedRefer2 = cfar_out[i, j + SpeedGuardUnitLen + 1:]
                num_of_point -= (RefUnitLen - len(speedRefer2))
            else:
                speedRefer1 = cfar_out[i, j - SpeedGuardUnitLen - RefUnitLen : j - SpeedGuardUnitLen]
                speedRefer2 = cfar_out[i, j + SpeedGuardUnitLen + 1 : j + SpeedGuardUnitLen + RefUnitLen + 1]
            noise_gate[i, j] = (rangeRefer1.sum() + rangeRefer2.sum() + speedRefer1.sum() + speedRefer2.sum()) / num_of_point
            current_val = cfar_out[i, j]
            if i > 0 and i < n_samples - 1:
                if j == 0:
                    if cfar_out[i - 1, j] < current_val and current_val > cfar_out[i + 1, j] and current_val > cfar_out[i, j + 1] and current_val > noise_gate[i, j] + cfra_threshold:
                        target.append([i, j])
                elif j == n_chirps - 1:
                    if cfar_out[i - 1, j] < current_val and current_val > cfar_out[i + 1, j] and current_val > cfar_out[i, j - 1] and current_val > noise_gate[i, j] + cfra_threshold:
                        target.append([i, j])
                else:
                    if cfar_out[i - 1, j] < current_val and current_val > cfar_out[i + 1, j] and current_val > cfar_out[i, j + 1] and current_val > cfar_out[i, j - 1] and current_val > noise_gate[i, j] + cfra_threshold:
                        target.append([i, j])
    return noise_gate, target
