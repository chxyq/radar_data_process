 #*- coding: UTF-8 -*-
 #encoding:utf-8
import numpy as np
import parameter


numTx = parameter.numTx - 1
numRx = parameter.numRx
Tc = parameter.chirpPeriod
lamda = parameter.lamda


def doppler_compensation(peak_val, speed):
    peak_val_after_compensate = np.zeros((numTx * numRx), dtype=complex)
    for i in range(numTx):
        peak_val_after_compensate[i * 16: (i + 1) * 16] = peak_val[i * 16: (i + 1) * 16] * np.exp(-1j * speed * np.pi * 4 * Tc / lamda * i)
    return peak_val_after_compensate
