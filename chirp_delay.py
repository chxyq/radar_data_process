 #*- coding: UTF-8 -*-
 #encoding:utf-8
 
import numpy as np
import cmath
import parameter


n_chirps = parameter.n_chirps
chirpPeriod = parameter.chirpPeriod
chirpPeriodDelay = parameter.chirpPeriodDelay
Vres = parameter.Vres
Vmax = parameter.Vmax

def calculate_speed(delay_tx, ori_tx, speed_index):
    phase1 = cmath.phase(delay_tx)
    phase2 = cmath.phase(ori_tx)
    frd = np.pi * 2 * (speed_index - n_chirps / 2) / n_chirps
    delta_phase = (phase2 - phase1) - frd * (chirpPeriod + chirpPeriodDelay) / (chirpPeriod * 13 + chirpPeriodDelay)
    while delta_phase > np.pi:
        delta_phase -= 2 * np.pi
    while delta_phase < -np.pi:
        delta_phase += 2 * np.pi
    temp_min = float('inf')
    for q in range(-4, 6):
        datatemp = delta_phase / 2 / np.pi - q * (chirpPeriod + chirpPeriodDelay) / (chirpPeriod * 13 + chirpPeriodDelay)
        while datatemp > 0.5:
            datatemp -= 1
        while datatemp <= -0.5:
            datatemp += 1

        #print(datatemp)
        if abs(datatemp) > temp_min:
            continue
        else:
            temp_min = abs(datatemp)
            iq = q
    v = (speed_index - n_chirps / 2) * Vres + iq * Vmax

    return v


