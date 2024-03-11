 #*- coding: UTF-8 -*-
 #encoding:utf-8

#base
n_samples = 512
n_chirps = 32
numTx = 13
numRx = 16

#calibration
calibrationInterp = 1
rangeResolution = 0.733
corner_reflecter_range = 5
search_index = 5

#cfar
RangeGuardUnitLen = 3
RefUnitLen = 5
SpeedGuardUnitLen = 2
cfra_threshold = 15

#doa
DoaTargetNum = 2
DoaTgtthr = 3
DBF_azi_range = 60
DBF_ele_range = 20
DBF_azi_step = 1200
DBF_ele_step = 200

#speed
lamda = 3.9e-3
chirpPeriod = 3.6e-5
chirpPeriodDelay = 1.2e-5
Vmax = lamda / 2 / (13 * chirpPeriod + chirpPeriodDelay)
Vres = Vmax / n_chirps