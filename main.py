 #*- coding: UTF-8 -*-
 #encoding:utf-8
import numpy as np
import matplotlib.pyplot as plt
from numpy import fft
from mpl_toolkits.mplot3d import Axes3D
import os
import scipy as sp
import io
import sys
from numpy.matlib import repmat
import matplotlib as mpl
import pandas as pd
import cmath
import struct
import seaborn as sns
import time
import detection
import data_parse
import doa
import doppler_process
import range_process
import parameter
import doppler_compensation
import calibration
import chirp_delay

n_samples = parameter.n_samples
n_chirps = parameter.n_chirps
numTx = parameter.numTx
numRx = parameter.numRx
calibrationInterp = parameter.calibrationInterp
RangeGuardUnitLen = parameter.RangeGuardUnitLen
RefUnitLen = parameter.RefUnitLen
SpeedGuardUnitLen = parameter.SpeedGuardUnitLen
cfra_threshold = parameter.cfra_threshold
rangeResolution = parameter.rangeResolution / calibrationInterp
DBF_azi_range = parameter.DBF_azi_range
DBF_ele_range = parameter.DBF_ele_range
DBF_azi_step = parameter.DBF_azi_step
DBF_ele_step = parameter.DBF_ele_step

def hann_local(len):
    win = (np.arange(1 , len + 1)) / float(len + 1)
    win = 0.5 - 0.5 * np.cos(2 * np.pi * win)
    return win

def draw_rdm_and_cfar_curve(cfar_curve, target, dopploer_out):
    x = np.linspace(0, n_samples * calibrationInterp - 1, n_samples * calibrationInterp) * rangeResolution
    rdm = abs(dopploer_out).sum(axis = 1)
    y = 20 * np.log10(rdm)
    plt.subplot(2, 2, 1)
    for i in target:
        #if i[1] == n_chirps / 2:
            #point += str(np.around(i[0] * 0.707, 2)) + "m;"
        plt.scatter(i[0] * rangeResolution , y[i[0], i[1]], s = 50)
    plt.plot(x, y)
    plt.plot(x, cfar_curve[:, int(n_chirps / 2)] + cfra_threshold)
    plt.subplot(2, 2, 3)
    plot = sns.heatmap(pd.DataFrame(y, index = np.int16(x)))
    for i in target:
        plt.scatter(i[1], i[0] , color = 'y' , s = 10)

def draw_point_cloud_result(point_cloud):
    col = ['range', 'speed', 'azi', 'ele', 'power']
    row = []
    for i in range(len(point_cloud)):
        tmp = "Target" + str(i)
        row.append(tmp)
    plt.subplot(1, 2, 2)
    plt.axis('off')
    #plt.set_fontsize(14)
    #plt.table().auto_set_font_sizae(False)
    plt.table(cellText = point_cloud[: 27], rowLabels = row[ : 27], colLabels = col, cellLoc = 'center', loc = 'center')
    #plt.show()

def main(path):
    full_data, datatype = data_parse.read_16bit_data(path)
    #full_data, datatype = data_parse.read_Orin_Asample(path)
    if datatype == 1:
        f_calibration = np.loadtxt("/home/tusimple/workspace/code/cn-aeg-radar-test/data_parse_for_algo/008_f_cal.txt", dtype=float)
        phase_angle_calibration = np.loadtxt("/home/tusimple/workspace/code/cn-aeg-radar-test/data_parse_for_algo/008_p_a_cal.txt", dtype=complex)
        full_data_cal = calibration.apply_calibration(full_data, f_calibration, phase_angle_calibration)
        rangeFFT_out = range_process.range_data_process(full_data_cal)
    elif datatype == 2:
        rangeFFT_out = full_data   
    elif datatype == 4:
        rangeFFT_out = range_process.range_data_process(full_data)
    print(datatype)
    '''x = np.linspace(0, 511, 512)
    plt.plot(x, abs(rangeFFT_out.mean(axis = 2)))
    plt.show()'''
    doppler_out = doppler_process.Doppler_FFT(rangeFFT_out)
    cfar_curve, target = detection.Cfar(doppler_out)
    point = ""
    point_cloud = []
    for obj in target:
        angle_str = ''
        speed = chirp_delay.calculate_speed(doppler_out[obj[0], 0, obj[1]], doppler_out[obj[0], 16, obj[1]], obj[1])
        data, angle_result, power = doa.draw_2DDBF_angle(doppler_compensation.doppler_compensation(doppler_out[obj[0] ,  16: , obj[1]], speed))
        if len(angle_result) == 0:
            continue
        for j in range(len(angle_result)):
            azi_angle = np.around((angle_result[j][0][0] - DBF_azi_step / 2) * DBF_azi_range * 2. / (DBF_azi_step + 1), 2)
            ele_angle = np.around((angle_result[j][1][0] - DBF_ele_step / 2) * DBF_ele_range * 2. / (DBF_ele_step + 1), 2)
            pc = [np.around(obj[0] * rangeResolution, 2), speed, azi_angle, ele_angle, np.around(20 * np.log10(power), 2)]
            angle_str += 'azi:' + str(azi_angle) + ",ele:" + str(ele_angle) + ','
            point_cloud.append(pc)

        plot = sns.heatmap(data)
        plt.title(angle_str)
        #if abs(np.around(obj[0] * rangeResolution, 2) - 150) < 4:
        #plt.savefig(path + '-' + str(np.around(obj[0] * rangeResolution , 2)) + "m.png")
        plt.clf()
    plt.figure(figsize = (12, 6))
    #plt.subplot(2, 2, 4)
    #plot = sns.heatmap(save_data)
    #plt.title(save_angle)
    np.savetxt(path + "pc.txt", point_cloud)
    draw_point_cloud_result(point_cloud)
    draw_rdm_and_cfar_curve(cfar_curve, target, doppler_out)
    #plt.show()
    plt.savefig(path + ".png", dpi = 300)
    plt.clf()

main("/home/tusimple/workspace/20240129/d100_s20/d100_s20_az0_el5/214160_frame.bin")
