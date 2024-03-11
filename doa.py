 #*- coding: UTF-8 -*-
 #encoding:utf-8
import numpy as np
import parameter
import pandas as pd

numTx = parameter.numTx - 1
numRx = parameter.numRx
DoaTargetNum = parameter.DoaTargetNum
DoaTgtthr = parameter.DoaTgtthr
DBF_azi_range = parameter.DBF_azi_range
DBF_ele_range = parameter.DBF_ele_range
DBF_azi_step = parameter.DBF_azi_step
DBF_ele_step = parameter.DBF_ele_step

#A samples
dd = 0.5
TX_XLabel = [24.8088 , 9.5639, 4.9772, 6.1222, 10.6648, 21.2900, 47.3147, 39.1872, 31.3306, 25.8167 ,40.0374 ,44.8059]
TX_YLabel = [7.0735, 0.0524, 0.1095, 44.3448, 44.0243, 43.4385, 0.0249, 0.0251, 0.064, 44.3523, 44.3567, 44.5298]
RX_XLabel = [0.0312, 0.0297, 1.9757, 5.4591, 11.1095, 0.0347, 0.8217, 0.0026, 38.544, 47.7889, 50.966, 50.9849, 50.9704, 50.3815, 50.9769, 29.4868]
RX_YLabel = [1.6866, 10.1348, 18.5543, 21.5546, 23.2009, 26.9323, 35.2838, 43.927, 21.1135, 17.658, 13.4477, 3.6501, 43.1119, 33.3586, 24.9838, 20.696]
#steer_vector = np.load("./A_steer_vector.bin.npy")

def draw_2DDBF_angle(peakValMat_cal):
    coordinary_list = []
    angle_result = []
    for i in range(numTx):
        for j in range(numRx):
            locx = TX_XLabel[i] + RX_XLabel[j]
            locy = TX_YLabel[i] + RX_YLabel[j]
            coordinary_list.append([i * numRx + j, locx, locy])
    
    azi = np.linspace(-DBF_azi_range, DBF_azi_range, DBF_azi_step + 1)
    ele = np.linspace(-DBF_ele_range, DBF_ele_range, DBF_ele_step + 1)
    #peakValMat_cal = np.ones(192, dtype=complex)

    ###Calculte Azi Angle###
    steer_vector = np.zeros((DBF_azi_step + 1, DBF_ele_step + 1, numTx * numRx), dtype=float)  
    for i in range(numTx * numRx):
        steer_vector[:, :, i] = np.sin((azi) / 180. * np.pi).reshape(DBF_azi_step + 1, 1) * np.cos(ele / 180. * np.pi).reshape(1, DBF_ele_step + 1) * coordinary_list[i][1] + np.sin(ele / 180. * np.pi).reshape(1, DBF_ele_step + 1) * coordinary_list[i][2]
    steer_vector = np.exp(-1j * np.pi * 2 * dd * steer_vector)
    #np.save("./A_steer_vector.bin", steer_vector)
    #steer_vector = np.load("./A_steer_vector.bin.npy")
    dbf_after = np.dot(steer_vector, peakValMat_cal)
    dbf_after_log = 20 * np.log10(abs(dbf_after) / np.max(abs(dbf_after)))
    data = pd.DataFrame((abs(dbf_after) / np.max(abs(dbf_after))).T)
    angle_result = []
    peak_set = np.sort(np.extract(dbf_after_log >= -float(DoaTgtthr), dbf_after_log))[::-1]
    angle_result = []
    num_over_thre = len(peak_set)
    for i in range(num_over_thre):
        if num_over_thre > 200:
            pass
            #break
        loc = np.where(dbf_after_log == peak_set[i])
        if len(loc[0]) >= 2:
            continue
        if loc[0] <= 5 or loc[0] >= DBF_azi_step - 6 or loc[1] <= 5 or loc[1] >=  DBF_ele_step - 6:
            if len(angle_result) == 0:
                pass
            else:
                continue
        if len(angle_result) == 0:         
            angle_result.append(loc)
            continue
        elif len(angle_result) == DoaTargetNum:
            break
        else:
            left = dbf_after_log[loc[0] - 1, loc[1]]
            right = dbf_after_log[loc[0] + 1, loc[1]]
            up = dbf_after_log[loc[0], loc[1] - 1]
            down = dbf_after_log[loc[0], loc[1] + 1]                       
            if peak_set[i] > left and peak_set[i] > right and peak_set[i] > up and peak_set[i] > down:
                angle_result.append(loc)
    return data, angle_result, np.max(abs(dbf_after))

    
