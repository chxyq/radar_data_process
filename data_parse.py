 #*- coding: UTF-8 -*-
 #encoding:utf-8
import numpy as np
import sys

def header_tailer_parse(data):
    pkt_len = int(hex(int(data[0]) & 0xffffffff), base = 16) / 2 ** 16
    disc = int(hex(int(data[0]) & 0x3ff), base = 16) / 2 ** 9
    pad = int(hex(int(data[0]) & 0x1ff), base = 16) / 2 ** 8
    eof = int(hex(int(data[0]) & 0x7f), base = 16) / 2 ** 6
    sof = int(hex(int(data[0]) & 0x1f), base = 16) / 2 ** 4
    datatype = int(hex(int(data[0]) & 0xf), base = 16)
    res_N = int(hex(int(data[1]) & 0xfff), base = 16)
    res_M = int(hex(int(data[1]) & 0xffffff), base = 16) / 2 ** 12
    res_K = int(hex(int(data[1]) & 0xffffffff), base = 16) / 2 ** 24
    frame_Id = data[2] + data[3] * 2 ** 32
    frame_timestamp = data[4] + data[5] * 10 ** 9
    tx_sys_time = data[6] + data[7] * 10 ** 9
    res_V = int(hex(int(data[8]) & 0xffffffff), base = 16) / 2 ** 16
    res_H = int(hex(int(data[8]) & 0xffff), base = 16)
    rf_frame_in_sys_time = data[10] + data[11] * 10 ** 9
    return datatype, int(res_N), int(res_M), int(res_K)

def read_16bit_data(filepath):
    data_i = np.fromfile(filepath, dtype='<i2')
    lens = int(len(data_i) / 4098)
    datatype = int(hex(int(data_i[0]) & 0xf), base = 16)

    #header
    head_data = data_i[: lens]
    datatype, n_samples, n_chirps, num_virtual_channel = header_tailer_parse(np.frombuffer(head_data, dtype = 'i4'))

    #tailer
    #tailer_data = data_i[len(data_i) - lens:]
    #datatype, n_samples, n_chirps, num_virtual_channel = header_tailer_parse(np.frombuffer(tailer_data, dtype = 'i4'))

    if lens == 1552 and n_chirps == 128:
        chirp_divided = 4
    elif lens == 1552 * 2 and n_chirps == 128:
        chirp_divided = 2
    else:
        chirp_divided = 1

    data16_i = data_i[lens: (int(lens / 2) + int(n_chirps / chirp_divided )* n_samples * (num_virtual_channel + 2)) * 2] 
    adcdata16_i = data16_i[::2] + 1j * data16_i[1::2]
    full_data_tmp16_i = adcdata16_i.reshape(num_virtual_channel + 2, n_samples, int(n_chirps / chirp_divided), order='F')
    full_data = np.zeros((num_virtual_channel, n_samples, n_chirps), dtype=complex)
    if datatype == 1 or datatype == 4:       
        for i in range(n_chirps):
            for j in range(n_samples):
                tmp_193 = full_data_tmp16_i[0, j, i].imag * 2 ** 16 + full_data_tmp16_i[0, j, i].real
                line_num = int(hex(int(tmp_193) & 0xffffff), base = 16)
                full_data[:, int(line_num % n_samples), int(line_num / n_samples)] = full_data_tmp16_i[2:, j, i]
    elif datatype == 2:
        data32_i = np.zeros(len(data16_i), dtype='i4')
        data32_i = data16_i * 2 ** 16
        data32_f = np.frombuffer(data32_i, dtype = "<f4")
        adcdata32_f = data32_f[::2] + 1j * data32_f[1::2]
        full_data_tmp32_f = adcdata32_f.reshape(num_virtual_channel + 2, n_samples, int(n_chirps / chirp_divided), order='F')
        for i in range(n_chirps):
            for j in range(n_samples):
                first_line = full_data_tmp16_i[0, j, int(i / chirp_divided)].imag * 2 ** 16 + full_data_tmp16_i[0, j, int(i / chirp_divided)].real
                if first_line > 0 :
                    line_vld = 0
                    continue
                else:
                    line_vld = 1
                line_num = int(hex(int(first_line) & 0xffffff), base = 16)
                full_data[:, int(line_num % n_samples), int(line_num / n_samples)] = full_data_tmp32_f[2:, j, int(i / chirp_divided)]
    else:
        print("unknow datatype:" + str(datatype))
        sys.exit()
    return np.transpose(full_data, (1, 0, 2)), datatype

def read_Orin_Asample(filepath):
    data_i = np.fromfile(filepath, dtype='<i4')
    lens = int(len(data_i) / 8194)
    #header
    head_data = data_i[: lens]
    datatype, n_samples, n_chirps, num_virtual_channel = header_tailer_parse(head_data)

    #tailer
    #tailer_data = data_i[len(data_i) - lens:]
    #datatype, n_samples, n_chirps, num_virtual_channel = header_tailer_parse(tailer_data)

    if lens == 772 and n_chirps == 128:
        chirp_divided = 4
    elif lens == 1544 and n_chirps == 128:
        chirp_divided = 2
    else:
        chirp_divided = 1

    data1_i = data_i[lens: (int(lens / 2) + int(n_chirps / chirp_divided )* n_samples * (num_virtual_channel + 1)) * 2]
    adcdata_i = data1_i[::2] + 1j * data1_i[1::2]
    full_data_tmp_i = adcdata_i.reshape(num_virtual_channel + 1, n_samples, int(n_chirps / chirp_divided ), order='F')
    full_data = np.zeros((num_virtual_channel, n_samples, n_chirps), dtype=complex)
    if datatype == 1 or datatype == 4:       
        for i in range(n_chirps):
            for j in range(n_samples):
                first_line = full_data_tmp_i[0, j, i].real
                data_len = int(hex(int(full_data_tmp_i[0, j, i].imag) & 0xffff), base = 16)
                line_num = int(hex(int(first_line) & 0xffffff), base = 16)
                full_data[:, int(line_num % n_samples), int(line_num / n_samples)] = full_data_tmp_i[1:, j, i]
    elif datatype == 2:
        data_f = np.fromfile(filepath, dtype='<f4')
        data1_f = data_f[lens: (int(lens / 2) + int(n_chirps / chirp_divided ) * n_samples * (num_virtual_channel + 1)) * 2]        
        adcdata_f = data1_f[::2] + 1j * data1_f[1::2]
        full_data_tmp_f = adcdata_f.reshape(num_virtual_channel + 1, n_samples, int(n_chirps / chirp_divided ), order='F')
        for i in range(n_chirps):
            for j in range(n_samples):
                first_line = full_data_tmp_i[0, j, int(i / chirp_divided )].real
                data_len = int(hex(int(full_data_tmp_i[0, j, int(i / chirp_divided )].imag) & 0xffff), base = 16)
                if first_line > 0 :
                    line_vld = 0
                    continue
                else:
                    line_vld = 1
                line_num = int(hex(int(first_line) & 0xffffff), base = 16)
                full_data[:, int(line_num % n_samples), int(line_num / n_samples)] = full_data_tmp_f[1:, j, int(i / chirp_divided )]
    else:
        print("unknow datatype:" + str(datatype))
        sys.exit()
    return np.transpose(full_data, (1, 0, 2)), datatype
