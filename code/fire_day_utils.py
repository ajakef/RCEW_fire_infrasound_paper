#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 16:06:45 2023

@author: jake
"""
import numpy as np
import cleanbf
import obspy

def get_semblance_ticks(station):
    return {'TOP':[0.03, 0.1, 0.3, 1], 
            'QST': [0.1, 0.3, 1], 
            'VPT':[0.1, 0.3, 1]}.get(station, [0.25, 0.5, 1])

def get_t1_t2(station):
    return {'CHKB':['15:00', '23:59'], 
            'VPT':['18:00','22:50']}.get(station, ['12:00', '23:59'])

def get_power_ticks(station):
    return {'VPT': [3, 4, 5, 6, 7,8]}.get(station, [-2,-1,0,1,2,3])


def beam_stack_spectrum(st, sx, sy, win_length_sec = 10):
    print(f'({sx:0.4f}, {sy:0.4f})')
    cross_spec, FT, freqs, dfn, dfd = cleanbf.calc_cross_spectrum(st, win_length_sec = win_length_sec)
    weight = cleanbf.calc_weights(cleanbf.make_steering_vectors(st, freqs, [sx], [sy]))[:,0,0,:]
    P = np.abs(np.einsum('ki,ijk,kj->k', weight.conj(), cross_spec, weight)) # i, j: stations; k: freq.
    semblance = P/np.einsum('iik->k', cross_spec)
    return {'freqs':freqs, 'power':P, 'semblance':semblance, 'N':len(st)}

def beam_stack_spectrogram(st, forward_az, slowness, win_length_sec = 40, welch_ratio = 4, overlap = 0):
    num_windows = calc_num_windows(st[0].stats.endtime - st[0].stats.starttime, win_length_sec, overlap)
    sx = np.sin(forward_az * np.pi/180) * slowness
    sy = np.cos(forward_az * np.pi/180) * slowness
    if (len(sx) > 1) and (len(sx) != num_windows):
        raise Exception(f'sx and sy are the wrong length (should be 1 or {num_windows})')
                                
    f_inputs = []
    for sxx, syy in zip(sx, sy):
        f_inputs.append([sxx, syy, win_length_sec / welch_ratio])
    return apply_function_windows(st, beam_stack_spectrum, win_length_sec, overlap = overlap, f_inputs = f_inputs)

def apply_function_windows(st, f, win_length_sec, overlap = 0.5, f_inputs = [], verbose = True):
    """
    Run an analysis (or suite of analyses) on overlapping windows for some dataset
    
    Parameters:
    -----------
    st : obspy.Stream
    Stream including data to divide into windows and analyze

    f : function
    Accepts single variable "st" (obspy.Stream), returns dictionary of results

    win_length_sec : float
    Length of analysis windows [seconds]

    overlap : float
    Proportion of a window that overlaps with the previous window [unitless, 0-1]
    f_inputs: list
    List of inputs to provide to f() in addition to st (1 for each time window)

    Returns:
    --------
    dictionary with following items:
    --t_mid (obspy.UTCDateTime): mean time of each window
    --all elements of the output of "f", joined into numpy arrays

    Note:
    -----
    For each time window, the stream is trimmed to fit, but not tapered, detrended, or otherwise 
    processed. If those steps are necessary, be sure they are carried out in f().
"""
    # f must input a stream and return a dict
    eps = 1e-6
    t1 = st[0].stats.starttime
    t2 = st[0].stats.endtime
    data_length_sec = t2 - t1
    #num_windows = 1 + int(np.ceil((data_length_sec - win_length_sec) / (win_length_sec * (1 - overlap)) - eps))
    num_windows = calc_num_windows(data_length_sec, win_length_sec, overlap)
    if (len(f_inputs) > 0) and (len(f_inputs) != num_windows):
        raise Exception(f'len(f_inputs) is not the same as num_windows {num_windows}')
    print(num_windows)
    for i in range(num_windows):
        if verbose and ((i % (num_windows//100)) == 0):
            print(f'{i} of {num_windows}')
        win_start = t1 + i*(data_length_sec - win_length_sec) / (num_windows-1)
        st_tmp = st.slice(win_start-eps, win_start + win_length_sec - eps, nearest_sample = False)
        if len(f_inputs) > 0:
            win_dict = f(st_tmp, *(f_inputs[i]))
        else:
            win_dict = f(st_tmp)
        if i == 0:
            output_dict = {key:[] for key in win_dict.keys()}
            output_dict['t_mid'] = []
        for key in win_dict.keys():
            output_dict[key].append(win_dict[key])
        output_dict['t_mid'].append(win_start + win_length_sec/2)
    output_dict = {key:np.array(output_dict[key]) for key in output_dict.keys()}
    return output_dict

def calc_num_windows(data_length_sec, win_length_sec, overlap):
    eps = 1e-3
    return 1 + int(np.floor((data_length_sec - win_length_sec) / (win_length_sec * (1 - overlap)) - eps))
