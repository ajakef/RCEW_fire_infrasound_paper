#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 16:06:45 2023

@author: jake
"""

def get_semblance_ticks(station):
    return {'TOP':[0.03, 0.1, 0.3, 1], 
            'QST': [0.1, 0.3, 1], 
            'VPT':[0.1, 0.3, 1]}.get(station, [0.25, 0.5, 1])

def get_t1_t2(station):
    return {'CHKB':['15:00', '23:59'], 
            'VPT':['18:00','22:50']}.get(station, ['12:00', '23:59'])

def get_power_ticks(station):
    return {'VPT': [3, 4, 5, 6, 7,8]}.get(station, [-2,-1,0,1,2,3])