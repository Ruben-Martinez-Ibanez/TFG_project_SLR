# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 19:12:44 2024

@author: ruben
"""

import math
import numpy as np
import pandas as pd
#from RadarLaserEqtnZhang.py import RadarLaserEqtnZhang
import matplotlib.pyplot as plt
import os
from Extract_Data_Functions import extract_data_CRD
import seaborn as sns;sns.set()
sns.set()

def RadarLaserEqtnZhang (lambda_laser, eta_q, E_t, A_r, rho, S, T, K_t, K_r, alfa, theta_t, R):

    #INPUT
    #Lambda_laser: longitud de onda del laser            (metros)
    #eta_q:        eficiencia cuántica del receptor      (adimensional)
    #E_t:          energía por pulso del laser           (julios)
    #eta_t:        eficiencia óptica de transmisión      (adimensional)
    #A_r:          área de la óptica de recepción        (metros2) 
    #rho:          reflectividad del satélite            (adimensional)
    #r_sat:        radio del satélite                    (metros)
    #theta_t       divergencia                           (radianes)
    #R:            distancia al satélite                 (metros)
    
    
    #eta_r:        eficiencia óptica de recepción        (adimensional)
    #Te:           transmisión atmosférica de ida        (adimensional)
    #Tc:           atenuación por cirros ida y vuelta    (adimensional)
    #R:            distancia al satélite                 (metros)
    #theta_t       divergencia                           (radianes)
    
    #OUTPUT
    #P_D:          probabilidad de detección
    #ns:           numero medio de fotones recibido por el fotodetector
    
    h            = 6.6260693 * 1E-34
    c            = 299792458
    
    numerador1 = lambda_laser * eta_q
    numerador2 = E_t * A_r * rho * S
    numerador3 = (T**2) * K_t * K_r * alfa
    denominador1 = h * c
    denominador2 = math.pi * (theta_t**2) * (R**4)

    n0 = numerador1 * numerador2 * numerador3 / denominador1 / denominador2 
    P = 1 - (math.exp(-n0))
    
    return(P, n0)

c = 299792458 # m/s

CRD_files = ['envisat_crd_20180131.frd', 'ajisai_20231208.fr2']
station_code: int = 7839

rho = [0.78, 1]
S   = [19.4975, math.pi*(2.15/2)**2] #m^2

for k in range(len(CRD_files)):
    day_t_CRD, fl_t_CRD, P0, T0, RH, lambda_laser, day_t_cal, obj_d_cal, del_t_cal, \
        MJD, YEAR, MONTH, DAY, satelitte = extract_data_CRD(CRD_files[k], station_code)
    range_CRD = []
    
    for i in range(len(day_t_CRD)):
        
        # Rango al objeto
        range_CRD.append(fl_t_CRD[i]*c/2)
    
    P_D_array = []
    n0_array = []
    
    for i in range(len(range_CRD)):       
        lambda_laser = 532 * 1E-9 #m
        eta_q        = 0.2
        E_t          = 0.06 #J
        A_r          = 0.5 #m^2
        T            = 0.6 
        K_t          = 0.6
        K_r          = 0.6
        alfa         = 10**(-13/10)
        theta_t      = math.radians(5/3600)
        R            = range_CRD[i] #800*1000
        # R            = np.linspace(800, 5000, num=100)*1000
        
        
        (P_D, n0) = RadarLaserEqtnZhang(lambda_laser, eta_q, E_t, A_r, rho[k], S[k], T, K_t, K_r, alfa, theta_t, R)
        P_D_array.append(P_D)
        n0_array.append(n0)
        
    print(f'{satelitte} Pr = {max(P_D_array):.3f}')

    plt.figure(figsize=(10,6)) 
    plt.plot(range_CRD, P_D_array, 'o')
    plt.xlabel('Distancia [m]', fontsize=20)
    plt.ylabel('Probabilidad', fontsize=20)
    #plt.title(f'P', fontsize=17)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)

plt.show()
plt.savefig('Probabilidad ecos láser', dpi=100.0)
