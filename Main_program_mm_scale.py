# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 12:46:30 2023

@author: ruben
"""

import numpy as np
import matplotlib.pyplot as plt
import math

from Extract_Data_Functions import extract_data_CRD
from Extract_Data_Functions import extract_data_CPF
from Extract_Data_Functions import extract_data_station
from CustomLagrange import lagrange
from transformations import ecef2lla
from transformations import geo2topo
from MariniMurrayCorrections import Marini_Murray

import seaborn as sns;sns.set()
sns.set()

#%% INPUTS

### ------ ENVISAT ------
## Datos Artículo ##
# 2013-07-12
# CRD_file: str = '7839_CRD_20130712_2040.frd'
# CPF_file: str = 'envisat_cpf_130711_6921.dlr'
# data_section = [74690, 74720]
# station_code: int = 7839

# 2013-07-14
# CRD_file: str = '7839_CRD_20130714_2106.frd'
# CPF_file: str = 'envisat_cpf_130711_6921.dlr'
# data_section = [200, 250]
# station_code: int = 7839

#######
 
# 2018-01-31 ### 7839
# CRD_file: str = 'envisat_crd_20180131.frd'
# CPF_file: str = 'envisat_cpf_180125_5251.dlr'
# data_section = [270, 370]
# station_code: int = 7839

# 2018-01-29 ### 7839
# CRD_file: str = 'envisat_crd_20180129.frd'
# CPF_file: str = 'envisat_cpf_180125_5251.dlr'
# data_section = [300, 350]
# station_code: int = 7839

# 2019-04-27 ### 7839
# CRD_file: str = 'envisat_20190930.frd'
# CPF_file: str = 'envisat_cpf_190929_7731.aas'
# data_section = [10, 30]
# station_code: int = 7839

# 2019-10-23 ### 7839 FUNCIONA!!!
# CRD_file: str = 'envisat_20191023.frd'
# CPF_file: str = 'envisat_cpf_191021_7941.aas'
# data_section = [190, 250]
# station_code: int = 7839

CRD_file: str = 'ajisai_20231208.fr2'
CPF_file: str = 'ajisai_cpf_231204_33801.dgf'
data_section = [300, 320]
station_code: int = 7839





stations_file: str = 'SLRF2014_POS+VEL_2030.0_180504.snx.txt'

'''
Código coordenadas Graz, Austria: 7839
Código coordenadas San Fernando, España: 7824
Código coordenadas Matera, Italia: 7941
Código coordenadas Borowiec, Polonia: 7811
'''

print('## 00 ## - Se han leído como Inputs los nombre de los ficheros:')
print(f'    - Archivo CRD: {CRD_file}')
print(f'    - Archivo CPF: {CPF_file}')
print(f'    - Archivo de las estaciones: {stations_file}')
print(f'    - Código de la estación: {station_code}')

#%% CRD & CPF data

day_t_CRD, fl_t_CRD, P0, T0, RH, lambda_laser, day_t_cal, obj_d_cal, del_t_cal, \
    MJD, YEAR, MONTH, DAY, satelitte = extract_data_CRD(CRD_file, station_code)

day_t_CPF, position_CPF = extract_data_CPF(CPF_file, MJD)

if MONTH < 10:
    if DAY < 10:
        fecha_datos = f'{YEAR}-0{MONTH}-0{DAY}'
    else:
        fecha_datos = f'{YEAR}-0{MONTH}-{DAY}'
else:
    if DAY < 10:
        fecha_datos = f'{YEAR}-{MONTH}-0{DAY}'
    else:
        fecha_datos = f'{YEAR}-{MONTH}-{DAY}'


print('## 01 ## - Se han leído los archivos y se han extraído los datos necesarios.')

#%% Coordenadas geocéntricas del satélite

c = 299792458 # m/s
day_t_interp = []

range_CRD = []

for i in range(len(day_t_CRD)):
    # Tiempo de interpolación
    day_t_interp.append(day_t_CRD[i] + fl_t_CRD[i]/2)
    
    # Rango al objeto
    range_CRD.append(fl_t_CRD[i]*c/2)

x_geo = position_CPF[0]
y_geo = position_CPF[1]
z_geo = position_CPF[2]

plt.figure(figsize=(10,7))
plt.plot(day_t_CRD, range_CRD, 'o')
plt.xlabel('Tiempo [s]',fontsize=20)
plt.ylabel('Rango [m]',fontsize=20)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.title(f'{satelitte} - {fecha_datos}',fontsize=17)
plt.savefig(f'0. Rango - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
plt.show()

print(f'Seguimiento de {day_t_CRD[-1]-day_t_CRD[0]:.0f} segundos')

print('## 02 ## - Se han calculado las coordenadas geocéntricas del satélite.')

#%% Coordenadas interpoladas del satélite (interpolador facilitado - FUNCIONA)

NMAX = len(day_t_CPF)
N = 9

x_obs = []
y_obs = []
z_obs = []

for i in range(len(day_t_CRD)):
    
    X_int = lagrange(day_t_CPF, x_geo, N, day_t_interp[i])
    Y_int = lagrange(day_t_CPF, y_geo, N, day_t_interp[i])
    Z_int = lagrange(day_t_CPF, z_geo, N, day_t_interp[i])
    
    x_obs.append(X_int)
    y_obs.append(Y_int)
    z_obs.append(Z_int)
    
print('## 03 ## - Se han determinado las coordenadas interpoladas del satélite \
con los datos de tiempo del CPF.')

#%% Coordenadas geocéntricas de la estación láser donde se toman las medidas

st_position, st_velocity, data_st_epoch = extract_data_station(str(station_code), stations_file)
x_st = st_position[0]
y_st = st_position[1]
z_st = st_position[2]

print('## 04 ## - Se han determinado las coordenadas de la estación (sin corrección):')
print(f'    x = {x_st:.3f} m')
print(f'    y = {y_st:.3f} m')
print(f'    z = {z_st:.3f} m')


#%% Distancia (rango) de la estación al satélite 
# r = sqrt[(x_st-x_obs)^2 + (y_st-y_obs)^2 + (z_st-z_obs)^2]

range_int = []
for i in range(len(x_obs)):
    range_int.append(math.sqrt((x_st - x_obs[i])**2 + \
                               (y_st - y_obs[i])**2 + \
                               (z_st - z_obs[i])**2))
print('## 05 ## - Se han calculado las distancias de la estación al satélite.')
      
#%% Coordenadas topocéntricas

latitud_st, longitud_st, altitud_st = ecef2lla(x_st, y_st, z_st)
coords_station =  x_st, y_st, z_st, latitud_st, longitud_st, altitud_st
azimuth, elevation = geo2topo(x_obs, y_obs, z_obs, coords_station, range_int)

print('## 06 ## - Se han calculado las coordenadas topocéntricas de la estación.')

#%% Correcciones (Extended ranging equation) - Marini Murray
'''
d = (1/2)cΔt + Δd_0 + Δd_s + Δd_b + Δd_r + η

(1/2)cΔt: = range_CRD
Δd_0: corrección de excentricidad en tierra = 0 (no se aplica) χ
Δd_s: corrección de excentricidad del satélite ?
Δd_b: retardo de la estación --> CRD (40) ?
Δd_r: corrección modelo Marini Murray = corr_range_MM
η: errores observacionales y sistemáticos restantes = (no se aplica) χ


P0 = 962.50 
T0 = 294.46
RH = 57.8
lambda_laser = 0.532 # micrometros
'''
corr_range_MM = []
for i in range(len(elevation)):
    corr_range_MM.append(Marini_Murray(latitud_st, altitud_st, elevation[i], \
                                       P0, T0, RH, lambda_laser))

print('## 07 ## - Se ha aplicado la corrección de Marini Murray.')
#%% Residuos

residuos = []
corr_residuos = []
min_day_t_part = math.inf
for i in range(len(range_CRD)):
    if day_t_CRD[i] < min_day_t_part:
        min_day_t_part = day_t_CRD[i]
    day_t_CRD[i] = day_t_CRD[i]-min_day_t_part
    residuos.append(range_CRD[i] - range_int[i])
    corr_residuos.append(residuos[i] - corr_range_MM[i])

plt.figure(figsize=(10,7))
plt.plot(day_t_CRD, corr_residuos, 'o')
plt.xlabel('Time [s]',fontsize=20)
plt.ylabel('Range residuals [m]',fontsize=20)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.ylim(-8, 7)
#plt.title(f'{satelitte} - {fecha_datos}',fontsize=17)
plt.savefig(f'1. Residuos sin tratar - {satelitte} - {fecha_datos}.jpg', dpi=150.0)
plt.show()

print('## 08 ## - Se han calcuado los residuos (plot 1).')


#%% OscilacionesEscalaMilimétrica
'''
x = day_t_CRD
y = residuos_detrend


xx = np.linspace(min(day_t_CRD), max(day_t_CRD), len(x))

yy = []
for i in range(len(xx)):
    yy.append(lagrange(x, y, 7, xx[i]))
    print(i)
plt.plot(xx, yy, 'o')

'''
intervalo_t = 20

subsection:int = int((max(data_section)-min(data_section))/intervalo_t)
day_t_part_full = []
residuos_mil_scale_full = []
min_day_t_part = math.inf

for i in range(subsection):
    data_subsection = [min(data_section) + i*intervalo_t, min(data_section) + (i+1)*intervalo_t]
    day_t_part = []
    residuos_part = []
    day_t_part_plot = []
    
    with open('Datos_Residuos_escala_métrica.txt', 'w') as r:
        r.write('Tiempo/s Residuos/m\n')
        for i in range(len(day_t_CRD)):
            if day_t_CRD[i] > data_subsection[0]:
                if day_t_CRD[i] < data_subsection[1]:
                    if day_t_CRD[i] < min_day_t_part:
                        min_day_t_part = day_t_CRD[i]
                    day_t_part.append(day_t_CRD[i]-min_day_t_part)
                    day_t_part_plot.append(day_t_CRD[i])
                    day_t_part_full.append(day_t_CRD[i]-min_day_t_part)
                    residuos_part.append(corr_residuos[i])
                    r.write(f'{day_t_CRD[i]} {corr_residuos[i]}\n')
    
    polinomio = np.polyfit(day_t_part, residuos_part, 5)
    residuos_pol = np.polyval(polinomio, day_t_part)
    
    residuos_mil_scale = []
    for i in range(len(day_t_part)):
        residuos_mil_scale.append((residuos_pol[i]-residuos_part[i])*1000)
        residuos_mil_scale_full.append((residuos_pol[i]-residuos_part[i])*1000)
    with open('Datos_Residuos_escala_milimétrica.txt', 'w') as r:
        r.write('Tiempo/s Residuos/mm\n')
        for i in range(len(day_t_part)):
            r.write(f'{day_t_part[i]} {residuos_mil_scale[i]}\n')
    
    #Plot Sección Escala Métrica
    plt.figure(figsize=(10,7))
    plt.plot(day_t_CRD, corr_residuos, 'o',  day_t_part_plot, residuos_pol, '.')
    plt.xlabel('Tiempo [s]',fontsize=20)
    plt.ylabel('Residuos [m]',fontsize=20)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.title(f'{satelitte} - {fecha_datos}',fontsize=17)
    plt.savefig(f'5. Sección de residuos en Escala métrica - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
    plt.show()
    
#residuos_mil_scale_detrend = signal.detrend(residuos_mil_scale)

# Plot Escala Milimétrica
plt.figure(figsize=(22*subsection,15))
plt.plot(day_t_part_full, residuos_mil_scale_full, 'o')
plt.xlabel('Time [s]',fontsize=45)
plt.ylabel('Range residuals [mm]',fontsize=45)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
plt.xlim(7, 9)
#plt.title(f'{satelitte} - {fecha_datos}',fontsize=37)
plt.savefig(f'6. Sección de residuos en Escala milimétrica - {satelitte} - {fecha_datos}.jpg', dpi=200.0)
plt.show()

#%% Ajuste Sinusoidal Escala milimétrica

y = []
x = []

for i in range(len(day_t_part_full)):
    x.append(day_t_part_full[i])
    y.append(10*np.sin(29.099*x[i]))

# Plot Escala milimétrica Ajuste sinusoidal
plt.figure(figsize=(22, 15))
plt.plot(day_t_part_full, residuos_mil_scale_full, 'o')
plt.plot(x, y, linewidth=8)
plt.xlabel('Time [s]',fontsize=45)
plt.ylabel('Range residuals [mm]',fontsize=45)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
plt.xlim(7, 9)
#plt.title(f'{satelitte} - {fecha_datos}',fontsize=37)
plt.savefig(f'Sección de residuos en Escala milimétrica sinusoidal - {satelitte} - {fecha_datos}.jpg', dpi=200.0)
plt.show()

#%%
res_peak = []
t_peak = []

part = [7.8, 8.08]

for i in range(len(day_t_part_full)):
    if day_t_part_full[i] > part[0] and day_t_part_full[i] < part[-1]:
        res_peak.append(residuos_mil_scale_full[i])
        t_peak.append(day_t_part_full[i])


tt = np.linspace(part[0], part[-1], 100)
p = np.polyfit(t_peak, res_peak, 6)
ptt = np.polyval(p, tt)

AC = np.inf
for i in range(len(tt)):
    if ptt[i] < AC:
        AC = ptt[i]
        t_AC = tt[i]

plt.figure(figsize=(10*subsection,15))
plt.plot(t_peak, res_peak, 'o')
plt.plot(tt, ptt, linewidth= 7.0)
plt.plot(t_AC, AC, 'ko', markersize= 20)
plt.text(t_AC+0.005, AC-0.7, 'AC', weight = 'heavy', size = 'xx-large')
plt.xlabel('Time [s]',fontsize=45)
plt.ylabel('Range residuals [mm]',fontsize=45)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
#plt.xlim(7.8, 8.05)
#plt.title(f'{satelitte} - {fecha_datos}',fontsize=37)
plt.savefig(f'Pico milimétrico y AC - {satelitte} - {fecha_datos}.jpg', dpi=200.0)
plt.show()

print(f'AC = {AC:.1f} mm.')


#%%

# plt.figure(figsize=(25,15))
# plt.plot(day_t_part_full, residuos_mil_scale_full, 'o')
# plt.xlabel('Tiempo [s]',fontsize=45)
# plt.ylabel('Residuos [mm]',fontsize=45)
# plt.yticks(fontsize=30)
# plt.xticks(fontsize=30)
# #plt.xlim(18, 22)
# #plt.ylim(-6,6)
# plt.title(f'{satelitte} - {fecha_datos}',fontsize=37)
# plt.savefig(f'6. Sección de residuos en Escala milimétrica {data_section[0]} - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
# plt.show()
#%%
# from scipy.signal import lombscargle

# w = np.linspace(10, 40, 1000) #frequencies
# a_lomb = lombscargle(day_t_part_full, residuos_mil_scale_full, w, normalize=True)

# plt.plot(w, a_lomb)
# plt.title(f'{satelitte} - {fecha_datos}')
# plt.xlabel('ω [rad/s]')
# plt.ylabel('Amplitud normalizada')
# plt.savefig(f'4. Algoritmo de Lomb - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
# plt.show()

# peak = max(a_lomb)
# for i in range(len(w)):
#     if peak == a_lomb[i]:
#         w_peak = w[i]

# print('## 7 ## - Se ha aplicado el algoritmo de Lomb y se ha calculado el periodo:')
# print(f'    ####   {YEAR}-{MONTH}-{DAY}   ####')
# print(f'    Frecuencia: ω = {w_peak:.4f} rad/s.')
# period = 2*np.pi/w_peak
# print(f'    Periodo: t = {period:.4f} s.')

# #%%
# x = np.arange(0, 50, 0.01)
# y = 6*np.sin(-w_peak*x)

# plt.figure(figsize=(22, 15))
# plt.plot(day_t_part_full, residuos_mil_scale_full, 'o')
# plt.plot(x-0.02, y, linewidth=8)
# plt.xlabel('Tiempo [s]',fontsize=45)
# plt.ylabel('Residuos [mm]',fontsize=45)
# plt.yticks(fontsize=30)
# plt.xticks(fontsize=30)
# #plt.xlim(7, 9)
# #plt.ylim(-6,6)
# plt.title(f'{satelitte} - {fecha_datos}',fontsize=37)
# plt.savefig(f'8. Ajuste sinusoidal a Escala milimétrica {data_section[0]} - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
# plt.show()
