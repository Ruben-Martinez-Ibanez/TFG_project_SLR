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
# data_section = [100, 150]
# station_code: int = 7839

# 2013-07-14
CRD_file: str = '7839_CRD_20130714_2106.frd'
CPF_file: str = 'envisat_cpf_130711_6921.dlr'
data_section = [100, 150]
station_code: int = 7839

#######
 
# 2018-01-31 ### 7839
# CRD_file: str = 'envisat_crd_20180131.frd'
# CPF_file: str = 'envisat_cpf_180125_5251.dlr'
# data_section = [340, 390]
# station_code: int = 7839

# 2018-01-29 ### 7839
# CRD_file: str = 'envisat_crd_20180129.frd'
# CPF_file: str = 'envisat_cpf_180125_5251.dlr'
# data_section = [100, 150]
# station_code: int = 7839

# 2019-04-27 ### 7839
# CRD_file: str = 'envisat_crd_20190427.frd'
# CPF_file: str = 'envisat_cpf_190427_6171.dlr'
# data_section = [65460, 65480]
# station_code: int = 7839

# 2018-09-25 ### 7810
# CRD_file: str = 'envisat_crd_20180925_18_00.frd'
# CPF_file: str = 'envisat_cpf_180924_7671.aas'
# data_section = [67680, 67750]
# station_code: int = 7810

### ------ TOPEX ------ 
# 2019-01-03 ### 7811
# CRD_file: str = 'topex_crd_20190103_01_00.frd'
# CPF_file: str = 'topex_tlecpf_190102_5021.aas'
# station_code: int = 7811

# 2022-03-20 ### 7941
# CRD_file: str = 'topex_crd_20220320_2145_00.frd'
# CPF_file: str = 'topex_cpf_220318_5771.aas'
# data_section = [78459, 78463]
# station_code: int = 7841

#2022-03-20 ### 7810
# CRD_file: str = 'envisat_crd_20170814_20_00.frd'
# CPF_file: str = 'envisat_cpf_170813_7251.aas'
# data_section = [73250, 73320]
# station_code: int = 7810

### ------ HY2A ------
# 2019-01-08 ### 7811
# CRD_file: str = 'hy2a_20200910.frd'
# CPF_file: str = 'hy2a_cpf_200908_7521.sha'
# station_code: int = 7839

### ------ TECHNOSAT ------
# 2023-04-07 ### 7840/7841
# CRD_file: str = 'technosat_20230501.fr2'
# CPF_file: str = 'technosat_cpf_230501_12101.dlr'
# station_code: int = 7841

### ------ IRNSS-1A ------
# 2023-04-07 ### 7827
# CRD_file: str = 'irnss1a_20221227.frd'
# CPF_file: str = 'irnss1a_cpf_221226_36001.isr'
# station_code: int = 7827



stations_file: str = 'SLRF2014_POS+VEL_2030.0_180504.snx.txt'

'''
Código coordenadas Graz, Austria: 7839
Código coordenadas San Fernando, España: 7824
Código coordenadas Matera, Italia: 7941
Código coordenadas Borowiec, Polonia: 7811
Código coordenadas 	Herstmonceux, Reino Unido: 7840
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


#%% Coordenadas geocéntricas de la estación corregidas
# Falta sacar el tiempo de los datos de la estación
# MDJ_st

# Velocidades en metros/año

# vx_st = st_velocity[0]
# vy_st = st_velocity[1]
# vz_st = st_velocity[2]


# t = (MJD[0]-MDJ_st)/365

# x_st +=  vx_st*t
# y_st +=  vy_st*t
# z_st +=  vz_st*t

# print('Coordenadas de la estación corregidas:')
# print(f'x = {x_st} m')
# print(f'y = {y_st} m')
# print(f'z = {z_st} m')
# print('')
# print('## 05 ##')

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

time0 = []
mint = min(day_t_CRD)
for i in range(len(range_CRD)):
    residuos.append(range_CRD[i] - range_int[i])
    corr_residuos.append(residuos[i] - corr_range_MM[i])
    time0.append(day_t_CRD[i]-mint)

plt.figure(figsize=(10,7)) 
plt.plot(day_t_CRD, corr_residuos, 'o')
plt.xlabel('Tiempo [s]', fontsize=20)
plt.ylabel('Residuos [m]', fontsize=20)
plt.title(f'{satelitte} - {fecha_datos}', fontsize=17)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig(f'1. Residuos sin tratar - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
plt.show()

print('## 08 ## - Se han calcuado los residuos (plot 1).')

#%% Eliminación de la tendencia - Método 1
# from scipy import signal
  
# residuos_detrend =  signal.detrend(corr_residuos)
# plt.plot(day_t_CRD, residuos_detrend, '.')
# plt.title(f'Residuos: {fecha_datos}')
# plt.xlabel('day time/s')
# plt.ylabel('Residuos/m')
# plt.savefig(f'2. Residuos sin tendencia - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
# plt.show()

# plt.plot(day_t_CRD, residuos_detrend, 'o', day_t_CRD, corr_residuos, 'o')
# plt.title(f'Residuos con y sin tendencia: {fecha_datos}')
# plt.xlabel('day time/s')
# plt.ylabel('Residuos/m')
# plt.savefig(f'3. Residuos con y sin tendencia - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
# plt.show()   

# print('## 09 ## - Se ha eliminado la tendencia en los residuos (plot 2).')

#%% Eliminación de la tendencia - Método 2 - Regresión lineal
# from sklearn.linear_model import LinearRegression
# model = LinearRegression().fit(np.reshape(day_t_CRD, (-1, 1)), np.reshape(corr_residuos, (-1, 1)))
# res_pred = model.predict(np.reshape(day_t_CRD, (-1,1))).flatten()

# residuos_detrend = corr_residuos - res_pred

# plt.plot(day_t_CRD, residuos_detrend, '.')
# plt.title(f'Residuos: {fecha_datos}')
# plt.xlabel('day time/s')
# plt.ylabel('Residuos/m')
# plt.savefig(f'2. Residuos sin tendencia - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
# plt.show()

# plt.plot(day_t_CRD, residuos_detrend, 'o', day_t_CRD, corr_residuos, 'o')
# plt.title(f'Residuos con y sin tendencia: {fecha_datos}')
# plt.xlabel('day time/s')
# plt.ylabel('Residuos/m')
# plt.savefig(f'3. Residuos con y sin tendencia - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
# plt.show()   

# print('## 09 ##')

#%% Eliminación de la tendencia - Método 3 - Regresión polinómica

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

poly = PolynomialFeatures(degree=4, include_bias=False)
poly_features = poly.fit_transform(np.reshape(day_t_CRD, (-1, 1)))
poly_reg_model = LinearRegression()
poly_reg_model.fit(poly_features, corr_residuos)
res_pred = poly_reg_model.predict(poly_features)

residuos_detrend = corr_residuos - res_pred
plt.figure(figsize=(10,7)) 
plt.plot(day_t_CRD, corr_residuos, 'o', day_t_CRD, res_pred, 'r')
plt.xlabel('Tiempo [s]', fontsize=20)
plt.ylabel('Residuos [m]', fontsize=20)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.title(f'{satelitte} - {fecha_datos}', fontsize=17)
plt.show()

plt.figure(figsize=(10,7)) 
plt.plot(day_t_CRD, residuos_detrend, 'o')
#plt.title(f'{satelitte} - {fecha_datos}', fontsize=17)
plt.xlabel('Tiempo [s]', fontsize=20)
plt.ylabel('Residuos [m]', fontsize=20)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
#plt.ylim(-2, 2)
plt.savefig(f'2. Residuos sin tendencia - {satelitte} - {fecha_datos}.jpg', dpi=200.0)
plt.show()

plt.figure(figsize=(10,7)) 
plt.plot(time0, residuos_detrend, 'o', time0, corr_residuos, 'o')
#plt.title(f'{satelitte} - {fecha_datos}', fontsize=17)
plt.xlabel('Time [s]', fontsize=20)
plt.ylabel('Range residuals [m]', fontsize=20)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
#plt.ylim(-8, 7)
plt.savefig(f'3. Residuos con y sin tendencia - {satelitte} - {fecha_datos}.jpg', dpi=150.0)
plt.show()

print('## 09 ##')

#%% Periodo - Algoritmo de Lomb

from scipy.signal import lombscargle

w = np.linspace(0.00001, 0.4, 1000) #frequencies
a_lomb = lombscargle(day_t_CRD, residuos_detrend, w, normalize=True)

plt.figure(figsize=(10,5)) 
plt.plot(w, a_lomb)
#plt.title(f'{satelitte} - {fecha_datos}', fontsize=17)
plt.xlabel('ω [rad/s]', fontsize=20)
plt.ylabel('Normalised amplitude', fontsize=20)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
#plt.xlim(0, 0.4)
plt.savefig(f'4. Algoritmo de Lomb - {satelitte} - {fecha_datos}.jpg', dpi=100.0)
plt.show()

peak = max(a_lomb)
for i in range(len(w)):
    if peak == a_lomb[i]:
        w_peak = w[i]

print('## 10 ## - Se ha aplicado el algoritmo de Lomb y se ha calculado el periodo:')
print(f'    ####   {YEAR}-{MONTH}-{DAY}   ####')
print(f'    Frecuencia: ω = {w_peak:.4f} rad/s.')
period = 2*np.pi/w_peak
print(f'    Periodo: t = {period:.4f} s.')

#%% Estimación del periodo a partir de datos del 2013
# Datos usados del artículo publicado en Diciembre de 2014 que tiene como título:
# 'Attitude and Spin Period of Space Debris Envisat Measured by Satellite Laser Ranging'

# Datos del artículo
t = 134.74 #s
YEAR_I = 2013
MONTH_I = 9
DAY_I = 25

#rate
r = 0.0367 #s/day

D_YEAR = YEAR-YEAR_I-1
DAYS_YEAR_I = (MONTH_I-12)*30 + 30 - DAY_I
DAYS_YEAR_F = (MONTH-1)*30 + DAY


#Cálculo del periodo
T = t + r*(D_YEAR*365 + DAYS_YEAR_I + DAYS_YEAR_F)
print(f'## 10.5 ## - Se ha estimado el periodo para {YEAR}-{MONTH}-{DAY} \
mediante datos del artículo:')
print(f'Periodo estimado: t` = {T:.3f} s')
print('ESTO ES SOLO PARA ENVISAT')


#%% Ángulos de incidencia - TEST
'''
# A = H·sin(IA)
A = [] #Amplitud
IA = [] #incident angle
H = 2.5 #m   Offset
res_ant = 0
res_sig = 0
for i in range(len(residuos_detrend)):
    #if res_i > res
    A.append(abs(residuos_detrend[i]))
    IA.append(np.arcsin(A[i]/H)*180/np.pi)

plt.plot(IA, residuos_detrend)
plt.show()

plt.plot(day_t_CRD, IA)
'''
#%% 
#picos
'''
#from scipy.signal import find_peaks

# residuos_abs = []
# for i in range(len(residuos_detrend)):
#     residuos_abs.append(abs(residuos_detrend[i]))

#curva = np.polyfit(x, y, deg)

#plt.plot(day_t_CRD, residuos_abs, '.')
tt = np.linspace(66900,67300, 10000)
residuos_p = np.polyfit(day_t_CRD, residuos_detrend, 10)
residuos_curve = np.polyval(residuos_p, tt)

plt.plot(tt, residuos_curve, '.', day_t_CRD, residuos_detrend, '.')
'''

#%%
# LASER vector
x_laser = []
y_laser = []
z_laser = []

for i in range(len(x_obs)):
    x_laser.append(x_obs[i]-x_st)
    y_laser.append(y_obs[i]-y_st)
    z_laser.append(z_obs[i]-z_st)

print('## 11 ## - ')

#%% OscilacionesEscalaMilimétrica

# residuos = []
# corr_residuos = []
# min_day_t_part = math.inf
# for i in range(len(range_CRD)):
#     if day_t_CRD[i] < min_day_t_part:
#         min_day_t_part = day_t_CRD[i]
#     day_t_CRD[i] = day_t_CRD[i]-min_day_t_part
#     residuos.append(range_CRD[i] - range_int[i])
#     corr_residuos.append(residuos[i] - corr_range_MM[i])

# min_day_t_part = math.inf

# #data_subsection = [min(data_section) + i*10, min(data_section) + (i+1)*10]
# day_t_part = []
# residuos_part = []

# for i in range(len(day_t_CRD)):
#     if day_t_CRD[i] > data_section[0]:
#         print(day_t_CRD[i], data_section[0])
#         if day_t_CRD[i] < data_section[1]:
#             if day_t_CRD[i] < min_day_t_part:
#                 min_day_t_part = day_t_CRD[i]
#             day_t_part.append(day_t_CRD[i]-min_day_t_part)
#             #day_t_part_full.append(day_t_CRD[i]-min_day_t_part)
#             residuos_part.append(corr_residuos[i])
#             #r.write(f'{day_t_CRD[i]} {corr_residuos[i]}\n')

# polinomio = np.polyfit(day_t_part, residuos_part, 5)
# residuos_pol = np.polyval(polinomio, day_t_part)

# residuos_mil_scale = []
# for i in range(len(day_t_part)):
#     residuos_mil_scale.append((residuos_pol[i]-residuos_part[i])*1000)
#     #residuos_mil_scale_full.append((residuos_pol[i]-residuos_part[i])*1000)
# with open('Datos_Residuos_escala_milimétrica.txt', 'w') as r:
#     r.write('Segundo_del_día Residuos/mm\n')
#     for i in range(len(day_t_part)):
#         r.write(f'{day_t_part[i]} {residuos_mil_scale[i]}\n')


# plt.plot(day_t_part, residuos_part, 'o',  day_t_part, residuos_pol, '.')
# plt.xlabel('day time/s')
# plt.ylabel('Residuos/m')
# plt.title(f'Sección escala métrica: {fecha_datos}')
# #plt.savefig(f'5. Sección de residuos en Escala métrica - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
# plt.show()
    
# #residuos_mil_scale_detrend = signal.detrend(residuos_mil_scale)
    

# plt.figure(figsize=(15,10))
# plt.plot(day_t_part, residuos_mil_scale, '.')
# plt.xlabel('day time/s',fontsize=30)
# plt.ylabel('Residuos/mm',fontsize=30)
# plt.yticks(fontsize=20)
# plt.xticks(fontsize=20)
# #plt.ylim(-10, 10)
# plt.title(f'Sección escala milimétrica: {fecha_datos}',fontsize=30)
# #plt.savefig(f'6. Sección de residuos en Escala milimétrica - {satelitte} - {fecha_datos}.jpg', dpi=500.0)
# plt.show()



    
