# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 14:22:18 2023

@author: ruben
"""

import numpy as np
import matplotlib.pyplot as plt
import math

import seaborn as sns;sns.set()
sns.set()

from Extract_Data_Functions import extract_data_CRD
from Extract_Data_Functions import extract_data_CPF
from Extract_Data_Functions import extract_data_station
from CustomLagrange import lagrange
from transformations import ecef2lla
from transformations import geo2topo
from MariniMurrayCorrections import Marini_Murray


### ------ TOPEX ------ ### Estación 7811 --- 2019
CRD_files = ['topex_crd_20190103_01_00.frd', 'topex_crd_20190130_18_00.frd', 'topex_crd_20190208_18_00.frd', \
              'topex_crd_20190227_17_00.frd', 'topex_crd_20190330_00_00.frd', 'topex_crd_20190418_22_00.frd', \
              'topex_crd_20190507_23_00.frd', 'topex_crd_20190529_22_00.frd', 'topex_crd_20190603_22_00.frd', \
              'topex_crd_20190726_01_00.frd', 'topex_crd_20190821_21_00.frd', 'topex_crd_20190831_01_00.frd', \
              'topex_crd_20190912_22_00.frd', 'topex_crd_20191005_19_00.frd', 'topex_crd_20191130_02_00.frd', \
              'topex_crd_20191209_21_00.frd', 'topex_crd_20191220_18_00.frd']
CPF_files = ['topex_tlecpf_190102_5021.aas', 'topex_tlecpf_190129_5291.aas', 'topex_tlecpf_190207_5381.aas', \
              'topex_tlecpf_190226_5571.aas', 'topex_tlecpf_190329_5881.aas', 'topex_tlecpf_190417_6071.aas', \
              'topex_tlecpf_190506_6261.aas', 'topex_tlecpf_190528_6481.aas', 'topex_tlecpf_190602_6531.aas', \
              'topex_tlecpf_190725_7061.aas', 'topex_tlecpf_190820_7321.aas', 'topex_tlecpf_190830_7421.aas', \
              'topex_tlecpf_190911_7541.aas', 'topex_tlecpf_191004_7771.aas', 'topex_tlecpf_191129_8331.aas', \
              'topex_tlecpf_191208_8421.aas', 'topex_tlecpf_191219_8531.aas']


### ------ ENVISAT ------ ### Estación 7810 --- 2017, 2018
# CRD_files = [ 'envisat_crd_20170217_20_00.frd', 'envisat_crd_20170219_20_00.frd', 'envisat_crd_20170410_19_00.frd', \
#               'envisat_crd_20170515_19_00.frd', 'envisat_crd_20170611_19_00.frd', 'envisat_crd_20170619_19_00.frd', \
#               'envisat_crd_20170716_19_00.frd', 'envisat_crd_20170814_20_00.frd', 'envisat_crd_20170903_19_00.frd', \
#               'envisat_crd_20171101_19_00.frd', 'envisat_crd_20171115_19_00.frd', \
#               'envisat_crd_20180124_19_00.frd', 'envisat_crd_20180213_07_00.frd', \
#               'envisat_crd_20180316_19_00.frd', 'envisat_crd_20180329_20_00.frd', 'envisat_crd_20180410_19_00.frd', \
#               'envisat_crd_20180421_19_00.frd', 'envisat_crd_20180427_18_00.frd', 'envisat_crd_20180517_19_00.frd', \
#               'envisat_crd_20180622_19_00.frd', \
#               'envisat_crd_20180826_18_00.frd', 'envisat_crd_20180905_19_00.frd', 'envisat_crd_20180925_18_00.frd', \
#               'envisat_crd_20181021_19_00.frd', 'envisat_crd_20181030_18_00.frd', \
#               'envisat_crd_20181129_18_00.frd', 'envisat_crd_20181213_06_00.frd', 'envisat_crd_20181220_18_00.frd', \
#               'envisat_crd_20181230_19_00.frd']

# CPF_files = [ 'envisat_cpf_170216_5471.aas', 'envisat_cpf_170216_5471.aas', 'envisat_cpf_170409_5991.aas', \
#               'envisat_cpf_170514_6341.aas', 'envisat_cpf_170609_6601.aas', 'envisat_cpf_170616_6671.aas', \
#               'envisat_cpf_170715_6961.aas', 'envisat_cpf_170813_7251.aas', 'envisat_cpf_170831_7431.aas', \
#               'envisat_cpf_171031_8041.aas', 'envisat_cpf_171114_8181.aas', \
#               'envisat_cpf_180123_5231.aas', 'envisat_cpf_180212_5431.aas', \
#               'envisat_cpf_180315_5741.aas', 'envisat_cpf_180328_5871.aas', 'envisat_cpf_180409_5991.aas', \
#               'envisat_cpf_180420_6101.aas', 'envisat_cpf_180426_6161.aas', 'envisat_cpf_180516_6361.aas', \
#               'envisat_cpf_180621_6721.aas', \
#               'envisat_cpf_180825_7371.aas', 'envisat_cpf_180904_7471.aas', 'envisat_cpf_180924_7671.aas', \
#               'envisat_cpf_181020_7931.aas', 'envisat_cpf_181029_8021.aas', \
#               'envisat_cpf_181128_8321.aas', 'envisat_cpf_181211_8451.aas', 'envisat_cpf_181219_8531.aas', \
#               'envisat_cpf_181229_8631.aas']


stations_file: str = 'SLRF2014_POS+VEL_2030.0_180504.snx.txt'
station_code: int = 7811

'''
Código coordenadas Graz, Austria: 7839
Código coordenadas San Fernando, España: 7824
Código coordenadas Matera, Italia: 7941
Código coordenadas Borowiec, Polonia: 7811
'''

#print('## 00 ## - Se han leído como Inputs los nombre de los ficheros:')

print(f'    - Archivo de las estaciones: {stations_file}')
print(f'    - Código de la estación: {station_code}')
print('     ###########')

frecuencias = []
periodos = []
d_year2017 = []
d_year2018 = []
d_year2019 = []
years = []
last_year = 0


for i in range(len(CRD_files)):
    print(f'    - Archivo CRD: {CRD_files[i]}')
    print(f'    - Archivo CPF: {CPF_files[i]}')
    # CRD & CPF data
    
    day_t_CRD, fl_t_CRD, P0, T0, RH, lambda_laser, day_t_cal, obj_d_cal, del_t_cal, \
    MJD, YEAR, MONTH, DAY, satelitte = extract_data_CRD(CRD_files[i], station_code)
    
    if YEAR != last_year:
        years.append(YEAR)
        last_year = YEAR
        
    if YEAR == 2017:
        d_year2017.append(MJD)
    if YEAR == 2018:
        d_year2018.append(MJD)
    if YEAR == 2019:
        d_year2019.append(MJD)
    
    day_t_CPF, position_CPF = extract_data_CPF(CPF_files[i], MJD)
    
    if MONTH < 10:
        fecha_datos = f'{YEAR}-0{MONTH}-{DAY}'
    else:
        fecha_datos = f'{YEAR}-{MONTH}-{DAY}'
    
    #print('## 01 ## - Se han leído los archivos y se han extraído los datos necesarios.')
    
    # Coordenadas geocéntricas del satélite
    
    c = 299792458 # m/s
    day_t_interp = []
    
    range_CRD = []
    
    for j in range(len(day_t_CRD)):
        # Tiempo de interpolación
        day_t_interp.append(day_t_CRD[j] + fl_t_CRD[j]/2)
        
        # Rango al objeto
        range_CRD.append(fl_t_CRD[j]*c/2)
    
    x_geo = position_CPF[0]
    y_geo = position_CPF[1]
    z_geo = position_CPF[2]
    
    #print('## 02 ## - Se han calculado las coordenadas geocéntricas del satélite.')
    
    # Coordenadas interpoladas del satélite (interpolador facilitado - FUNCIONA)
    
    NMAX = len(day_t_CPF)
    N = 9
    
    x_obs = []
    y_obs = []
    z_obs = []
    
    for k in range(len(day_t_CRD)):
        
        X_int = lagrange(day_t_CPF, x_geo, N, day_t_interp[k])
        Y_int = lagrange(day_t_CPF, y_geo, N, day_t_interp[k])
        Z_int = lagrange(day_t_CPF, z_geo, N, day_t_interp[k])
        
        x_obs.append(X_int)
        y_obs.append(Y_int)
        z_obs.append(Z_int)
        
    #print('## 03 ## - Se han determinado las coordenadas interpoladas del satélite \
    #con los datos de tiempo del CPF.')
    
    # Coordenadas geocéntricas de la estación láser donde se toman las medidas
    
    st_position, st_velocity, data_st_epoch = extract_data_station(str(station_code), stations_file)
    x_st = st_position[0]
    y_st = st_position[1]
    z_st = st_position[2]
    
    #print('## 04 ## - Se han determinado las coordenadas de la estación (sin corrección):')
    print(f'    x = {x_st:.3f} m')
    print(f'    y = {y_st:.3f} m')
    print(f'    z = {z_st:.3f} m')
    
    
    # Distancia (rango) de la estación al satélite 
    # r = sqrt[(x_st-x_obs)^2 + (y_st-y_obs)^2 + (z_st-z_obs)^2]
    
    range_int = []
    for l in range(len(x_obs)):
        range_int.append(math.sqrt((x_st - x_obs[l])**2 + \
                                   (y_st - y_obs[l])**2 + \
                                   (z_st - z_obs[l])**2))
    #print('## 05 ## - Se han calculado las distancias de la estación al satélite.')
          
    # Coordenadas topocéntricas
    
    latitud_st, longitud_st, altitud_st = ecef2lla(x_st, y_st, z_st)
    coords_station =  x_st, y_st, z_st, latitud_st, longitud_st, altitud_st
    azimuth, elevation = geo2topo(x_obs, y_obs, z_obs, coords_station, range_int)
    
    #print('## 06 ## - Se han calculado las coordenadas topocéntricas de la estación.')
    
    # Correcciones (Extended ranging equation) - Marini Murray
    '''
    d = (1/2)cΔt + Δd_0 + Δd_s + Δd_b + Δd_r + η
    
    (1/2)cΔt: = range_CRD
    Δd_0: corrección de excentricidad en tierra = 0 (no se aplica) χ
    Δd_s: corrección de excentricidad del satélite ?
    Δd_b: retardo de la estación --> CRD (40) ?
    Δd_r: corrección modelo Marini Murray = corr_range_MM
    η: errores observacionales y sistemáticos restantes = (no se aplica) χ
    
    '''
    P0 = 962.50 
    T0 = 294.46
    RH = 57.8
    lambda_laser = 0.532 # micrometros
    
    corr_range_MM = []
    for n in range(len(elevation)):
        corr_range_MM.append(Marini_Murray(latitud_st, altitud_st, elevation[n], \
                                           P0, T0, RH, lambda_laser))
    
    #print('## 07 ## - Se ha aplicado la corrección de Marini Murray.')
    # Residuos
    
    residuos = []
    corr_residuos = []
    for r in range(len(range_CRD)):
        residuos.append(range_CRD[r] - range_int[r])
        corr_residuos.append(residuos[r] - corr_range_MM[r])
    

    
    #print('## 08 ## - Se han calcuado los residuos (plot 1).')
    ##################
    
    
    # Eliminación de la tendencia - Método 1
    # from scipy import signal
      
    # residuos_detrend =  signal.detrend(corr_residuos)
    # plt.plot(day_t_CRD, residuos_detrend, '.')
    # plt.title(f'Residuos: {fecha_datos}')
    # plt.xlabel('day time/s')
    # plt.ylabel('Residuos/m')
    # plt.show()
    
    # plt.plot(day_t_CRD, residuos_detrend, 'o', day_t_CRD, corr_residuos, 'o')
    # plt.title(f'Residuos con y sin tendencia: {fecha_datos}')
    # plt.xlabel('day time/s')
    # plt.ylabel('Residuos/m')
    # plt.show()   
    #print('## 09 ## - Se ha eliminado la tendencia en los residuos (plot 2).')
    ##############
    
    # Eliminación de la tendencia - Método 2
    # from sklearn.linear_model import LinearRegression
    # model = LinearRegression().fit(np.reshape(day_t_CRD, (-1, 1)), np.reshape(corr_residuos, (-1, 1)))
    # res_pred = model.predict(np.reshape(day_t_CRD, (-1,1))).flatten()
    
    # residuos_detrend = corr_residuos - res_pred
    
    # plt.plot(day_t_CRD, residuos_detrend, '.')
    # plt.title(f'Residuos: {fecha_datos}')
    # plt.xlabel('day time/s')
    # plt.ylabel('Residuos/m')
    # plt.show()
    
    # plt.plot(day_t_CRD, residuos_detrend, 'o', day_t_CRD, corr_residuos, 'o')
    # plt.title(f'Residuos con y sin tendencia: {fecha_datos}')
    # plt.xlabel('day time/s')
    # plt.ylabel('Residuos/m')
    # plt.show()   
    #################
    
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.linear_model import LinearRegression
    
    poly = PolynomialFeatures(degree=2, include_bias=False)
    poly_features = poly.fit_transform(np.reshape(day_t_CRD, (-1, 1)))
    poly_reg_model = LinearRegression()
    poly_reg_model.fit(poly_features, corr_residuos)
    res_pred = poly_reg_model.predict(poly_features)
    
    residuos_detrend = corr_residuos - res_pred
    
    plt.plot(day_t_CRD, corr_residuos, '.', day_t_CRD, res_pred, 'r')
    plt.xlabel('day time/s')
    plt.ylabel('Residuos/m')
    plt.title(f'{i+1}. Residuos sin tratar: {fecha_datos}')
    plt.show()
    
    plt.plot(day_t_CRD, residuos_detrend, '.')
    plt.title(f'Residuos: {fecha_datos}')
    plt.xlabel('day time/s')
    plt.ylabel('Residuos/m')
    plt.show()
    
    plt.plot(day_t_CRD, residuos_detrend, 'o', day_t_CRD, corr_residuos, 'o')
    plt.title(f'Residuos con y sin tendencia: {fecha_datos}')
    plt.xlabel('day time/s')
    plt.ylabel('Residuos/m')
    plt.show()  
    #print('## 09 ##')
    ####################
            
    # Periodo - Algoritmo de Lomb
    
    from scipy.signal import lombscargle
    
    w = np.linspace(0.00001, 2, 1000) #frequencies
    a_lomb = lombscargle(day_t_CRD, residuos_detrend, w, normalize=True)
    
    plt.plot(w, a_lomb)
    plt.title(f'Algoritmo de Lomb: {fecha_datos}')
    plt.xlabel('w [rad/s]')
    plt.ylabel('Amplitud normalizada')
    plt.show()
    
    peak = max(a_lomb)
    for i in range(len(w)):
        if peak == a_lomb[i]:
            w_peak = w[i]  
    
    #print('## 10 ## - Se ha aplicado el algoritmo de Lomb y se ha calculado el periodo:')
    print(f'    ####   {YEAR}-{MONTH}-{DAY}   ####')
    print(f'    Frecuencia: w = {w_peak:.4f} rad/s.')
    period = 2*np.pi/w_peak
    print(f'    Periodo: t = {period:.4f} s.')
    
    frecuencias.append(w_peak)
    periodos.append(period)

#%%
'''
for j in years:
    if years[j] == 2019:
        for i in range(len(dias)):
            dias[i] = dias[i] - 58483
    
    if years[j] == 2018:
        for i in range(len(dias)):
            dias[i] = dias[i] - 58118
    if years[j] == 2017:
        for i in range():
'''
dias_year = []
for i in range(len(d_year2017)):
    d_year2017[i] = d_year2017[i] - 57753
    dias_year.append(d_year2017[i])
for i in range(len(d_year2018)):
    d_year2018[i] = d_year2018[i] - 58118
    if years[0] == 2017:
        d_year2018[i] += 365
    dias_year.append(d_year2018[i])
for i in range(len(d_year2019)):
    d_year2019[i] = d_year2019[i] - 58483
    if years[0] == 2017:
        d_year2018[i] += 730
    if years[0] == 2018:
        d_year2018[i] += 365
    dias_year.append(d_year2019[i])
    
    # TOPEX 2019 --> 58483 
    # ENVISAT 2018 --> 58118
    
with open(f'Evolution {satelitte}.txt', 'w') as f:
    for i in range(len(periodos)):
        f.write(f'{dias_year[i]} {periodos[i]}\n')
    
#%%

t = [0, 365*len(years)]
p = np.polyfit(dias_year, periodos, 1)
pt = np.polyval(p, t)

variation_period = (pt[1]-pt[0])/(365*len(years))

if variation_period < 0:
    print(f'El periodo decrece {abs(variation_period):.5f} s/día')
else:
    print(f'El periodo crece {variation_period:.5f} s/día')

plt.figure(figsize=(10,5)) 


# if YEAR==2019:
#     plt.xticks([1, 1.5, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335], \
#                ['\n\n2019', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Ago', 'Sept', 'Oct', 'Nov', 'Dec'])
# else:
#     plt.xticks([1, 1.5, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, \
#                 366, 366.5, 397, 425, 456, 486, 517, 547, 578, 609, 639, 670, 700], \
#                ['\n\n2017', 'Jan.', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Ago', 'Sept', 'Oct', 'Nov', 'Dec', \
#                 '\n\n2018', 'Jan.', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Ago', 'Sept', 'Oct', 'Nov', 'Dec'])
if YEAR==2019:
    plt.xticks([1, 1.5, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335], \
               ['\n\n2019', 'Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sept', 'Oct', 'Nov', 'Dic'])
else:
    plt.xticks([1, 1.5, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, \
                366, 366.5, 397, 425, 456, 486, 517, 547, 578, 609, 639, 670, 700], \
               ['\n\n2017', 'Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sept', 'Oct', 'Nov', 'Dic', \
                '\n\n2018', 'Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sept', 'Oct', 'Nov', 'Dic'])

plt.plot(dias_year, periodos, 'o', t, pt)
#plt.title(f'Periodos de giro - {satelitte}')
#plt.xlabel('Date', fontsize=15)
plt.ylabel('Periodo [s]', fontsize=15)
#plt.ylim(0, 500)
plt.ylim(8, 12)
plt.xlim(0, 360*len(years)+5)


plt.savefig(f'Periodos de giro - {satelitte}.jpg', dpi=200.0)
plt.show()

#%%
Periodo_medio = sum(periodos)/len(periodos)
sigma = np.sqrt((sum((periodos-Periodo_medio)**2))/len(periodos))
print(f'Periodo medio: {Periodo_medio:.2f} ± {sigma:.2f} s.')
 