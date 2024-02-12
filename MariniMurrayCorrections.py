# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 12:50:05 2023

@author: ruben
"""
import math
import numpy as np

def Marini_Murray(latitud, altitud, elevation, P0, T0, RH, lambda_laser):
    '''
    Función que aplica la corrección troposférica del láser mediante la 
    formulación de Marini y Murray. Con ello se tienen en cuenta, en la 
    propagación del láser, parámetros atmosféricos como la temperatura, la
    presión, la humedad relativa, etc., que pueden tener gran relevancia 
    para ángulos mayores de 10º.
    
             f(𝜆)            A + B
    Δd_r = --------- · ----------------------
            f(𝜙,H)               B/(A + B)
                       sin E + --------------
                                sin E + 0.01
    _________________________________________________________________
    A = 0.002357·P0 + 0.000141·e0
                                                P0^2        2
    B = (1.084×10^−8)·P0·T0·K + (4.734×10^−8)· ------ · ---------
                                                T0      (3 - 1/K)
    K = 1.163 − 0.00968·cos(2𝜙) − 0.00104·T0 + 0.00001435·P0
    _________________________________________________________________
    # Parámetro de frecuencia láser (f_lamda):    
        
                    0.0164   0.000228
    f(𝜆) = 0.9650 + ------ + --------
                     𝜆^2       𝜆^4
                     
    # Función del sitio láser (f_lat):     
        
    f(𝜙,H) = 1 - 0.0026·cos(2𝜙) - 0.00031·H
    
    # Presión de vapor del agua en el lugar del láser (mbar - 100 Pa)
    
          RH               7.5·(T0 - 273.15)
    e0 = ---- · 6.11·10^-----------------------
         100             273.3 + (T0 - 273.15)
    _________________________________________________________________
        
    Parameters
    ----------
    latitud: float
        𝜙: latitud de la estación (rad)
    altitud: float
        altitud de la estación (km) --> H: altitud (m)
    elevation: list[float]
        elevación del satélite (º) --> E: elevación (rad)
    P0: float
        presión atmosférica en el lugar del láser (100 Pa)
    T0: float
        temperatura atmosférica en el lugar del láser (K)
    RH: float
        humedad relativa en el lugar del láser
    lambda_laser: float
        𝜆: longitud de onda del láser (um)
        
    Returns
    -------
    d_r: list[float]
        Δd_r: corrección del rango (m)
        
    '''
    phi = latitud
    H = altitud/1000
    E = elevation*math.pi/180
    
    cos2phi = np.cos(2*phi)
    
    e0 = (RH/100)*6.11*10**((7.5*(T0 - 273.15))/(273.3 + (T0 - 273.15)))
    A = 0.002357*P0 + 0.000141*e0
    K = 1.163 - 0.00968*cos2phi - 0.00104*T0 + 0.00001435*P0
    B = 1.084*10**(-8)*P0*T0*K + 4.734*10**(-8)*(P0**2/T0)*2/(3 - 1/K)
    f_lambda = 0.9650 + (0.0164/lambda_laser**2) + (0.000228/lambda_laser**4)
    f_lat = 1 - 0.0026*cos2phi - 0.00031*H
    
    sinE = np.sin(E)
    d_r = (f_lambda/f_lat)*(A + B)/(sinE + (B/(A + B))/(sinE + 0.01))
    return d_r
    