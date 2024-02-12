# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 15:57:32 2022

@author: ruben
"""
import numpy as np

#%% CRD (Consolidated laser Range Data format)

def extract_data_CRD(file_name: str, station_code:int):
    '''
    Función para extraer los datos necesarios del archivo CRD 
    (Consolidated laser Range Data format)
    
    Parameters
    ----------
    file_name : str
        Nombre del archivo CRD
    station_code : int
        DESCRIPTION.

    Returns
    -------        
    day_t_CRD: array[float]
        Day time - Momento del día en que se han tomado las medidas, en segundos.
    fl_t_CRD: array[float]
        Flight time - Tiempo de vuelo del láser
    day_t_cal: array[float]
        Day time calibration - NO SE USA
    obj_d_cal: array[float]
        Object distance calibration - NO SE USA
    del_t_cal: array[float]
        Delayed time calibration - NO SE USA
    P0: float
        Presión atmosférica en el lugar del láser (100 Pa)
    T0: float
        Temperatura atmosférica en el lugar del láser (K)
    RH: float
        Humedad relativa en el lugar del láser
    lambda_laser: float
        Longitud de onda del láser (um)
    MJD: float
        Modified Julian Date - Día Juliano Modificado
    YEAR: int
        Año en el que se tomaron los datos usados
    MONTH: int
        Mes en el que se tomaron los datos usados
    DAY: int
        Día en el que se tomaron los datos usados

    '''
    #ruta_de_acceso = 'C:/Users/ruben/OneDrive/Documentos/4ºFísica/2ºCuatrimestre/TFG/Trabajo/SLR_TFG/Archivos CRD y CPF/'
    ## DATOS de las MEDIDAS##
    # Tiempo del día (en segundos):
    day_time: np.array[float] = []
    # Tiempo de vuelo (en segundos):
    flight_time: np.array[float] = []
    
    ## DATOS ATOMSFÉRICOS
    # Presión superficial
    press: np.array[float] = []
    # Temperatura superficial
    temp: np.array[float] = []
    # Humedad relativa
    rel_hum: np.array[float] = []
    
    
    ##DATOS de CALIBRACIÓN
    # Tiempo del día (en segundos):
    day_time_calibration: np.array[float] = []
    # Distancia al blanco (en metros):
    object_distance_calibration: np.array[float] = []
    # Retardo medido en la calibración (en picosegundos):
    delayed_time_calibration: np.array[float] = [] 
    
    P0 = 962.50 
    T0 = 294.46
    RH = 57.8
    lambda_laser = 0.532 # micrometros
    MJD = 0
    YEAR, MONTH, DAY = 0, 0, 0
    satellite: str = 'satellite_name'
    
    ACTION = 0
    try:
        with open(file_name) as input_file:
            with open('CRD_test.txt', 'w') as f:
                for line in input_file:
                    parameters = line.split()
                    
                    if parameters[0].upper() == 'H2':
                        if int(parameters[2]) == station_code:
                            ACTION = 1
                            
                    if ACTION == 1:
                        if parameters[0].upper() == 'C2':
                            lambda_laser = float(parameters[4])/1000 #micrómetros
                        
                        # MEDIDAS
                        if parameters[0] == '10':
                            day_time.append(float(parameters[1]))
                            flight_time.append(float(parameters[2]))
                            f.write(f'{float(parameters[0])} {float(parameters[1])} {float(parameters[2])}\n')                        
                        
                        # DATOS ATMOSFÉRICOS
                        if parameters[0] == '20':
                            press.append(float(parameters[2]))
                            temp.append(float(parameters[3]))
                            rel_hum.append(float(parameters[4]))
                            
                        # CALIBRACIÓN
                        if parameters[0] == '40':
                            day_time_calibration.append(float(parameters[1]))
                            #object_distance_calibration.append(float(parameters[6]))
                            delayed_time_calibration.append(float(parameters[7]))
                        
                        # SATÉLITE
                        if parameters[0].upper() == 'H3':
                            satellite = parameters[1]
                            
                        # FECHA
                        if parameters[0].upper() == 'H4':
                            YEAR = int(parameters[2])
                            MONTH = int(parameters[3])
                            DAY = int(parameters[4])
                            
                            
                        if parameters[0].upper() == 'H8':
                            ACTION = 0
                            break
        
        # Día Juliano Modificado
        M = MONTH
        Y = YEAR
        if M == 1 or M == 2:
            M = MONTH + 12
            Y = YEAR - 1
            
        A = int(Y/100)
        B = 2 - A + int(A/4)
        JD = int(365.25*(Y + 4716)) + int(30.6001*(M + 1)) + DAY + B - 1524.5
        MJD = JD - 2400000.5
        
        
    except IOError:
           print("There's no file CRD")
           
    return day_time, flight_time, P0, T0, RH, lambda_laser, day_time_calibration, \
           object_distance_calibration, delayed_time_calibration, MJD, YEAR, MONTH, \
           DAY, satellite

#%% CPF (Consolidated Prediction Format)

def extract_data_CPF(file_name: str, MJD:float):
    '''    
    Función para extraer los datos necesarios del archivo CPF
    (Consolidated Prediction Format)    
    
    Parameters
    ----------
    file_name : str
        Nombre del archivo CPF
    MJD : float
        Modified Julian Date - Día Juliano Modificado

    Returns
    -------
    day_t_CPF: array[float]
        Day time - Momento del día en que se han tomado las medidas, en segundos.
    position_CPF: array[array[float]]
        Matriz de posiciones calculadas. Es un vector que contiene 3 vectores, 
          uno por cada coordenada geocéntrica.

    '''
    # Tiempo del día (en segundos):
    day_time: np.array[float] = []
    # Posición x:
    x: np.array[float] = []
    # Posición y:
    y: np.array[float] = []
    # Posición z:
    z: np.array[float] = []
    position: np.array = []

    try:
        with open(file_name) as input_file:
            with open('CPF_test.txt', 'w') as f:
                for line in input_file:
                    parameters = line.split()
                    if parameters[0] == '10':
                        if float(parameters[2]) == MJD:
                            #print(parameters[2])
                            day_time.append(float(parameters[3]))
                            x.append(float(parameters[5]))
                            y.append(float(parameters[6]))
                            z.append(float(parameters[7]))
                            f.write(f'{float(parameters[2])} {float(parameters[3])} {float(parameters[4])} {float(parameters[5])} {float(parameters[6])} {float(parameters[7])}\n')

        position: np.array = [x, y, z]
    except IOError:
        print("There's no file CPF")
    return day_time, position

#%% STATIONS

def extract_data_station(station_code:str, file_name: str):
    '''
    Función para extraer los datos necesarios de la estación

    Parameters
    ----------
    station_code : str
        Código específico de la estación que ha medido los datos empleados
    file_name : str
        Nombre del archivocon la información de las distintas estación de SLR

    Returns
    -------
    [x, y, z]: list[float]
        Lista con las coordenadas de posición de la estación.
    [vx, vy, vz]: list[float]
        Lista con las coordenadas de velocidad de la estación 
        (para fines de corrección)
    Mean_data_epoch : str
        ¡¡¡No sé para qué sirve!!!

    '''
    with open(file_name) as input_file:
        Type_Data: int = 0
        x, y, z = 0, 0, 0
        vx, vy, vz = 0, 0, 0
        
        for line in input_file:

            parameters = line.split()
            
            if parameters[0] == '+SOLUTION/EPOCHS':
                Type_Data = 1
                
            if parameters[0] == '+SOLUTION/ESTIMATE':
                Type_Data = 2
                
            if Type_Data == 1:
                if parameters[0] == station_code:
                    Mean_data_epoch = parameters[6]
                    print(Mean_data_epoch)
                    Type_Data = 0
                    
            if Type_Data == 2:
                try:
                    if parameters[2] == station_code:
                        if parameters[3] == 'A':
                            if parameters[1] == 'STAX':
                                x = float(parameters[8])
                            elif parameters[1] == 'STAY':
                                y = float(parameters[8])
                            elif parameters[1] == 'STAZ':
                                z = float(parameters[8])
                            elif parameters[1] == 'VELX':
                                vx = float(parameters[8])
                            elif parameters[1] == 'VELY':
                                vy = float(parameters[8])
                            elif parameters[1] == 'VELZ':
                                vz = float(parameters[8])
                                Type_Data = 0
                except:
                    None
    return [x, y, z], [vx, vy, vz], Mean_data_epoch


