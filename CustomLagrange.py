import numpy as np


def lagrange(X, Y, N: int, Xp):
    '''
    Función interpoladora de Lagrange. Se emplea la formulación de Lagrange 
    para interpolar un punto en una lista de coordenadas deseada. Para ello
    se aporta la coordenada x del punto a interpolar, el grado deseado del 
    polinomio de Lagrange que se empleará y las coordenadas x e y de los
    puntos iniciales que se emplean para interpolar.
    
    Parameters
    ----------
    X : array[float]
        Vector con las Coordenadas en x de los Puntos iniciales
    Y : array[float]
        Vector con las Coordenadas en y de los Puntos iniciales
    N : int
        Grado del polinomio interpolador de Lagrange
    XP : array[float]
        Coordenada en x del Punto a interpolar

    Returns
    -------
    YP : array[float]
        Coordenada en y del Punto interpolado
    '''
    NPOINT = N + 1
    NMAX = len(X)
    Yp = 0

    # Nos aseguramos que el punto a interpolar esté dentro del intervalo deseado
    if Xp < X[0]:
        Xp = X[0]
    if Xp > X[NMAX - 1]: 
        Xp = X[NMAX - 1]
        
    # Obtiene el valor y el índice del punto en X mayor que XP y más cercano
    index = {}
    for k, val in enumerate(X):
        if(val >= Xp):
            index = k + 1
            break

    I0 = index - NPOINT/2

    if I0 <= 0:
        I0 = 1
    if I0 + N > NMAX:
        I0 = NMAX - N
        
    # Interpolación
    for i in range(NPOINT):
        p = 1
        for j in range(NPOINT):
            if j != i:
                p = p*(Xp-X[int(j-1+I0)])/(X[int(i+I0)]-X[int(j+I0)])
        Yp = Yp + p*Y[int(i - 1 + I0)]
    return Yp
