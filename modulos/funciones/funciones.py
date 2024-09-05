'''Funciones Módulo 1 - Análisis Estático Lineal por Diego Alejandro Noriega Barbosa, MSc'''

#%% Librerias
import pandas as pd
import numpy as np
import math

# Definir una función de formato personalizada
def formato_personalizado(x):
    if pd.isna(x):  # Manejar NaN
        return ''
    elif int(x) == x:  # Si el número es entero
        return '{:.0f}'.format(x)
    else:  # Si el número tiene decimales
        return '{:.2f}'.format(x)  # Puedes ajustar el número de decimales aquí

# Aplicar la función a todo el DataFrame
pd.options.display.float_format = formato_personalizado

#%% Función lectura de información

def lectura_informacion(path_excel:str)->dict:
    '''Función para la lectura de información'''

    # Información del documento de excel de entrada
    informacion_excel = pd.ExcelFile(path_excel,engine='openpyxl')

    # Nombre de las hojas de excel
    hojas_documento_excel = informacion_excel.sheet_names

    # Inicialización de diccionario
    informacion_estructural = {}

    for hoja in hojas_documento_excel:
        # Lee la información de cada hoja a partir de la libreria de pandas (como dataframes)
        informacion_estructural[hoja] = pd.read_excel(path_excel,sheet_name=hoja)

    # Retorna el diccionario
    return informacion_estructural

#%% Función matriz de rigidez de elementos
def Ke(elementoTemporal:pd.Series):
    '''Cálcula la matriz de rigides de un elemento'''

    # Datos iniciales
    elemento = elementoTemporal['T']
    A = elementoTemporal['A']
    I = elementoTemporal['I']
    L = elementoTemporal['L']
    theta = elementoTemporal['An']
    E = elementoTemporal['E']


    # Variables calculation
    if elemento == 1:
        alpha = 12; beta = 4; gamma = 6
    
    elif elemento == 4:
        alpha = 0; beta = 0; gamma = 0

    else:
        alpha = 3; beta = 3; gamma = 3

    # Angle transformation 
    theta = math.radians(theta)

    # Stiffness coefficients
    k1 = E*A/L * math.cos(theta)**2 + alpha*E*I/(L**3)*math.sin(theta)**2
    k2 = E*A/L*math.sin(theta)**2 + alpha*E*I/(L**3)*math.cos(theta)**2
    k3 = beta*E*I/L
    k4 = (E*A/L-alpha*E*I/(L**3)) * math.sin(theta)*math.cos(theta)
    k5 = gamma*E*I/(L**2) * math.cos(theta)
    k6 = gamma*E*I/(L**2) * math.sin(theta)

    # Element stiffness
    ke = np.array([ 
        [k1, k4, -k6, -k1, -k4, -k6],
        [k4, k2, k5, -k4, -k2, k5],
        [-k6, k5, k3, k6, -k5, 0.5*k3],
        [-k1, -k4, k6, k1, k4, k6],
        [-k4, -k2, -k5, k4, k2, -k5],
        [-k6, k5, 0.5*k3, k6, -k5, k3]
        ])

    if elemento==2:
        ke[2,:] = 0; ke[:,2] = 0

    elif elemento==3:
        ke[5,:] = 0; ke[:,5] = 0
    
    elif elemento==4:
        ke[5,:] = 0; ke[:,5] = 0; ke[2,:] = 0; ke[:,2] = 0

    return ke

#%% Función matriz de transformación lambda
def matriz_lambda(theta:float):
    '''Cálcula la matriz de transformación lambda'''
    lambda_m = np.zeros([6,6])

    # Angle transformation 
    theta = math.radians(theta)

    # Sin and cos
    c = math.cos(theta) ; s = math.sin(theta)

    # Lambda matrix formation
    lambda_m[np.ix_([0,1,2],[0,1,2])] = np.array([ [c,s,0],[-s,c,0],[0,0,1]])
    lambda_m[np.ix_([3,4,5],[3,4,5])] = np.array([ [c,s,0],[-s,c,0],[0,0,1]])

    return lambda_m

#%% Función vector de fuerzas
def vector_po(cargas:pd.Series,elementoTemporal:dict)->np.array:
    '''Cálculo vector de fuerzas'''
    # Inicialización vector de fuerzas y datos
    po = np.zeros([6,1])

    # Datos iniciales
    elemento = elementoTemporal['T']
    A = elementoTemporal['A']
    I = elementoTemporal['I']
    L = elementoTemporal['L']
    theta = elementoTemporal['An']
    E = elementoTemporal['E']

    # Postensado
    b = -2/L *math.sqrt(cargas['e1']-cargas['e2'])*(math.sqrt(cargas['e1']-cargas['e2']) + math.sqrt(cargas['e3']-cargas['e2']))
    a = (cargas['e3']-cargas['e1']-b*L)/(L**2)
    c = cargas['e1'] 


    if elemento == 3:
        po[1,0] = 5*cargas['wy']*L/8 
        po[2,0] = cargas['wy']*L**2/8 
        po[4,0] = 3*cargas['wy']*L/8 

    return po

def vector_Po(m_lambda:np.array,po:np.array)->np.array:
    '''Cálculo vector Po'''
    # Matriz lambda inversa
    lambda_transpuesta = m_lambda.T

    return lambda_transpuesta@po


def division_matrices(GdLs:int,GdLr:int,matriz_rigidez_estructura:np.array)->dict:
    '''División de matrices'''

    matrices_k = {}

    matrices_k['kll'] = matriz_rigidez_estructura[0:GdLs,0:GdLs]
    matrices_k['klr'] = matriz_rigidez_estructura[0:GdLs,GdLs:]
    matrices_k['krr'] = matriz_rigidez_estructura[GdLs:,GdLs:]
    matrices_k['krl'] = matriz_rigidez_estructura[GdLs:,0:GdLs]

    return matrices_k
