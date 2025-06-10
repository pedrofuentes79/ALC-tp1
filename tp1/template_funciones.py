import numpy as np
import scipy

def construye_adyacencia(D,m): 
    # Función que construye la matriz de adyacencia del grafo de museos
    # D matriz de distancias, m cantidad de links por nodo
    # Retorna la matriz de adyacencia como un numpy.
    D = D.copy()
    l = [] # Lista para guardar las filas
    for fila in D: # recorriendo las filas, anexamos vectores lógicos
        l.append(fila<=fila[np.argsort(fila)[m]] ) # En realidad, elegimos todos los nodos que estén a una distancia menor o igual a la del m-esimo más cercano
    A = np.asarray(l).astype(int) # Convertimos a entero
    np.fill_diagonal(A,0) # Borramos diagonal para eliminar autolinks
    return(A)

def calculaLU(A):
    # matriz es una matriz de NxN
    # Retorna la factorización LU a través de una lista con dos matrices L y U de NxN.
    # Completar! Have fun
    
    n = A.shape[0]
    # inicializo L como la identidad. 
    # Aca me voy a guardar los coeficientes para triangular la matriz U
    L = np.identity(n)
    # inicializo U como la matriz original, 
    # la que voy a triangular y me van a quedar 0s por debajo de la diagonal
    U = A.copy()

    for i in range(n):
        # si el pivot actual es 0, busco la primera fila con valor distinto de 0
        pivot = U[i,i]

        # Genero ceros por debajo de A[i][i]
        for j in range(i+1,n):
            # En L[j][i] guardo el coeficiente por el que multiplico la fila i para obtener la fila j
            coeficiente = U[j,i]/pivot
            L[j,i] = coeficiente

            # Genero ceros por debajo de A[i][i]
            # A la fila j le resto la fila i multiplicada por el coeficiente
            U[j,:] = U[j,:] - coeficiente*U[i,:]
    return L,U


def construye_matriz_K(A):
    # Función que construye la matriz K de grado de los nodos
    # A: Matriz de adyacencia
    # Retorna la matriz K como un numpy.
    
    # Sumo las filas de A para obtener el grado de cada nodo
    suma_filas_A = np.sum(A, axis=1)
    
    # Construyo la matriz diagonal K con los grados
    K = np.diag(suma_filas_A)
    
    return K


def calcula_matriz_C(A): 
    # Función para calcular la matriz de trancisiones C
    # A: Matriz de adyacencia
    # Retorna la matriz C

    # Matriz de grado de los nodos 
    K = construye_matriz_K(A)    
    
    # calculo C como A^T * K^-1
    C = A.T @ np.linalg.inv(K)
    return C

    
def calcula_pagerank(A,alfa):
    # Función para calcular PageRank usando LU
    # A: Matriz de adyacencia
    # alpha: coeficientes de damping
    # Retorna: Un vector p con los coeficientes de page rank de cada museo
    C = calcula_matriz_C(A)
    N = A.shape[0]
    M = (N/alfa) * (np.identity(N) - (1-alfa)*C)

    # Calculamos descomposición PA = LU para M
    L, U = calculaLU(M)

    # El sistema es M p = b => P M p = P b => L U p = P b
    b = np.ones(N)

    # Paso 1: Resolver L y = P b para y (forward substitution)
    y = scipy.linalg.solve_triangular(L, b, lower=True, unit_diagonal=True)

    # Paso 2: Resolver U p = y para p (backward substitution)
    p = scipy.linalg.solve_triangular(U, y, lower=False)

    return p

def calcula_matriz_C_continua(D):
    D_copia = D.astype(float, copy=True)
    np.fill_diagonal(D_copia, 1.0)

    with np.errstate(divide='ignore'): # Ignorar warnings de división por cero si ocurren
        A = 1.0 / D_copia
    # devuelvo la diagonal a 0s, ya que la distancia de un museo a sí mismo es 0
    np.fill_diagonal(A, 0.0)

    # Calcular la suma de las inversas para cada columna
    column_sums = np.sum(A, axis=1)

    K = np.diag(column_sums)
    Kinv = np.linalg.inv(K)
    C = A.T @ Kinv
    return C

def calcula_B(C, cantidad_de_visitas):
    # Recibe la matriz C de transiciones, y calcula la matriz B que representa la relación entre el total de visitas y el número inicial de visitantes
    # suponiendo que cada visitante realizó cantidad_de_visitas pasos
    # C: Matriz de transiciones
    # cantidad_de_visitas: Cantidad de pasos en la red dado por los visitantes. Indicado como r en el enunciado
    # Retorna:Una matriz B que vincula la cantidad de visitas w con la cantidad de primeras visitas v

    # inicializamos B como la matriz identidad para que la primera multiplicación sea C
    N = C.shape[0]
    B = np.identity(N)
    C_acumulada = C.copy()

    for i in range(1, cantidad_de_visitas):
        B += C_acumulada

        # Si no es la ultima iteración, calculo la siguiente potencia de C
        if i < cantidad_de_visitas - 1:
            C_acumulada = C_acumulada @ C

    return B