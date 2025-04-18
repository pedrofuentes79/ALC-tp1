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
    if A.shape[0] != A.shape[1]:
        raise ValueError("La matriz de entrada debe ser cuadrada")

    P = np.identity(n)
    L = np.identity(n)
    U = A.copy()

    for i in range(n):
        # --- Pivoteo Parcial ---
        # Encontrar fila con máximo valor absoluto en la columna i (desde la diagonal hacia abajo)
        pivot_row_index = i + np.argmax(np.abs(U[i:, i]))

        if i != pivot_row_index:
            # Intercambiar filas en U (optimizado: solo desde columna i en adelante)
            U[[i, pivot_row_index], i:] = U[[pivot_row_index, i], i:]
            # Intercambiar filas en P
            P[[i, pivot_row_index], :] = P[[pivot_row_index, i], :]
            # Intercambiar filas en L (solo la parte ya calculada, columnas 0 a i-1)
            L[[i, pivot_row_index], :i] = L[[pivot_row_index, i], :i]

        # --- Chequeo de singularidad ---
        valor_pivot = U[i, i]
        # Usar una tolerancia pequeña para evitar división por cero
        if np.abs(valor_pivot) < 1e-12:
            raise ValueError("La matriz es singular (o casi singular).")

        # --- Eliminación Gaussiana ---
        # Colocar 1 en la diagonal de L
        L[i, i] = 1.0

        # Calcular multiplicadores y eliminar elementos debajo del pivote
        for j in range(i + 1, n):
            coeficiente = U[j, i] / valor_pivot
            L[j, i] = coeficiente  # Guardar coeficiente en L

            # Actualizar fila j de U (optimizado: solo desde columna i)
            U[j, :] -= coeficiente * U[i, :]
            # Fijar explícitamente a cero por estabilidad numérica
            U[j, i] = 0.0

    return P, L, U

def calcula_matriz_C(A): 
    # Función para calcular la matriz de trancisiones C
    # A: Matriz de adyacencia
    # Retorna la matriz C
    
    # extraigo un arreglo con las sumas de cada fila de A
    suma_filas_A = np.sum(A, axis=1)

    # armo la matriz diagonal K con los valores de sumas_filas_A
    K = np.diag(suma_filas_A)
    Kinv = np.linalg.inv(K)
    
    # calculo C como A^T * K^-1
    C = A.T @ Kinv
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
    P, L, U = calculaLU(M)

    # El sistema es M p = b => P M p = P b => L U p = P b
    b = np.ones(N)
    Pb = P @ b

    # Paso 1: Resolver L y = P b para y (forward substitution)
    y = scipy.linalg.solve_triangular(L, Pb, lower=True, unit_diagonal=True)

    # Paso 2: Resolver U p = y para p (backward substitution)
    p = scipy.linalg.solve_triangular(U, y, lower=False)

    return p

def calcula_matriz_C_continua(D): 
    # Función para calcular la matriz de trancisiones C
    # D: Matriz de distancias
    # Retorna la matriz C en versión continua
    
    D = D.copy()
    np.fill_diagonal(D, 1.0)  # para evitar division por cero
    F = 1/D
    np.fill_diagonal(F,0)

    # Calcula inversa de la matriz K, que tiene en su diagonal la suma por filas de F 
    # Vemos que la suma de las filas de F es lo que esta en el denominador de la formula de Cji
    suma_filas_F = np.sum(F, axis=1)
    
    # Evita division por cero
    # Es decir, deja el valor en 0 si la suma de la fila es 0
    Kinv = np.zeros_like(D)
    indices_sin_cero = np.where(suma_filas_F != 0)[0]

    # Calculo la inversa (1/suma_filas_F) para los indices que no son cero.
    # No hace falta hacer np.linalg.inv porque es una matriz diagonal.
    Kinv[indices_sin_cero, indices_sin_cero] = 1.0 / suma_filas_F[indices_sin_cero]

    # Calcula C multiplicando Kinv y F
    C = Kinv @ F
    return C

def calcula_B(C, cantidad_de_visitas):
    # Recibe la matriz C de transiciones, y calcula la matriz B que representa la relación entre el total de visitas y el número inicial de visitantes
    # suponiendo que cada visitante realizó cantidad_de_visitas pasos
    # C: Matriz de transiciones
    # cantidad_de_visitas: Cantidad de pasos en la red dado por los visitantes. Indicado como r en el enunciado
    # Retorna:Una matriz B que vincula la cantidad de visitas w con la cantidad de primeras visitas v

    # inicializamos B como la matriz identidad para que la primera multiplicación sea C
    B = np.zeros((C.shape[0], C.shape[0]))
    for i in range(1, cantidad_de_visitas):
        B += C**i
    
    return B
