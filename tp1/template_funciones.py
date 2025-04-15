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
    L, U = calculaLU(M) # Calculamos descomposición LU a partir de C y d
    b = np.ones(N) 
    Up = scipy.linalg.solve_triangular(L,b,lower=True) # Primera inversión usando L
    p = scipy.linalg.solve_triangular(U,Up) # Segunda inversión usando U
    return p

def calcula_matriz_C_continua(D): 
    # Función para calcular la matriz de trancisiones C
    # A: Matriz de adyacencia
    # Retorna la matriz C en versión continua
    D = D.copy()
    F = 1/D
    np.fill_diagonal(F,0)
    Kinv = ... # Calcula inversa de la matriz K, que tiene en su diagonal la suma por filas de F 
    C = ... # Calcula C multiplicando Kinv y F
    return C

def calcula_B(C,cantidad_de_visitas):
    # Recibe la matriz T de transiciones, y calcula la matriz B que representa la relación entre el total de visitas y el número inicial de visitantes
    # suponiendo que cada visitante realizó cantidad_de_visitas pasos
    # C: Matirz de transiciones
    # cantidad_de_visitas: Cantidad de pasos en la red dado por los visitantes. Indicado como r en el enunciado
    # Retorna:Una matriz B que vincula la cantidad de visitas w con la cantidad de primeras visitas v
    B = np.eye(C.shape[0])
    for i in range(cantidad_de_visitas-1):
        # Sumamos las matrices de transición para cada cantidad de pasos
        pass
    return B


if __name__ == "__main__":
    # Test de la función calculaLU
    for _ in range(10):
        n = np.random.randint(2, 10)
        A = np.random.rand(n, n)
        L, U = calculaLU(A)
        assert np.allclose(A, L@U)


    # Test de calculaLU con una matriz que fuerce a cambiar el pivot
    #     
    # ESTE TEST FALLA.
    A = np.array([[1,1,3],[1,1,6],[7,8,9]])
    L, U = calculaLU(A)
    print(np.allclose(A,L@U))


    # ESTE TEST ANDA (si manualmente cambiamos las filas de orden)
    A = np.array([[1,1,3],[7,8,9],[1,1,6]])
    L, U = calculaLU(A)
    print(np.allclose(A,L@U))
    

