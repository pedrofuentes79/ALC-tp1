import numpy as np
from scipy.linalg import solve_triangular
from template_funciones import construye_matriz_K
from template_funciones import calculaLU
from typing import List


simetrizar_A = lambda A: (A + A.T) // 2 


def calcula_L(A:np.ndarray) -> np.ndarray:
    """
    Calcula la matríz L, a partir de la matriz A matriz de adyacencia.    
    """
    A = simetrizar_A(A)
    K = construye_matriz_K(A)
    L = K - A
    return L


def construye_matriz_P(A:np.ndarray) -> np.ndarray:
    """
    Construye la matriz P que cumple:
    Pij = k_i * k_j /2E   
    """ 
    k = np.sum(A, axis=1)  # <-- grados
    suma_doble_E = np.sum(k) / 2   # <- constante 2E
    P = np.outer(k, k) / suma_doble_E  
    return P


def calcula_R(A:np.ndarray) -> np.ndarray:
    """
    Calcula la matríz R, a partir de la matriz A, que es una matriz de adyacencia.
    """    
    A = simetrizar_A(A)
    P = construye_matriz_P(A)
    R = A - P
    return R


def calcula_lambda(L,v):
    """
    Calcula la cantidad de conexiones entre los dos grupos.
    Usamos la fórmula:
    λ = 1/4 * s^T * L * s
    Para 
    - L: matriz laplaciana
    - s: vector de partición     
    """
    s = np.sign(v)  # <-- S_i = signo(v_i)
    _lambda = 1/4 * s.T @ L @ s 
    return _lambda


def calcula_Q(R,v):
    """
    Calcula la cantidad de conexiones entre los dos grupos.
    Usamos la fórmula:
    Q = 1/4E * s^T * R * s
    Para
    - R: matriz de conexiones
    - s: vector de partición
    - E: número de conexiones totales
    """
    s = np.sign(v)  # <-- S_i = signo(v_i)
    E = np.sum(R) / 2  # <-- número de conexiones totales
    Q = 1/4/E * s.T @ R @ s 
    return Q    


def metpot(M:np.ndarray, DIST_CONVERGENCIA: float=1e-6, MAX_STEPS=10**6) -> tuple:
    """
    Método de la potencia para encontrar el autovalor de mayor módulo y su correspondiente autovector.
    Parámetros:
    - M: matriz cuadrada de la que se desea encontrar el autovalor y autovector.
    """
    n = M.shape[0]
    x0 = np.random.rand(n) 
    x0 /= np.linalg.norm(x0) 

    for _ in range(MAX_STEPS):
        Ax = M @ x0
        x1 = Ax / np.linalg.norm(Ax)
        if np.allclose(x0, x1, atol=DIST_CONVERGENCIA):
            break
        x0 = x1

    autovector = x0
    autovalor = autovector.T @ M @ autovector  # <-- cociente de Rayleigh

    return autovalor, autovector    


def deflaciona(M: np.ndarray) -> np.ndarray:
    """ 
    Recibe una matriz M, calcule su primer autovector y autovalor, y calcule la matríz M deflacionada.
    """
    autovalor, autovector = metpot(M)
    M_deflacionada = M - autovalor * np.outer(autovector, autovector) / np.linalg.norm(autovector)**2
    
    return M_deflacionada


def metpotI(M: np.ndarray, mu:float) -> tuple:
    """
    Método para obtener el autovalor más chico de M + mu*I y su correspondiente autovector. 
    """
    M2 = M + mu * np.eye(M.shape[0])
    L, U = calculaLU(M2)
    M2_inv = calcula_inversa_con_LU(L, U)    
    autovalor, autovector = metpot(M2_inv)
    
    return 1/autovalor, autovector


def calcula_inversa_con_LU(L: np.ndarray, U: np.ndarray) -> np.ndarray:
    """
    Calcula la inversa de una matriz utilizando su factorización LU.
    """
    n = L.shape[0]
    I = np.eye(n)
    L_inversa = solve_triangular(L, I, lower=True)
    U_inversa = solve_triangular(U, I, lower=False)

    return U_inversa @ L_inversa


def metpotI2(M:np.ndarray, mu:float) -> tuple:
    """
    Método para obtener el segundo autovalor más chico de M + mu*I y su correspondiente autovector.
    """  

    # Calculamos los autovalores y autovectores más chicos de M2 = M + mu*I
    # usando el método de la potencia y luego deflacionamos la matríz para que 
    # el siguiente autovalor más chico no sea el mismo que el anterior.
    n = M.shape[0]
    M2 = M + mu * np.eye(n)

    M2_inv = calcula_inversa_con_LU(*calculaLU(M2))
    M_inv_deflacionada = deflaciona(M2_inv)  
    
    # Aplicamos el método de la potencia a la matriz deflacionada
    autovalor, autovector = metpot(M_inv_deflacionada)

    return 1/autovalor, autovector

def laplaciano_iterativo(A: np.ndarray, niveles: int) -> list:
    """
    Calcula las particiones de una red de forma recursiva usando bisección espectral.
    (Enunciado: "partir la red en dos grupos, y luego repetir el algoritmo en cada
    uno de los subgrupos encontrados.")

    Con niveles=k, se obtienen 2^k particiones.
    """
    
    def particionar(indices_nodos: np.ndarray, nivel_restante: int) -> List[List[int]]:
        """
        Función interna y recursiva que realiza la partición.
        """
        # CASO BASE
        ## Ya llegamos al k-esimo nivel, o la comunidad es trivial.
        if nivel_restante == 0 or len(indices_nodos) <= 1:
            return [list(indices_nodos)] if len(indices_nodos) > 0 else []

        # PASO RECURSIVO
        # 1. Crear el subgrafo para la comunidad actual.
        ## Usamos np.ix_ para armar el subgrafo que solo tiene esos nodos.
        sub_A = A[np.ix_(indices_nodos, indices_nodos)]
        if sub_A.shape[0] < 2: return [list(indices_nodos)]
        
        # 2. Partir la red en dos grupos usando el Laplaciano.
        ## Calculamos el v2, el 2do autovector mas chico de L
        sub_L = calcula_L(sub_A)
        _, v2 = metpotI2(sub_L, 1e-9) # si pasamos mu=0, la matriz podria ser singular! Asi los autovalores son todos != 0


        ## Dividimos los nodos del subgrafo en dos comunidades.
        ## np.where nos devuelve los indices dentro de v2 donde es positivo o negativo.
        indices_subgrafo_com1 = np.where(v2 > 0)[0]
        indices_subgrafo_com2 = np.where(v2 <= 0)[0]
        
        ## Mapeamos los indices locales del subgrafo de vuelta a los indices originales.
        ## Al pasar indices_nodos[indices_subgrafo_com1] estamos devolviendo los indices
        ## de los nodos que pertenecen a la comunidad 1, y lo mismo para 2.
        comunidad1 = indices_nodos[indices_subgrafo_com1]
        comunidad2 = indices_nodos[indices_subgrafo_com2]
        
        # Si una comunidad quedó vacía, devolvemos la comunidad como está
        if len(comunidad1) == 0 or len(comunidad2) == 0:
            return [list(indices_nodos)]

        return particionar(comunidad1, nivel_restante - 1) + \
               particionar(comunidad2, nivel_restante - 1)

    # Primera llamada, le pasamos los indices de todos los nodos del grafo.
    indices_iniciales = np.arange(A.shape[0])
    return particionar(indices_iniciales, niveles)


    
    