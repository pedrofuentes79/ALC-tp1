import numpy as np
from scipy.linalg import solve_triangular
from template_funciones import construye_matriz_K
from template_funciones import calculaLU
from typing import List


simetrizar_A = lambda A: np.ceil((A + A.T) / 2)


def calcula_L(A: np.ndarray) -> np.ndarray:
    """
    Calcula la matríz L, a partir de la matriz A matriz de adyacencia.
    """
    A = simetrizar_A(A)
    K = construye_matriz_K(A)
    L = K - A
    return L


def calcula_P(A: np.ndarray) -> np.ndarray:
    """
    Construye la matriz P que cumple:
    Pij = k_i * k_j /2E
    """
    k = np.sum(A, axis=1)  # <-- grados
    suma_doble_E = np.sum(k)  # <- constante 2E
    P = np.outer(k, k) / suma_doble_E
    return P


def calcula_R(A: np.ndarray) -> np.ndarray:
    """
    Calcula la matríz R, a partir de la matriz A, que es una matriz de adyacencia.
    """
    A = simetrizar_A(A)
    P = calcula_P(A)
    R = A - P
    return R


def calcula_lambda(L, v):
    """
    Calcula la cantidad de conexiones entre los dos grupos.
    Usamos la fórmula:
    λ = 1/4 * s^T * L * s
    Para
    - L: matriz laplaciana
    - s: vector de partición
    """
    s = np.sign(v)  # <-- S_i = signo(v_i)
    _lambda = 1 / 4 * s.T @ L @ s
    return _lambda


def calcula_Q(R, v, n_aristas):
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
    Q = 1 / (4 * n_aristas) * s.T @ R @ s
    return Q


def metpot(M: np.ndarray, DIST_CONVERGENCIA: float = 1e-5, MAX_STEPS=10**5) -> tuple:
    """
    Método de la potencia para encontrar el autovalor de mayor módulo y su correspondiente autovector.
    Parámetros:
    - M: matriz cuadrada de la que se desea encontrar el autovalor y autovector.
    """
    n = M.shape[0]
    anterior = np.random.rand(n)
    anterior /= np.linalg.norm(anterior)

    for _ in range(MAX_STEPS):
        Ax = M @ anterior

        actual = Ax / np.linalg.norm(Ax)

        if (
            min(np.linalg.norm(anterior - actual), np.linalg.norm(anterior + actual))
            < DIST_CONVERGENCIA
        ):
            break

        anterior = actual

    autovector = actual
    autovalor = autovector.T @ M @ autovector  # <-- cociente de Rayleigh

    return autovalor, autovector


def deflaciona(M: np.ndarray, autovalor: float, autovector: np.ndarray) -> np.ndarray:
    """
    Recibe una matriz M, su primer autovector y autovalor, y calcule la matríz M deflacionada.
    """
    M_deflacionada = (
        M
        - autovalor * np.outer(autovector, autovector) / np.linalg.norm(autovector) ** 2
    )

    return M_deflacionada


def metpotI(M: np.ndarray, mu: float) -> tuple:
    """
    Método para obtener el autovalor más chico de M + mu*I y su correspondiente autovector.
    """

    # Invertimos la matriz M + mu*I
    M2 = M + mu * np.eye(M.shape[0])
    L, U = calculaLU(M2)
    M2_inv = calcula_inversa_con_LU(L, U)

    # Aplicamos el método de la potencia a la matriz inversa
    autovalor, autovector = metpot(M2_inv)

    return 1 / autovalor, autovector


def calcula_inversa_con_LU(L: np.ndarray, U: np.ndarray) -> np.ndarray:
    """
    Calcula la inversa de una matriz utilizando su factorización LU.
    """
    n = L.shape[0]
    I = np.eye(n)
    L_inversa = solve_triangular(L, I, lower=True)
    U_inversa = solve_triangular(U, I, lower=False)

    return U_inversa @ L_inversa


def metpotI2(M: np.ndarray, mu: float) -> tuple:
    """
    Método para obtener el segundo autovalor más chico de M + mu*I y su correspondiente autovector.
    """

    # Calculamos los autovalores y autovectores más chicos de M2 = M + mu*I
    # usando el método de la potencia y luego deflacionamos la matríz para que
    # el siguiente autovalor más chico no sea el mismo que el anterior.
    n = M.shape[0]
    M2 = M + mu * np.eye(n)

    M2_inv = calcula_inversa_con_LU(*calculaLU(M2))
    autovalor_inv, autovector_inv = metpot(M2_inv)
    M_inv_deflacionada = deflaciona(M2_inv, autovalor_inv, autovector_inv)

    # Aplicamos el método de la potencia a la matriz deflacionada
    autovalor, autovector = metpot(M_inv_deflacionada)

    return (1 / autovalor) - mu, autovector


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
        if sub_A.shape[0] < 2:
            return [list(indices_nodos)]

        # 2. Partir la red en dos grupos usando el Laplaciano.
        ## Calculamos el v2, el 2do autovector mas chico de L
        sub_L = calcula_L(sub_A)
        _, v2 = metpotI2(
            sub_L, 1e-9
        )  # si pasamos mu=0, la matriz podria ser singular! Asi los autovalores son todos != 0

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

        return particionar(comunidad1, nivel_restante - 1) + particionar(
            comunidad2, nivel_restante - 1
        )

    # Primera llamada, le pasamos los indices de todos los nodos del grafo.
    indices_iniciales = np.arange(A.shape[0])
    return particionar(indices_iniciales, niveles)


def buscar_autovalor_positivo(M: np.ndarray) -> float:
    """
    Busca el autovalor positivo de M.
    """
    A = M.copy()

    for _ in range(M.shape[0]):
        if np.linalg.norm(A) < 1e-9:  # Matriz casi nula
            break

        autovalor, v1 = metpot(A)

        if autovalor > 0:
            return autovalor, v1

        A = deflaciona(A, autovalor, v1)

    return -1, None


def modularidad_iterativo(A: np.ndarray) -> List[List[int]]:
    """
    Divide la red en comunidades usando la heurística de modularidad.
    Empieza con una única comunidad con todos los nodos y va dividiendo en comunidades.
    La heuristica de division es: calcular la modularidad de dividir la comunidad actual en dos. 
    Digamos que la comunidad actual es Q1. Dividimos en Q1a y Q1b.
    Para efectivamente guardarnos este split como el mejor dentro de las comunidades existentes, 
    tenemos que chequear que la modularidad de Q1a + la de Q1b sea mayor que la de Q1.
    Es decir mod(Q1a) + mod(Q1b) > mod(Q1). Si es asi, nos guardamos el split.
    Ojo, puede suceder que haya un Q2 cuyos Q2a y Q2b sean un mejor split, 
    es decir, que su aporte a la mejora de modularidad sea mayor que el de Q1.
    Formalmente, esto seria que mod(Q2a) + mod(Q2b) - mod(Q2) > mod(Q1a) + mod(Q1b) - mod(Q1).
    Si esto es asi, entonces Q2 es el mejor split y lo aplicamos.
    
    Repetimos el proceso hasta que no haya ningún split que mejore la modularidad, ya que puede suceder que
    mod(Qa) + mod(Qb) - mod(Q) < 0, (para todo Q comunidad dentro de las comunidades existentes).
    """
    # 1. Pre-cálculos iniciales sobre el grafo completo
    R = calcula_R(A)
    n_aristas = np.sum(A) / 2

    # 2. Estado inicial: una única comunidad con todos los nodos.
    comunidades = [list(range(A.shape[0]))]

    while True:
        mejor_aumento_de_modularidad = 0
        mejor_split_info = None

        # 3. Iterar sobre todas las comunidades actuales para encontrar el mejor split
        for i, comunidad_actual in enumerate(comunidades):
            if len(comunidad_actual) <= 1:
                continue  # las comunidades triviales no las consideramos

            indices_actuales = np.array(comunidad_actual)

            # Busco el subgrafo actual
            sub_R = R[np.ix_(indices_actuales, indices_actuales)]

            # Estoy en la comunidad i. Cual es la modularidad actual? 
            # Sin hacer ningun split, entonces tomo vector de unos.
            vector_unos = np.ones(indices_actuales.shape[0])
            modularidad_actual = calcula_Q(sub_R, vector_unos, n_aristas)

            # Si la comunidad es trivial, no la consideramos
            if sub_R.shape[0] < 2:
                continue

            # La heurística busca el autovalor MÁS POSITIVO.
            # Para ello, deflacionamos la matriz repetidamente hasta encontrar uno,
            # o hasta que hayamos probado todos los autovalores posibles.
            autovalor, v1 = buscar_autovalor_positivo(sub_R)

            # Si encontramos un autovalor positivo, calculamos su delta_Q
            if autovalor > 0:
                # Usamos la sub_R original para el cálculo de Q, pero con el autovector v1 encontrado
                delta_Q = calcula_Q(sub_R, v1, n_aristas)
            else:
                # Si no, significa que no hay divisiones que mejoren la modularidad.
                # pongo -inf para que el aumento de modularidad no se considere.
                delta_Q = -np.inf

            # Si hacer el split mejora la modularidad respecto de no hacerlo, lo guardamos.
            aumento_de_modularidad = delta_Q - modularidad_actual
            if aumento_de_modularidad > mejor_aumento_de_modularidad:
                mejor_aumento_de_modularidad = aumento_de_modularidad

                indices_locales_c1 = np.where(v1 > 0)[0]
                indices_locales_c2 = np.where(v1 <= 0)[0]

                comunidad1 = list(indices_actuales[indices_locales_c1])
                comunidad2 = list(indices_actuales[indices_locales_c2])

                if comunidad1 and comunidad2:
                    mejor_split_info = {
                        "indice_comunidad_a_dividir": i,
                        "nuevas_comunidades": [comunidad1, comunidad2],
                    }

        # 4. Si el mejor split encontrado mejora la modularidad (delta_Q > 0), lo aplicamos.
        if mejor_split_info and mejor_aumento_de_modularidad > 0:
            indice_a_dividir = mejor_split_info["indice_comunidad_a_dividir"]

            # Actualizamos la lista de comunidades
            comunidades.pop(indice_a_dividir)
            comunidades.extend(mejor_split_info["nuevas_comunidades"])
        else:
            # 5. Si no se encontró ningún split que mejore Q, cortamos.
            break

    return comunidades
