from linalg import (
    metpot,
    calcula_inversa_con_LU,
    calculaLU,
    metpotI2,
    laplaciano_iterativo,
    modularidad_iterativo,
)
import numpy as np


def test_met_pot():
    """
    Test para la función metpot.
    """
    A = np.array([[1, 2, 3], [2, 1, 4], [3, 4, 1]])
    eigenvalues, eigenvectors = np.linalg.eig(A)
    _eigenvalue, _eigenvector = metpot(A)
    assert np.allclose(eigenvectors.T[0], _eigenvector)
    assert np.allclose(eigenvalues[0], _eigenvalue)


def test_inversa_con_lu():
    """
    Test para la función que calcula la inversa de una matriz usando LU.
    """
    A = np.array([[1, 2, 3], [0, 1, 4], [5, 6, 0]])
    A_inv_np = np.linalg.inv(A)
    A_inv_lu = calcula_inversa_con_LU(*calculaLU(A))
    assert np.allclose(A_inv_np, A_inv_lu)


def test_met_pot_I2():
    """
    Test para la función metpotI2.
    """
    # Usamos una bien condicionada. Si no, no pasa el test :D
    A = np.array([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 2.0]])
    mu = 0.1
    eigenvalues, eigenvectors = np.linalg.eig(A + mu * np.eye(A.shape[0]))
    assert all(eigenvalues > 0)

    _eigenvalue, _eigenvector = metpotI2(A, mu)

    assert np.isclose(np.abs(np.dot(eigenvectors.T[1], _eigenvector)), 1)
    assert np.allclose(eigenvalues[1], _eigenvalue)


# Este es el grafo de la figura 1 del enunciado
grafo_figura_1 = np.array(
    [
        [0, 1, 1, 1, 0, 0, 0, 0],
        [1, 0, 1, 1, 0, 0, 0, 0],
        [1, 1, 0, 1, 0, 0, 0, 0],
        [1, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 1],
        [0, 0, 0, 0, 1, 0, 1, 1],
        [0, 0, 0, 0, 1, 1, 0, 1],
        [0, 0, 0, 0, 1, 1, 1, 0],
    ]
)


def test_laplaciano_iterativo():
    """
    Test para la función laplaciano_iterativo.
    """
    A = grafo_figura_1
    niveles = 2  # entonces son 4 comunidades
    particiones = laplaciano_iterativo(A, niveles)
    assert len(particiones) <= 2**niveles
    assert all(len(comunidad) > 0 for comunidad in particiones)

    for i, comunidad in enumerate(particiones):
        print(f"Comunidad {i}: {comunidad}")


def test_modularidad_iterativo():
    """
    Test para la función modularidad_iterativo.
    """
    A = grafo_figura_1
    comunidades = modularidad_iterativo(A)
    print(comunidades)


if __name__ == "__main__":
    np.random.seed(99)

    test_met_pot()
    test_inversa_con_lu()
    test_met_pot_I2()
    test_laplaciano_iterativo()
    test_modularidad_iterativo()
    print("Test passed successfully.")
