
import numpy as np
import networkx as nx 
import matplotlib.pyplot as plt

def construye_adyacencia(D,m): 
    D = D.copy()
    l = []     
    for fila in D: # recorriendo las filas, anexamos vectores lógicos
        l.append(fila<=fila[np.argsort(fila)[m]] ) # En realidad, elegimos todos los nodos que estén a una distancia menor o igual a la del m-esimo más cercano
    A = np.asarray(l).astype(int) # Convertimos a entero
    np.fill_diagonal(A,0) # Borramos diagonal para eliminar autolinks
    return(A)


def construir_grafo(m, D, museos):
    A = construye_adyacencia(D, m)
    G = nx.from_numpy_array(A) # Construimos la red a partir de la matriz de adyacencia
    G_layout = {
        i:v for i,v in enumerate
        (zip(museos.to_crs("EPSG:22184").get_coordinates()['x'],museos.to_crs("EPSG:22184").get_coordinates()['y']))
    }
    return G, G_layout

def describePartitions(comunidades, m, A, museos, barrios):
    logSummary(comunidades, museos)
    drawMapWithPartitions(comunidades, m, A, museos, barrios)

def logSummary(comunidades, museos):
    museos['comunidad'] = comunidades
    museos_comunidades = museos.groupby('comunidad').size().reset_index(name='cantidad_museos')
    print(museos_comunidades)

def drawMapWithPartitions(comunidades, m, A, museos, barrios):
    G, G_layout = construir_grafo(m, A, museos)
    fig, ax = plt.subplots(figsize=(15, 15)) # Visualización de la red en el mapa
    barrios.to_crs("EPSG:22184").boundary.plot(color='gray', ax=ax)
    node_colors = ['red' if comunidades[i] == 1 else 'blue' for i, _ in enumerate(G.nodes())]
    nx.draw_networkx(
        G,
        G_layout,
        ax=ax,
        node_size=30,
        node_color=node_colors, # Use the generated color list
        with_labels=False # Ensure no labels are drawn by default draw_networkx
    )
    plt.title(f"Comunidad de Museos (m={m})", fontsize=20)
    plt.show()