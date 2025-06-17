## Ejercicio 5

#### <u>Introducción: </u>
En la primera parte nos interesaba ver de que manera pesaba la importancia de cada museo, lo cual se medía en cantidad de visitas, esto a su vez reflejaba de que forma los museos vecinos se iban a ver impactados por el peso del museo en cuestión y viceversa. 

En el segundo TP, buscamos agrupar a los museos a partir de una asignación óptima que indique a qué comunidad pertenece cada uno, empleando métodos espectrales basados en autovalores y autovectores.

#### <u>¿Qué métodos usamos y qué vimos?</u> 

##### ***TP 1***

Para el cálculo del ranking se incorporó un factor $\alpha$ que representa la probabilidad de que un visitante acceda a un museo alejado del circuito habitual. Se observó que tanto la cantidad de vecinos (museos cercanos) como el valor de $\alpha$ influyen significativamente en el ranking: al aumentar los vecinos, el ranking del museo mejoraba; en cambio, al incrementar $\alpha$, el ranking tendía a distribuirse más equitativamente, reduciendo la concentración en museos céntricos. En síntesis, contar con más vecinos incrementa la probabilidad de visita y, por ende, mejora el posicionamiento del museo.

Además, se comprobó que, con una matriz de transición bien condicionada, es posible estimar de forma razonable la cantidad total de visitantes, incluso en presencia de errores en los datos de entrada.

##### ***TP 2***

Dado que el ranking generaba una partición implícita de los museos, se exploró la posibilidad de agruparlos en comunidades, enfocándose en la bisección de la red. Se emplearon dos enfoques: corte mínimo, que busca separar con la menor cantidad de conexiones posibles, y modularidad, que optimiza la cohesión interna de los grupos a través de particiones sucesivas.

En ambos casos, se buscó construir un vector que indique la pertenencia de cada museo a una comunidad. Para el corte mínimo, se utilizó el segundo autovalor más pequeño de la matriz laplaciana, cuyo autovector minimiza el número de enlaces entre comunidades. En el caso de la modularidad, se consideró el mayor autovalor positivo de la matriz $R$, cuyo autovector maximiza la modularidad de la partición.

#### <u>¿Qué otras cosas se pueden hacer con estos datos? </u> 

- En el ***tp1*** habíamos usado que la distancia entre museos era el peso de las aristas de la matríz $D$. Otra alternativa que quizas refleja mejor la probabilidad de pasar de un museo a otro, es usar el tiempo de viaje en transporte publico entre cada punto (teniendo en cuenta la ruta mínima y el tráfico).
- En el ***tp2*** podríamos intentar una clusterización de los museos usando otros métodos como K-Nearest-Neighbors, que vimos en el laboratorio, y compararlo con los resultados que nos dieron.



