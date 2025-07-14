# Ejercicio 1
### Inciso a
Probar que $\mathbb{1}$, el vector de $n$ unos, es autovector de $R$ y $L$, dar su autovalor, y explicar qué agrupación de la red representa ese autovector.

Veamos primero que es autovector de $L$.
Veamos que 
$$L_{ij} = \begin{cases}
        \sum_{h=1}^{n}{A_{ih}} - A_{ij} & \text{si } i=j \\ 
        -A_{ij} & \text{si } i\neq j 
\end{cases}$$

Reducimos esta expresion, usando que $A_{ii} = 0$, para cualquier $i$.

$$L_{ij} = \begin{cases}
        \sum_{h=1}^{n}{A_{ih}} & \text{si } i=j \\ 
        -A_{ij} & \text{si } i\neq j 
\end{cases}$$

Veamos que cumple la ecuacion

$$L \cdot \mathbb{1} = \lambda_1 \cdot \mathbb{1}$$

Esto es lo mismo que ver un sistema de ecuaciones de $n$ ecuaciones. Miremos una arbitraria, con el indice $k$

$$\sum_{i=1}^{n} L_{ki} = \lambda_1$$

$$\sum_{i=1 \land i\neq k}^{n} L_{ki} + L_{kk} = \lambda_1$$

$$\sum_{i=1 \land i\neq k}^{n}{-A_{ki}} + \sum_{h=1}^{n}{A_{kh}} = \lambda_1$$

$$\sum_{h=1}^{n}{A_{kh}} = \lambda_1 + \sum_{i=1 \land i\neq k}^{n}{A_{ki}}$$

Si hacemos un cambio de indices, vemos que de un lado nos queda solo el $kk$ y del otro $\lambda_1$

$$\sum_{h=1}^{n}{A_{kh}} = \lambda_1 + \sum_{h=1 \land h\neq k}^{n}{A_{ki}}$$

$$A_{kk} = \lambda_1$$

Luego, vemos que, como esta era una ecuacion genérica, esto vale para cualquier $k$. Como sabemos que la diagonal de la matriz de adyacencias $A$ es 0 (pues nunca un museo se conecta consigo mismo), concluimos que:

$$\lambda_1 = 0$$


Veamos lo mismo para la matriz $R = A - P$, donde $P_{ij} = \frac{k_i \cdot k_j}{2E}$. 
En esta ecuacion, E es la cantidad total de aristas, es decir, $2E = \sum_{i=1}^{n}{\sum_{j=1}^{n}{A_{ij}}}$.
Ademas, $k_i = K_{ii} = \sum_{j=1}^{n}{A_{ij}}$

Veamos entonces la demostracion de que $\mathbb{1}$ es autovector de $R$, y con qué autovalor. Lo llamamos $\lambda_2$.

$$R \cdot\mathbb{1} = \lambda_2 \cdot \mathbb{1}$$

Vamos a usar la misma idea de la demostración anterior. Podemos escribir esto como un sistema de $n$ ecuaciones, y tomamos una genérica, la $k$-ésima

$$
\sum_{j=1}^{n}{R_{kj}} = \lambda_2
$$

$$
\sum_{j=1}^{n}{A_{kj} - \frac{k_k \cdot k_j}{2E}} = \lambda_2
$$

$$
\sum_{j=1}^{n}{A_{kj} = \lambda_2 + \sum_{j=1}^{n}{\frac{k_k \cdot k_j}{2E}}} 
$$

$$
\sum_{j=1}^{n}{A_{kj} = \lambda_2 + k_k \cdot \sum_{j=1}^{n}\frac{k_j}{2E}} 
$$

Recordamos que $k_k = \sum_{h=1}^{n}{A_{kh}}$. Luego:

$$
\sum_{j=1}^{n}{A_{kj} = \lambda_2 + (\sum_{h=1}^{n}{A_{kh}}) \cdot \sum_{j=1}^{n}{\frac{\sum_{h=1}^{n}{A_{jh}}}
{2E}}} 
$$

$$
\sum_{j=1}^{n}{A_{kj} = \lambda_2 + \frac{\sum_{h=1}^{n}{A_{kh}}}{2E}\sum_{j=1}^{n}({\sum_{h=1}^{n}{A_{jh}}})} 
$$

$$
\sum_{j=1}^{n}{A_{kj} = \lambda_2 + \frac{\sum_{h=1}^{n}{A_{kh}}}{2E} \sum_{j=1}^{n}\sum_{h=1}^{n}{A_{jh}}}
$$

Recordemos que $2E = \sum_{i=1}^{n}{\sum_{j=1}^{n}{A_{ij}}}$. Luego:

$$
\sum_{j=1}^{n}{A_{kj} = \lambda_2 + \sum_{h=1}^{n}{A_{kh}}}
$$

Si renombramos el indice $h$ al indice $j$, vemos que nos queda:

$$\lambda_2 = 0$$

Luego, queda demostrado que el vector $\mathbb{1}$ es autovector de las matrices $R$ y $L$, con el mismo autovalor $\lambda_1 = \lambda_2 = 0$.

La interpretacion del autovector $\mathbb{1}$ es la asignación de todos los museos al mismo grupo, ya que la asignación a un grupo o a otro es con $1$ o $-1$

### Inciso b
Usamos la sugerencia del enunciado, primero demostrando que, para una matriz simétrica arbitraria $M$, con dos autovectores $v1$, $v2$ con $\lambda_1 \neq \lambda_2$ autovalores asociados. Veamos que ambas ecuaciones, $v_1^t M v_2$ y $v_2^t M v_1$ son un escalar, por lo tanto su transpuesta es la misma. Esto lo podemos verificar asi, sabiendo que M es simétrica.

$$v_1^t M v_2 = v_2^t M^t v_1 = v_2^t M v_1$$ 

Luego, ambas expresiones son el mismo escalar.

$$v_1^t M v_2 = v_2^t M v_1$$

$$v_1^t \lambda_2 v_2 = v_2^t \lambda_1 v_1$$

$$\lambda_2 v_1^t v_2 = \lambda_1 v_2^t  v_1$$

$$\lambda_2 v_1^t v_2 - \lambda_1 v_2^t  v_1 = 0$$

Usamos que $v_1^t v_2 = v_2^t v_1$

$$(\lambda_2 - \lambda_1)(v_1^t v_2) = 0$$

Luego, alguno de los dos terminos tiene que ser cero. 

Como ya dijimos que los autovalores son distintos, no queda otra que concluir.

$$v_1^t v_2 = 0$$

Luego, tenemos que verificar que esto vale para las matrices $L$ y $R$. 
Quiero ver que $L^t = L$

$$L^t = K^t - A^t = A$$

Esto vale trivialmente puesto que $A$ es simétrica y $K$ es diagonal.

Quiero ver que $R^t = R$

$$R^t = A^t - P^t$$

Como ya vimos $A$ es simétrica, queda ver que lo sea $P$. Queremos ver que:

$$P_{ij} = P_{ji}$$

$$\frac{k_i \cdot k_j}{2E} = \frac{k_j \cdot k_i}{2E} $$

Lo cual es trivialmente cierto. Luego, queda demostrado que la propiedad vale para $R$ y $L$.



### Inciso c
Para este ejercicio vamos a usar la propiedad recién demostrada. Tomamos una matriz simétrica $M$ que puede ser $R$ o $L$. También usamos que $\mathbb{1}$ es autovector de $R$ y de $L$ con autovalor 0. Como sabemos que el autovalor de $v \neq 0$, podemos usar la propiedad de que el producto de dos vectores con autovalores distintos de una matriz simétrica es igual a 0. 

$$v_1^t v_2 = 0$$


En este caso particular, la ecuación queda así.

$$v \cdot \mathbb{1} = 0$$

$$v_1 + v_2 + \dots + v_n = 0$$

$$\sum_{i=1}^{n}{v_i} = 0$$

Tal como se quería demostrar.

# Ejercicio 2
### Inciso a

Mostramos que
1) autovalores de $M + \mu I$ son $\gamma_i = \lambda_i + \mu$
2) El autovector asociado a $\gamma_i$ es $v_i$
3) Si $\mu + \lambda_i \neq 0 \forall i$, entonces $M + \mu I$ es inversible

Sabemos que, para $v_i \neq 0$, por ser autovector, vale que

$$M v_i = \lambda_i v_i$$

Luego, desarollemos la siguiente expresion

$$(M + \mu I) v_i = \lambda_i v_i + \mu v_i = (\lambda_i + \mu) v_i$$

Luego, si tomamos $\gamma_i = \lambda_i + \mu$, vemos que

$$(M + \mu I) v_i = \gamma_i v_i$$

Luego, $v_i$ es autovector de $M + \mu I$ con autovalor $\gamma_i$. Quedan demostrados los puntos 1 y 2.

Veamos el punto 3. Tomamos $i$ arbitraria. Queremos ver que

$$\mu + \lambda_i \neq 0 \implies M + \mu I \text{ es inversible}$$

$$\gamma_i \neq 0 \implies M + \mu I \text{ es inversible}$$

Luego, esto vale por la propiedad de que una matriz es inversible si y solo si todos sus autovalores son distintos de 0. Luego, si $\gamma_i \neq 0 \forall i$, entonces $M + \mu I$ es inversible.

### Inciso b

#### Parte 1
Para $\mu \gt 0$, mostrar que $L + \mu I$ es inversible, sabiendo que los autovalores de $L$ son no negativos. Recordamos tambien que $L = K - A$, con $K_{ii} = \sum_{j=1}^{n}{A_{ij}}$

Usando el inciso anterior, vemos que podemos definir los autovalores de $L + \mu I$ como $\gamma_i = \lambda_i + \mu$, siendo $\lambda_i$ el i-esimo autovalor de $L$ (ordenados descendientemente).

Luego, sabemos que $\mu \gt 0$ y $\lambda_i \geq 0$. Luego, $\gamma_i \gt 0$. Entonces, como $L + \mu I$ tiene todos sus autovalores positivos, es inversible.

#### Parte 2
Mostrar que aplicar el metodo de la potencia a $(L + \mu I)^{-1}$ converge al autovector de $L$ asociado a su autovalor mas chico si se parte de una semilla adecuada. 

Veamos que, en general, para una matriz $A$ inversible y diagonalizable, vale que $A v_i = \lambda_i v_i$. Luego, vale que $A^{-1} v_i = \frac{1}{\lambda_i} v_i$.

Ya sabemos que $L + \mu I$ es inversible, por el inciso anterior. Para ver que es diagonalizable, usamos el teorema espectral, que nos dice que por ser simétrica real, es diagonalizable. Esto lo vemos usando que $L$ es simetrica real. Luego, $L + \mu I$ tambien lo es. Luego, $L + \mu I$ es diagonalizable. Finalmente, vemos que los autovectores de $L + \mu I$ son los mismos que los de $L$, con autovalores asociados iguales a los de $L$ mas $\mu$. Entonces, por lo que vimos recién, los autovectores de $(L + \mu I)^{-1}$ son los mismos que los de $L$, con autovalores $\frac{1}{\lambda_i + \mu}$.

Veamos los autovalores de $L$ con un orden particular. Supongamos que $\lambda_1 \geq \lambda_2 \geq \dots \geq \lambda_n \geq 0$.

Tomemos $A = (L + \mu I)^{-1}$.

Sabiendo que, para una matriz $A$, el metodo de la potencia converge al autovector de $A$ asociado al mayor autovalor (son todos positivos asi que no importa el modulo aca), vemos que, si partimos de una semilla adecuada, converge al autovector de $A$ asociado al mayor autovalor, que es $\frac{1}{\lambda_n + \mu}$, ya que para maximizar esa expresion, minimizamos el denominador. Y el minimo es $\lambda_n$, que justamente es el autovalor mas chico de $L$.

Las condiciones de "semilla adecuada" refieren a que el vector inicial para el metodo de la potencia no sea ortogonal a $v_n$.

El caso de que haya solamente un autovector asociado al autovalor mas chico sucede cuando la multiplicidad de ese autovalor es 1. 

### Inciso c
Tenemos una matriz $M \in \mathbb{R}^{n \times n}$ simétrica, y queremos ver que $M - \lambda_1 \frac{v_1 v_1^t}{v_1^t v_1}$ tiene los mismos autovectores que $M$, pero el autovalor asociado a $v_1$ es 0.

Como $M$ es simétrica real, admite una b.o.n de autovectores. Por lo tanto, vemos que la expresión $v_1^t v_1$ es lo mismo que $||v_1||_2^2 = 1$, pues $v_1$ pertenece a una b.o.n. Luego, escribimos directamente la expresión $M - \lambda_1 v_1 v_1^t$.

Veamos que se cumple lo que pide el ejercicio

$$(M - \lambda_1 v_1 v_1^t)  v_i = \gamma_i v_i$$

Separamos en casos. Primero vemos $i=0$ y luego $i>0$

$$(M - \lambda_1 v_1 v_1^t)  v_1 = \gamma_1 v_1$$

$$M  v_1 - \lambda_1 v_1 v_1^t  v_1 = \gamma_1 v_1$$

$$ \lambda_1 v_1 - \lambda_1 v_1 = \gamma_1 v_1$$

$$0 = \gamma_1 v_1$$

Luego, vemos que $\gamma_1 = 0$, pues $v_1 \neq 0$ por ser autovector. Veamos $i>0$

$$(M - \lambda_1 v_1 v_1^t) v_i = \gamma_i v_i$$

$$M v_i - \lambda_1 v_1 v_1^t v_i = \gamma_i v_i$$

$$\lambda_i v_i - \lambda_1 v_1 (v_1^t v_i) = \gamma_i v_i$$

Vemos que $v_1^t v_i = 0$, pues $v_1$ y $v_i$ son autovectores distintos, que pertenecen a una b.o.n.

Luego, queda que:

$$\lambda_i v_i = \gamma_i v_i$$

$$\gamma_i = \lambda_i$$

Entonces, vale que los autovectores de $M$ son los mismos que los de $M - \lambda_1 v_1 v_1^t$, con autovalores asociados iguales a los de $M$ menos $\lambda_1$, que es lo que queríamos demostrar.
