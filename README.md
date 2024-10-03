# README MARRONERO

El programa de forma estandar trabaja asumiendo que habrá dos variables del modelo que queremos *simular*: una variable continua y la otra discreta. 

Esto se puede cambiar *hardcodeando* y poniendo explicitamene las variables dentro del array ``VARS``. Para añadir variables continuas, añadir en VARS strings con el formato ``c{q_size}``. En caso de variables discretas, añadir ``d{q_size}``, siendo ``q_size`` el tamaño en bits de cada variable.

## Como ejecutar
Esto es easy. Todo el programa está contenido en el fichero **main.py**, por lo que para ejecutar solo es necesario llamar al interprete: ``python main.py C_Q_SIZE D_Q_SIZE R_INDEX OUT_FILE``.

En caso de ejecutar en el cluster: ``bnd -exec python main.py C_Q_SIZE D_Q_SIZE R_INDEX OUT_FILE``.

### Parámetros de las ejecuciones
El programa pide **4 parámetros** al usuario:

1. **C_Q_SIZE**:  Tamaño en qbits de la variable continua. Rango entre [1, inf)
2. **D_Q_SIZE**: Tamaño en qbits de la variable discreta. Rango entre [0, inf)
3. **R_INDEX**: Indica el resultado final que buscamos de la matriz. Rango entre [0, X] siendo ``X=(2^C_Q_SIZE + 2^D_Q_SIZE)``. El valor ``0`` indica que queremos que el circuito devuelva la equiprobabilidad entre todos los estados. Un valor ``i`` en el rango [1,X], indica que buscamos que el circuito converga a una solución que tenga toda la probabilidad asignada al estado ``i-1``.

    >**EXCEPCIÓN**: Si ``D_Q_SIZE=0``, entonces ``X=2^C_Q_SIZE``.

4. **OUT_FILE**: Nombre del fichero de salida. No está implementado pero el programa lo pide igualmente. Para obtener los ``print``s desviar la salida estandar a fichero: ``python main.py ... > out_file``.

## Notas:

El minimizador exige que las variables vengan descritas en un array unidimensional de floats. Para poder representar una matriz y de numeros complejos creo un array unidimensional de floats de 2*(m_size, m_size) siendo m_size el tamaño de la matriz que buscamos. De esta forma, la primera y segunda  mitad del array representa las partes reales e imaginarias de los valores de la matriz.

El pasar de float a complejo y viceversa se hace dentro de las funciones que se pasan al optimizador.

---

He ido tocando y ahora el minimizador casca por algo de los tipos. He estado tocando el tema de los números complejos dentro de la matriz y ahora peta.

