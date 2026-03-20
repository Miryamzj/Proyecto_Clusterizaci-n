Introducción

El clustering es una estrategia de análisis no supervisado que permite agrupar objetos con base en su similitud o disimilitud. En bioinformática, este tipo de análisis resulta especialmente útil para explorar relaciones entre secuencias biológicas y evaluar si los patrones de agrupamiento recuperan una organización congruente con categorías biológicas conocidas, como la función molecular o la pertenencia taxonómica. En este proyecto, el clustering jerárquico se aplicó a un conjunto de proteínas bacterianas pertenecientes a la familia de las enterotoxinas, utilizando como base una matriz de disimilitud derivada de bit scores obtenidos mediante BLASTP.

Para este análisis se seleccionaron 37 secuencias de proteínas descargadas de UniProtKB, todas correspondientes a enterotoxinas bacterianas. La búsqueda se realizó con el término (reviewed:true) enterotoxin, con el objetivo de trabajar únicamente con secuencias curadas manualmente, lo que reduce la probabilidad de incluir anotaciones erróneas o proteínas mal caracterizadas. Además, se aplicó el filtro taxonómico para el superreino Bacteria, ya que el interés del proyecto fue comparar representantes bacterianos bien definidos. Como criterio adicional, se restringió la longitud de las proteínas a un intervalo de 201 a 400 aminoácidos, con el fin de mantener una comparabilidad estructural razonable entre las secuencias analizadas y evitar que diferencias extremas de tamaño afectaran la interpretación del clustering.

El conjunto final incluyó proteínas de distintos géneros bacterianos y de diferentes tipos de enterotoxinas, lo cual permitió evaluar si los agrupamientos obtenidos reflejaban de alguna manera la taxonomía o el tipo de toxina. Entre los organismos representados se incluyeron Vibrio cholerae, con la subunidad A de la toxina del cólera; Staphylococcus aureus, con diversas enterotoxinas de tipo A, B, C-1, C-2, D, E, G y H; Escherichia coli, con cadenas A de enterotoxinas termolábiles; Streptococcus pyogenes, con exotoxinas pirogénicas de tipo A, C, G y H; y Clostridium perfringens, con la cadena B de una enterotoxina termolábil. Esta diversidad fue útil porque permitió poner a prueba si el clustering era capaz de recuperar agrupamientos biológicamente coherentes a partir de secuencias homólogas, pero taxonómicamente distribuidas.

En un análisis de clustering, una de las primeras consideraciones importantes es precisamente la selección del conjunto de entrada. Si las secuencias no son homólogas, si tienen tamaños demasiado distintos o si el número de elementos es excesivo, el dendrograma puede volverse difícil de interpretar y los agrupamientos pueden reflejar ruido más que relaciones biológicas significativas. Por ello, en este proyecto se procuró trabajar con un conjunto moderado de secuencias, con longitudes comparables y provenientes de proteínas funcionalmente relacionadas.

Otra consideración fundamental es la medida de similitud utilizada. En este trabajo se empleó el bit score derivado de BLASTP, ya que este valor resume la calidad del alineamiento entre dos secuencias y permite comparar pares de proteínas en términos de similitud local. Sin embargo, como los métodos de clustering jerárquico requieren una matriz de distancias o disimilitudes, fue necesario transformar la matriz de similitud en una matriz de disimilitud. Para ello, los bit scores se normalizaron al rango entre 0 y 1, y posteriormente se aplicó la transformación 

$$d_{i,j} = 1 - B_{i,j}$$
	​
donde $B_{i,j}$

 representa la similitud normalizada entre las proteínas $i$ y $𝑗$.

 
También es importante reconocer que el resultado de un clustering jerárquico depende del método de enlace utilizado. Aunque todos los métodos parten de la misma matriz de distancias, cada uno define de manera distinta cómo se calcula la distancia entre clusters, lo que puede producir árboles con estructuras muy diferentes. Por esta razón, en este proyecto se compararon cuatro métodos: single, average, complete y ward. Esta comparación no solo permite observar diferencias visuales en los dendrogramas, sino también discutir cuál de ellos produce agrupamientos más claros, cuáles parecen más congruentes con la taxonomía de las proteínas y cuál captura mejor la estructura jerárquica de los datos.

Finalmente, es importante señalar que un dendrograma obtenido por clustering no equivale automáticamente a un árbol filogenético formal. Aunque ambos representan relaciones entre secuencias, el clustering agrupa con base en similitud y no bajo un modelo evolutivo explícito. Aun así, comparar la estructura de los dendrogramas con la taxonomía esperada de las proteínas es una forma útil de evaluar si el análisis recupera patrones biológicamente razonables.

El objetivo de este proyecto fue implementar en R un flujo de análisis para realizar clustering jerárquico de proteínas enterotóxicas bacterianas a partir de resultados de BLASTP all-vs-all, comparar distintos métodos de agrupamiento y evaluar cuál de ellos representó de manera más informativa la estructura del conjunto de datos.

## Comparación de los árboles obtenidos

Para comparar el efecto del método de agrupamiento sobre la estructura del dendrograma, se generaron árboles mediante cuatro estrategias de clustering jerárquico: **single**, **average**, **complete** y **Ward.D2**.

![Dendrogramas obtenidos con cuatro métodos de clustering jerárquico](..\results\1r_results_prueba_k\figures\dendrogramas_4_metodos.png)


Al comparar los cuatro árboles, se observó que el método **single** fue el menos informativo. Este dendrograma presentó una estructura típica de encadenamiento, en la que muchas secuencias se van agregando progresivamente a un grupo grande sin formar particiones bien delimitadas. Como consecuencia, la interpretación de los clusters resulta menos clara y la comparación con la taxonomía se dificulta.


El método **average** mostró una organización intermedia. Aunque el árbol presentó algo más de estructura que single, todavía se observó un cluster grande que reúne secuencias de distintos grupos, por lo que la separación entre conjuntos biológicamente relacionados no fue completamente clara. El método **complete** generó agrupamientos más compactos que average y single, pero también fragmentó algunas secuencias y dejó varios taxones distribuidos en más de un cluster, lo que redujo su congruencia con la clasificación taxonómica esperada.

El árbol obtenido con **Ward.D2** fue el más informativo. Visualmente, este método produjo la estructura más ordenada y compacta, con tres grupos mejor definidos que en los otros métodos. Además, fue el que mostró una separación más clara entre subconjuntos de proteínas relacionadas, lo que facilitó su interpretación biológica. Aunque el dendrograma de Ward no recuperó de manera perfecta todos los grupos taxonómicos como clusters independientes, sí fue el que más se aproximó a una organización congruente con la procedencia de varias de las proteínas analizadas.

En conjunto, los cuatro métodos mostraron que existe cierta estructura en los datos, pero no todos recuperan esa estructura con la misma claridad. De los árboles obtenidos, el **más informativo fue el de Ward.D2**, mientras que el **menos informativo fue el de single**.

En cuanto a la congruencia con la taxonomía, la inspección visual de las etiquetas mostró que **ninguno de los cuatro árboles fue completamente congruente con la taxonomía de todas las proteínas**, ya que varios grupos quedaron mezclados o repartidos entre distintos clusters dependiendo del método. Sin embargo, **Ward.D2 fue el árbol que mostró la mejor aproximación a la taxonomía**, mientras que **single** fue el que mostró menor congruencia.

## Agglomerative Coefficient

Para complementar la comparación visual de los dendrogramas, se calculó el **Agglomerative Coefficient** para cada método. Este coeficiente permite evaluar qué tan fuerte es la estructura jerárquica detectada por el algoritmo; valores más altos indican agrupamientos más definidos.

Los valores obtenidos fueron:

- **single**: 0.6688
- **average**: 0.6745
- **complete**: 0.6745
- **ward**: 0.8409

De los cuatro métodos evaluados, el árbol con el **Agglomerative Coefficient más alto** fue el de **ward**, con un valor de **0.8409**. Este resultado coincide con la inspección visual de los dendrogramas, ya que Ward también fue el método que produjo la estructura más clara y compacta. En cambio, el método **single** presentó el valor más bajo, lo que refuerza la idea de que fue el menos informativo para este conjunto de proteínas.

## Conclusión

En este proyecto se implementó un flujo de análisis en **R** para realizar clustering jerárquico de proteínas bacterianas enterotóxicas a partir de resultados de **BLASTP all-vs-all**. Como parte del proceso, se seleccionaron y curaron previamente las secuencias, se construyó una matriz de similitud con los bit scores, se normalizó y se transformó en una matriz de disimilitud para aplicar distintos métodos de agrupamiento.

Se compararon los métodos **single**, **average**, **complete** y **Ward**, tanto de forma visual como mediante el **Agglomerative Coefficient**. En conjunto, los resultados mostraron que **Ward** fue el método más informativo, ya que produjo agrupamientos más claros y presentó el coeficiente más alto, mientras que **single** fue el menos informativo. Aunque ninguno de los árboles fue completamente congruente con la taxonomía, el proyecto permitió demostrar que la elección del método de clustering influye de manera importante en la interpretación biológica de los datos.
