# Proyecto_Clusterizaci-n

# Proyecto Módulo Clustering

Proyecto para el módulo de **clustering jerárquico** de la Licenciatura en Ciencias Genómicas.  
El objetivo fue construir un flujo de análisis en **R** para agrupar proteínas homólogas a partir de una matriz de disimilitud derivada de resultados de **BLASTP all-vs-all**, y comparar distintos métodos de agrupamiento jerárquico.

## Objetivo

Desarrollar un programa en R que:

- lea los resultados de una comparación all-vs-all entre proteínas homólogas,
- construya una matriz de similitud con base en **bit scores**,
- convierta dicha matriz en una matriz de **disimilitud**,
- aplique **clustering jerárquico** con distintos métodos,
- exporte los dendrogramas en formato **Newick**,
- compare los árboles generados,
- y calcule el **Agglomerative Coefficient** para cada método.

## Preparación previa de los datos

Antes de ejecutar el script principal, se realizó una preparación manual del conjunto de secuencias:

- selección de un conjunto de proteínas homólogas en formato FASTA,
- uso de no más de 100 secuencias para facilitar la interpretación,
- verificación de longitudes comparables entre proteínas,
- modificación de los identificadores para incluir información taxonómica,
- ejecución de **BLASTP all-vs-all** para obtener los bit scores entre cada par de secuencias.

El archivo FASTA utilizado fue:

- `enterotoxinas.fasta.gz`

La salida de BLAST se guardó en:

- `results/blast/resultados_blast.txt`

## Estructura del repositorio

```text
Proyecto_Clusterizaci-n/
├── results/
│   ├── blast/
│   │   ├── db/
│   │   └── resultados_blast.txt
│   ├── clusters/
│   │   ├── clusters_average_k3.tsv
│   │   ├── clusters_complete_k3.tsv
│   │   ├── clusters_single_k3.tsv
│   │   └── clusters_ward_k3.tsv
│   ├── figures/
│   │   ├── cluster_plot_complete_k3.png
│   │   ├── cluster_plot_ward_k3.png
│   │   ├── correlacion_baker.png
│   │   ├── dendrogramas_4_metodos.png
│   │   ├── k_optimo_gap.png
│   │   ├── k_optimo_silhouette.png
│   │   └── k_optimo_wss.png
│   ├── matrices/
│   │   ├── disimilitud.tsv
│   │   ├── similitud_bitscore.tsv
│   │   └── similitud_normalizada.tsv
│   ├── trees/
│   │   ├── dendrogram_average.nwk
│   │   ├── dendrogram_complete.nwk
│   │   ├── dendrogram_single.nwk
│   │   └── dendrogram_ward.nwk
│   └── agglomerative_coefficients.tsv
├── scripts/
│   └── 01_clustering.R
├── .gitignore
├── LICENSE
├── Proyecto_Clusterizaci-n.Rproj
├── ProyectoModuloClustering.pdf
├── README.md
└── enterotoxinas.fasta.gz