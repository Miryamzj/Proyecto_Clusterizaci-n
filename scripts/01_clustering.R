# ============================================================
# Proyecto Módulo Clustering
# Clustering jerárquico a partir de BLASTP all-vs-all 
# ============================================================

# Librerías 
library(cluster) 
library(ape)
library(factoextra)
library(dendextend)
library(corrplot)

# Rutas
blast_path <- "results/blast/resultados_blast.txt"

# Crear carpetas de salida
dir.create("results/matrices", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/trees",    recursive = TRUE, showWarnings = FALSE)
dir.create("results/clusters", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Leer salida BLAST
# ============================================================

blast <- read.table(blast_path, header = FALSE, stringsAsFactors = FALSE)

# outfmt 6 típico tiene 12 columnas:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
colnames(blast) <- c("query","subject","pident","aln_len","mismatch","gapopen",
                     "q_start","q_end","s_start","s_end","evalue","bit_score")

# Revisamos rápidamente
cat("Filas (hits):", nrow(blast), "\n")
cat("Columnas:", ncol(blast), "\n")

# ============================================================
# Matriz de similitud (bit score)
# ============================================================

# IDs únicos de proteínas
unique_seqs <- unique(c(blast$query, blast$subject))
cat("Proteínas únicas:", length(unique_seqs), "\n")

# Matriz cuadrada inicializada en 0
sim <- matrix(0,
              nrow = length(unique_seqs),
              ncol = length(unique_seqs),
              dimnames = list(unique_seqs, unique_seqs))

# Llenar la matriz con bit scores
# (si hay más de un hit para el mismo par, nos quedamos con el más alto)
for (i in 1:nrow(blast)) {
  q <- blast$query[i]
  s <- blast$subject[i]
  b <- blast$bit_score[i]
  
  if (b > sim[q, s]) sim[q, s] <- b
  if (b > sim[s, q]) sim[s, q] <- b
}

# Guardar matriz de similitud sin normalizar
write.table(sim, "results/matrices/similitud_bitscore.tsv",
            sep = "\t", quote = FALSE)

# ============================================================
# Normalización [0,1] y disimilitud
# ============================================================

max_value <- max(sim)
cat("Máximo bitscore:", max_value, "\n")

sim_norm <- sim / max_value
diag(sim_norm) <- 1  # similitud de i con i es 1

diss <- 1 - sim_norm

# Guardar matrices
write.table(sim_norm, "results/matrices/similitud_normalizada.tsv",
            sep = "\t", quote = FALSE)

write.table(diss, "results/matrices/disimilitud.tsv",
            sep = "\t", quote = FALSE)

# Convertir a objeto dist para clustering
d <- as.dist(diss)

# ============================================================
# Clustering con 4 métodos 
# ============================================================

csin <- hclust(d, method = "single")
cave <- hclust(d, method = "average")
ccom <- hclust(d, method = "complete")
cwar <- hclust(d, method = "ward.D2")

# k solo para visualizar y cortar árboles
k <- 3

# Guardar dendrogramas + rectángulos 
png("results/figures/dendrogramas_4_metodos.png", width = 2200, height = 900)
par(mfrow = c(2,2), mar = c(8,3,3,1))

plot(csin, hang = -1, main = "Single", xlab = "", sub = "", cex = 0.65)
rect.hclust(csin, k = k, border = 1:16)

plot(cave, hang = -1, main = "Average", xlab = "", sub = "", cex = 0.65)
rect.hclust(cave, k = k, border = 1:16)

plot(ccom, hang = -1, main = "Complete", xlab = "", sub = "", cex = 0.65)
rect.hclust(ccom, k = k, border = 1:16)

plot(cwar, hang = -1, main = "Ward.D2", xlab = "", sub = "", cex = 0.65)
rect.hclust(cwar, k = k, border = 1:16)

dev.off()

# Cortar árboles (guardar clusters) (como el ejemplo)
cl_single   <- cutree(csin, k = k)
cl_average  <- cutree(cave, k = k)
cl_complete <- cutree(ccom, k = k)
cl_ward     <- cutree(cwar, k = k)

write.table(data.frame(id = names(cl_single),   cluster = cl_single),
            "results/clusters/clusters_single_k4.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(id = names(cl_average),  cluster = cl_average),
            "results/clusters/clusters_average_k4.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(id = names(cl_complete), cluster = cl_complete),
            "results/clusters/clusters_complete_k4.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(id = names(cl_ward),     cluster = cl_ward),
            "results/clusters/clusters_ward_k4.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================
# Exportar árboles en Newick
# ============================================================

write.tree(as.phylo(csin), file = "results/trees/dendrogram_single.nwk")
write.tree(as.phylo(cave), file = "results/trees/dendrogram_average.nwk")
write.tree(as.phylo(ccom), file = "results/trees/dendrogram_complete.nwk")
write.tree(as.phylo(cwar), file = "results/trees/dendrogram_ward.nwk")

# ============================================================
# Correlación entre dendrogramas (Baker) + corrplot
# ============================================================

trees <- dendlist(
  Ward     = as.dendrogram(cwar),
  Single   = as.dendrogram(csin),
  Average  = as.dendrogram(cave),
  Complete = as.dendrogram(ccom)
)

baker <- cor.dendlist(trees, method = "baker")

png("results/figures/correlacion_baker.png", width = 900, height = 800)
corrplot(baker, method = "circle", type = "lower",
         tl.col = "black", addCoef.col = "white",
         number.cex = 1.0, col.lim = c(0, 1))
dev.off()

# ============================================================
# Visualización de clusters 
# ============================================================

png("results/figures/cluster_plot_complete_k4.png", width = 900, height = 700)
print(fviz_cluster(list(data = as.data.frame(diss), cluster = cl_complete),
                   ellipse.type = "convex", geom = "point"))
dev.off()

png("results/figures/cluster_plot_ward_k4.png", width = 900, height = 700)
print(fviz_cluster(list(data = as.data.frame(diss), cluster = cl_ward),
                   ellipse.type = "convex", geom = "point"))
dev.off()

# ============================================================
# Estimar k óptimo 
# ============================================================

png("results/figures/k_optimo_wss.png", width = 900, height = 700)
print(fviz_nbclust(as.data.frame(diss), FUN = hcut,
                   hc_method = "ward.D2", method = "wss", k.max = 10))
dev.off()

png("results/figures/k_optimo_silhouette.png", width = 900, height = 700)
print(fviz_nbclust(as.data.frame(diss), FUN = hcut,
                   hc_method = "ward.D2", method = "silhouette", k.max = 10))
dev.off()

png("results/figures/k_optimo_gap.png", width = 900, height = 700)
print(fviz_nbclust(as.data.frame(diss), FUN = hcut,
                   hc_method = "ward.D2", method = "gap_stat", k.max = 10))
dev.off()

# ============================================================
# Agnes 
# ============================================================

ac_single   <- agnes(d, method = "single")$ac
ac_average  <- agnes(d, method = "average")$ac
ac_complete <- agnes(d, method = "complete")$ac
ac_ward     <- agnes(d, method = "ward")$ac

ac_tab <- data.frame(
  method = c("single","average","complete","ward"),
  agglomerative_coefficient = c(ac_single, ac_average, ac_complete, ac_ward)
)

print(ac_tab)

write.table(ac_tab, "results/agglomerative_coefficients.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nListo. Se generaron:\n",
    "- results/matrices/\n",
    "- results/figures/ (dendrogramas, baker, k óptimo, cluster plots)\n",
    "- results/trees/*.nwk\n",
    "- results/clusters/\n",
    "- results/agglomerative_coefficients.tsv\n")