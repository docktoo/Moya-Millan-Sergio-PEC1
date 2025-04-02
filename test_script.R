# Cargar las librerías necesarias
library(SummarizedExperiment)
library(readr)
library(ggplot2)
library(e1071)  # Para calcular la asimetría
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap) 
# Leer el archivo CSV
data <- read_csv("~/Moya-Millan-Sergio-PEC1/human_cachexia.csv")

# Separar los metadatos y la matriz de datos de metabolitos
metadata <- data %>% select(`Patient ID`, `Muscle loss`)
metabolites <- data %>% select(-`Patient ID`, -`Muscle loss`)

# Asignar nombres de filas para que coincidan en ambos objetos
rownames(metadata) <- metadata$`Patient ID`
rownames(metabolites) <- metadata$`Patient ID`

# Crear los gráficos de la distribución de los metabolitos 
metabolites_long <- metabolites %>%
  gather(key = "Metabolite", value = "Value")  

# Gráfico de boxplot 
ggplot(metabolites_long, aes(x = Metabolite, y = Value)) +
  geom_boxplot(fill = "lightblue", color = "black", alpha = 0.7) +
  xlab("Metabolitos") +
  ylab("Valor del metabolito") +
  ggtitle("Distribución de Metabolitos") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  

# Aplicar la transformación logarítmica a todos los metabolitos 
metabolites_log <- log2(metabolites)

# Crear los gráficos de la distribución de los metabolitos después de la transformación logarítmica
metabolites_log_long <- as.data.frame(metabolites_log) %>%
  gather(key = "Metabolite", value = "Value")  

# Gráfico de boxplot después de la transformación
ggplot(metabolites_log_long, aes(x = Metabolite, y = Value)) +
  geom_boxplot(fill = "lightgreen", color = "black", alpha = 0.7) +
  xlab("Metabolitos") +
  ylab("Valor transformado") +
  ggtitle("Distribución de Metabolitos después de log-transformación") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  

# Realizar el escalado automático de los metabolitos
# Estándar: centrado y escalado (promedio = 0, desviación estándar = 1)
metabolites_scaled <- scale(metabolites_log)

# Transponer la matriz de metabolitos para que las columnas sean las muestras
metabolites_t_scaled <- t(metabolites_scaled)

# Crear el objeto SummarizedExperiment
ses <- SummarizedExperiment(
  assays = list(counts = as.matrix(metabolites_t_scaled)),  # Convertir metabolitos escalados en matriz
  colData = DataFrame(metadata)  # Convertir metadatos en DataFrame
)

# Mostrar el objeto SummarizedExperiment
ses
str(ses)

# Extraer los datos de los ensayos y la información de los grupos
counts <- assays(ses)$counts
group <- colData(ses)$Muscle.loss

# Crear un data frame para almacenar los resultados de las pruebas t
results <- data.frame(
  Feature = rownames(counts),
  t_statistic = numeric(nrow(counts)),
  p_value = numeric(nrow(counts))
)

# Realizar pruebas t para cada característica
for (i in 1:nrow(counts)) {
  feature_data <- counts[i, ]
  t_test <- t.test(feature_data ~ group)
  results$t_statistic[i] <- t_test$statistic
  results$p_value[i] <- t_test$p.value
}

# Ajustar los valores p para el control de la tasa de falsos descubrimientos (FDR)
results$adjusted_p_value <- p.adjust(results$p_value, method = "fdr")

# Mostrar los resultados
print(results)
# Filtrar los 10 metabolitos más significativos por valor p ajustado
top_10_metabolites <- results[order(results$adjusted_p_value), ][1:10, ]

# Ver los primeros 10 metabolitos
print(top_10_metabolites)


# Extraer los datos de los metabolitos (matriz de expresión)
metabolite_data <- assay(se)

# Transponer la matriz de metabolitos para que las filas sean metabolitos y las columnas muestras
metabolite_data_transposed <- t(metabolite_data)

# Calcular la matriz de correlación entre las filas (ahora metabolitos)
cor_matrix <- cor(metabolite_data_transposed, method = "pearson")

# Generar el heatmap para las correlaciones entre los metabolitos
pheatmap(cor_matrix,
         cluster_rows = TRUE,      # Agrupar metabolitos (filas)
         cluster_cols = FALSE,     # No agrupar muestras (columnas)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Paleta de colores
         main = "Mapa de calor de correlaciones entre metabolitos",
         display_numbers = FALSE,  # No mostrar los números de las correlaciones
         fontsize = 10,            # Tamaño de la fuente
         scale = "none")           # No escalar los 

summary(cor_matrix)
# Convertir la matriz de correlación en un data frame
cor_df <- as.data.frame(as.table(cor_matrix))

# Filtrar la diagonal principal y duplicados
cor_df <- cor_df[cor_df$Var1 != cor_df$Var2, ]

# Ordenar por valores absolutos de correlación (de mayor a menor)
cor_df <- cor_df[order(abs(cor_df$Freq), decreasing = TRUE), ]

# Ver los pares de metabolitos más correlacionados
head(cor_df, 10)  # Muestra los 10 pares más correlacionados

# Calcular la matriz de distancias euclidianas entre los metabolitos
distance_matrix <- dist(metabolite_data_transposed, method = "euclidean")

# Realizar el agrupamiento jerárquico usando el método de enlace completo
hc <- hclust(distance_matrix, method = "complete")

# Visualizar el dendrograma (agrupamiento jerárquico)
plot(hc, main = "Dendrograma de Agrupamiento Jerárquico", 
     xlab = "Metabolitos", ylab = "Distancia Euclidiana", 
     cex = 0.7)  # Ajusta el tamaño de los labels si es necesario

# Calcular la suma de los cuadrados dentro de los clusters (WSS) para diferentes valores de K
set.seed(42)
wss <- sapply(1:10, function(k) kmeans(metabolite_data_transposed, centers = k, nstart = 10)$tot.withinss)

wss
# Graficar el método del codo
plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Número de Clusters K", ylab = "Suma de los cuadrados dentro del cluster",
     main = "Método del Codo")

# Realizar PCA sobre la matriz de metabolitos
pca <- prcomp(t(counts))  # Transponer 'counts' para que las muestras sean las observaciones

# Crear un data frame con las coordenadas PCA y las etiquetas de los clusters
pca_data <- data.frame(pca$x)
pca_data$Cluster <- factor(clusters)  # Asignar los clusters como factor

# Graficar el PCA
library(ggplot2)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "PCA con 4 Clusters", x = "Componente Principal 1", y = "Componente Principal 2") +
  theme_minimal()
