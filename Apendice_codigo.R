knitr::opts_chunk$set(echo = TRUE, warning = FALSE)

## 
## install.packages("readxl")
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install("SummarizedExperiment")
## 


# Cargar las librerías
library(readxl)
library(SummarizedExperiment)



file_path <- "GastricCancer_NMR.xlsx"
data_df <- read_excel(file_path, sheet = "Data")
metabolites_df <- read_excel(file_path, sheet = "Peak")



data_df <- data_df[data_df$SampleType != "QC", ]



metabolite_columns <- grep("^M", names(data_df), value = TRUE)
data_matrix <- as.matrix(data_df[, metabolite_columns])


data_matrix <- t(data_matrix)



colData <- DataFrame(
    SampleID = data_df$SampleID,
    SampleType = data_df$SampleType,
    Class = data_df$Class
)

rowData <- DataFrame(
    Name = metabolites_df$Name,
    Label = metabolites_df$Label,
    Perc_missing = metabolites_df$Perc_missing,
    QC_RSD = metabolites_df$QC_RSD
)

se <- SummarizedExperiment(
    assays = list(counts = data_matrix),
    rowData = rowData,
    colData = colData,
    metadata = list(
        description = "Columns M1 ... M149 describe metabolite concentrations. Column SampleType indicates whether the sample was a pooled QC or a study sample. Column Class indicates the clinical outcome observed for that individual: GC = Gastric Cancer, BN = Benign Tumor, HE = Healthy Control."
    )
)

# Filtrar los metabolitos según los criterios de calidad
se_filtered <- se[
  rowData(se)$Perc_missing < 10 & rowData(se)$QC_RSD < 20, 
  ]

print(se_filtered)

## install.packages("ggplot2")

# Cargar las librerías necesarias
library(SummarizedExperiment)
library(ggplot2)

# Extraer la matriz de datos y la información de los grupos
data_matrix <- assay(se_filtered, "counts")
class_labels <- colData(se_filtered)$Class


# Aplicar el Test de Bartlett a cada metabolito
bartlett_results <- apply(data_matrix, 1, function(x) {
  bartlett.test(x ~ class_labels)$p.value
})

# Mostrar los resultados de los primeros 10 metabolitos
bartlett_results_df <- data.frame(
  Metabolite = rownames(data_matrix),
  P_Value = bartlett_results
)
print(bartlett_results_df)


# Aplicar ANOVA a cada metabolito
anova_results <- apply(data_matrix, 1, function(x) {
  summary(aov(x ~ class_labels))[[1]][["Pr(>F)"]][1]
})

# Mostrar los resultados de los primeros 10 metabolitos
anova_results_df <- data.frame(
  Metabolite = rownames(data_matrix),
  P_Value = anova_results
)
print(anova_results_df)


# Función para imputar valores faltantes con la mediana de cada fila (metabolito)
data_matrix <- assay(se_filtered, "counts")
data_matrix <- apply(data_matrix, 1, function(x) {
  x[is.na(x) | is.infinite(x)] <- median(x, na.rm = TRUE)
  return(x)
})
data_matrix <- t(data_matrix) 


# Cargar la librería para el PCA
library(ggplot2)

# Transponer la matriz de datos para que las muestras sean filas y los metabolitos columnas
pca_data <- t(data_matrix)

# Realizar el PCA
pca_result <- prcomp(pca_data, scale. = TRUE)  # scale. = TRUE para normalizar los datos

# Crear un data frame con los resultados del PCA
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Class = class_labels
)

# Gráfico del PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "Análisis de Componentes Principales (PCA)",
    x = "Componente Principal 1",
    y = "Componente Principal 2"
  ) +
  scale_color_manual(values = c("red", "blue", "green"))


## install.packages("dplyr")
## 

library(dplyr)

# Extraer los datos
data_matrix <- assay(se_filtered, "counts")
class_labels <- colData(se_filtered)$Class

# Inicializar un data frame para almacenar los resultados
results <- data.frame(Metabolite = rownames(data_matrix), p_GC_HE = NA, p_GC_BN = NA)

# Realizar las comparaciones entre grupos
for (i in 1:nrow(data_matrix)) {
  # Extraer concentraciones del metabolito actual
  metabolite_data <- data_matrix[i, ]
  
  # Comparación entre GC y HE
  results$p_GC_HE[i] <- wilcox.test(metabolite_data[class_labels == "GC"], 
                                     metabolite_data[class_labels == "HE"])$p.value
  
  # Comparación entre GC y BN
  results$p_GC_BN[i] <- wilcox.test(metabolite_data[class_labels == "GC"], 
                                     metabolite_data[class_labels == "BN"])$p.value
}

# Aplicar corrección por múltiples comparaciones (FDR de Benjamini-Hochberg)
results <- results %>%
  mutate(
    p_adj_GC_HE = p.adjust(p_GC_HE, method = "BH"),
    p_adj_GC_BN = p.adjust(p_GC_BN, method = "BH")
  )

# Mostrar los primeros resultados
print(results)


# Calcular fold change
results <- results %>%
  mutate(
    fold_change_GC_HE = apply(data_matrix, 1, function(x) {
      mean(x[class_labels == "GC"], na.rm = TRUE) / mean(x[class_labels == "HE"], na.rm = TRUE)
    }),
    fold_change_GC_BN = apply(data_matrix, 1, function(x) {
      mean(x[class_labels == "GC"], na.rm = TRUE) / mean(x[class_labels == "BN"], na.rm = TRUE)
    })
  )

# Mostrar los primeros resultados
print(results)


# Instala mixOmics si no lo tienes
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics")

library(mixOmics)


# Extraer la matriz de datos y etiquetas de clase
data_matrix <- t(assay(se_filtered, "counts"))  # Transponer para que las muestras sean filas
class_labels <- colData(se_filtered)$Class

# Convertir las etiquetas de clase a un factor
class_labels <- as.factor(class_labels)


# Ejecutar PLS-DA
plsda_result <- plsda(X = data_matrix, Y = class_labels, ncomp = 2)  # ncomp = 2 componentes

# Graficar el PLS-DA
plotIndiv(plsda_result, group = class_labels, legend = TRUE, 
          title = "PLS-DA de Metabolitos", ellipse = TRUE, comp = c(1, 2))


# Validación cruzada con mixOmics
set.seed(625745221)  # Para reproducibilidad
cv_result <- perf(plsda_result, validation = "Mfold", folds = 5, progressBar = TRUE, nrepeat = 10)

# Ver resultados de validación cruzada
print(cv_result)


# Ver tasa de error
print(cv_result$error.rate)


# Tasa de error por clase
print(cv_result$error.rate.class)


# Graficar resultados de validación cruzada
plot(cv_result)


# Extraer y visualizar las cargas para interpretar los metabolitos que más contribuyen
loadings <- plsda_result$loadings$X  # Cargas para las variables

# Mostrar los primeros metabolitos relevantes
top_metabolites <- rownames(loadings)[order(abs(loadings[, 1]), decreasing = TRUE)[1:10]]
print(top_metabolites)


# Crear un data frame con IDs y Labels de los metabolitos en rowData
metabolite_info <- data.frame(
  ID = rownames(rowData(se_filtered)),
  Label = rowData(se_filtered)$Label
)

# Filtrar para obtener solo los metabolitos en top_metabolites y asegurar el orden correcto
top_metabolite_info <- metabolite_info[match(top_metabolites, metabolite_info$ID), ]

# Mostrar los resultados correctamente alineados
print(top_metabolite_info)


