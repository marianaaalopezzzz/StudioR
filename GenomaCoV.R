---
title: "Evidencia2.1_A01712807"
author: "Mariana López Colín"
date: "2025-05-04"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
install.packages("ape")

```

```{r}
coronavirus <- c(
  "SARS-CoV-2",      
  "SARS-CoV-1",          
  "MERS-CoV",   
  "HCoV-OC43",    
  "HCoV-HKU1",      
  "HCoV-229E",        
  "RaTG13",      
  "Pangolin-CoV",      
  "HCoV-NL63",         
  "IBV"       
)
```

```{r}
library(rentrez)
library(ape) # Para el análisis filogenético
library(ggplot2) # Para la gráfica
library(dplyr) # Para manipulación de datos

# Lista de variantes de SARS-CoV-2 y otros coronavirus humanos
virus_lista <- c(
  "SARS-CoV-2",      
  "SARS-CoV-1",          
  "MERS-CoV",   
  "HCoV-OC43",    
  "HCoV-HKU1",      
  "HCoV-229E",        
  "RaTG13",      
  "Pangolin-CoV",      
  "HCoV-NL63",         
  "IBV"
)

descargar_y_analizar_secuencias <- function(lista_virus, retmax = 5) {
  resultados <- list()
  secuencias_fasta <- list() # Para almacenar secuencias para el análisis filogenético
  
  for (nombre_virus in lista_virus) {
    cat("Procesando virus:", nombre_virus, "\n")
    
    termino_busqueda <- paste0(nombre_virus, " complete genome")
    cat("Buscando:", termino_busqueda, "\n")
    
    resultado <- entrez_search(db = "nucleotide", term = termino_busqueda, retmax = retmax)
    
    if (length(resultado$ids) == 0) {
      cat("No se encontraron resultados para", nombre_virus, "\n")
      next
    }
    
    fasta <- entrez_fetch(db = "nucleotide", id = resultado$ids[1], rettype = "fasta")
    cat("✅ Secuencia descargada exitosamente\n")
    
    nombre_archivo <- gsub("[^A-Za-z0-9]", "_", nombre_virus)
    archivo_salida <- paste0(nombre_archivo, ".fasta")
    writeLines(fasta, archivo_salida)
    cat("Secuencia guardada en", archivo_salida, "\n")
    
    lineas <- unlist(strsplit(fasta, "\n"))
    secuencia <- paste(lineas[-1], collapse = "")
    bases <- strsplit(secuencia, "")[[1]]
    
    total_bases <- length(bases)
    
    # Cálculo del porcentaje de pares GC
    if (total_bases > 1) {
      pares_consecutivos <- paste0(bases[-length(bases)], bases[-1])
      num_gc_pares <- sum(pares_consecutivos %in% c("GC", "CG")) # Contar GC y CG
      total_pares <- total_bases - 1
      porcentaje_gc_pares <- 100 * num_gc_pares / total_pares
    } else {
      porcentaje_gc_pares <- NA
    }
    
    cat("📏 Longitud de la secuencia:", total_bases, "bases\n")
    cat("🧬 Porcentaje de pares GC:", porcentaje_gc_pares, "%\n")
    
    resultados[[nombre_virus]] <- list(
      fasta = fasta,
      longitud = total_bases,
      gc_porcentaje = porcentaje_gc_pares
    )
    
    secuencias_fasta[[nombre_virus]] <- secuencia # Guardar secuencia para filogenia
  }
  
  return(invisible(list(resultados = resultados, secuencias = secuencias_fasta)))
}

# Ejecutar la función
resultados_completos <- {descargar_y_analizar_secuencias(virus_lista)
}
resultados <- resultados_completos$resultados
secuencias_fasta <- resultados_completos$secuencias
```

```{r grafico_longitudes, echo=FALSE, message=FALSE, warning=FALSE}
longitudes_df <- data.frame(
  Virus = names(resultados),
  Longitud = sapply(resultados, function(x) x$longitud)
)

ggplot(longitudes_df, aes(x = Virus, y = Longitud, fill = Virus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Longitud de Secuencias de Coronavirus",
       x = "Virus",
       y = "Longitud (bases)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  coord_flip()
```

```{r filogenia, echo=FALSE, message=FALSE, warning=FALSE}
hamming_distance <- function(seq1, seq2) {
  if (nchar(seq1) != nchar(seq2)) {
    stop("Las secuencias deben tener la misma longitud para la distancia de Hamming.")
  }
  sum(strsplit(seq1, "")[[1]] != strsplit(seq2, "")[[1]])
}

num_virus <- length(secuencias_fasta)
matriz_distancia_aproximada <- matrix(NA, nrow = num_virus, ncol = num_virus)
rownames(matriz_distancia_aproximada) <- names(secuencias_fasta)
colnames(matriz_distancia_aproximada) <- names(secuencias_fasta)

for (i in 1:num_virus) {
  for (j in 1:num_virus) {
    if (i == j) {
      matriz_distancia_aproximada[i, j] <- 0
    } else {
      seq1 <- secuencias_fasta[[i]]
      seq2 <- secuencias_fasta[[j]]
      min_length <- min(nchar(seq1), nchar(seq2))
      similitud <- sum(strsplit(substr(seq1, 1, min_length), "")[[1]] ==
                         strsplit(substr(seq2, 1, min_length), "")[[1]])
      distancia_aproximada <- 1 - (similitud / min_length)
      matriz_distancia_aproximada[i, j] <- distancia_aproximada
    }
  }
}

distancia_filogenetica_aproximada <- as.dist(matriz_distancia_aproximada)
arbol_aproximado <- ape::nj(distancia_filogenetica_aproximada)
plot(arbol_aproximado, main = "Análisis Filogenético Aproximado de Coronavirus", type = "radial")
```

```{r grafico_final_longitudes, echo=FALSE, message=FALSE, warning=FALSE}
lengths_df <- data.frame(
  virus = names(secuencias_fasta),
  length = sapply(secuencias_fasta, nchar)
)

ggplot(lengths_df, aes(x = reorder(virus, length), y = length, fill = virus)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Longitud del Genoma de Coronavirus", x = "Virus", y = "Número de Bases") +
  theme(legend.position = "none")
```
