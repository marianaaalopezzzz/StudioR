---
title: "Evidencia 1 | Análisis inicial"
subtitle: "Mariana López Colín A01712807"
output: html_document

---

```{r setup, include=FALSE}

```{r}
variantes <- c(
  "B.1.1.7",      
  "B.1.351",          
  "P.1",   
  "B.1.617.2",    
  "BA.2",      
  "BA.5",        
  "XBB",      
  "BA.2.86",      
  "P.2",         
  "B.1.526"       
)
library(rentrez)

descargar_y_analizar_varias_secuencias <- function(lista_variantes, retmax = 5) {
  
  resultados <- list()
  
  for (nombre_variante in lista_variantes) {
    

    cat("Procesando variante:", nombre_variante, "\n")
    
    termino_busqueda <- paste0("SARS-CoV-2 ", nombre_variante, " complete genome OR ", nombre_variante)
    cat("Buscando:", termino_busqueda, "\n")
    
    resultado <- entrez_search(db = "nucleotide", term = termino_busqueda, retmax = retmax)
    
    if (length(resultado$ids) == 0) {
      cat("No se encontraron resultados para", nombre_variante, "\n")
      next
    }
    
    fasta <- entrez_fetch(db = "nucleotide", id = resultado$ids[1], rettype = "fasta")
    cat("✅ Secuencia descargada exitosamente\n")
    

    nombre_archivo <- gsub("[^A-Za-z0-9]", "_", nombre_variante)
    archivo_salida <- paste0(nombre_archivo, ".fasta")
    writeLines(fasta, archivo_salida)
    cat("Secuencia guardada en", archivo_salida, "\n")

    lineas <- unlist(strsplit(fasta, "\n"))
    secuencia <- paste(lineas[-1], collapse = "")
    bases <- strsplit(secuencia, "")[[1]]
    

    total_bases <- length(bases)
    if (total_bases > 1) {
      pares_consecutivos <- paste0(bases[-length(bases)], bases[-1])
      num_gc_pares <- sum(pares_consecutivos == "GC")
      total_pares <- total_bases - 1
      porcentaje_gc_pares <- 100 * num_gc_pares / total_pares
    } else {
      porcentaje_gc_pares <- NA
    }
    

    cat("📏 Longitud de la secuencia:", total_bases, "bases\n")
    cat("🧬 Porcentaje de pares GC:", porcentaje_gc_pares, "%\n")
    

    conteo_bases <- table(factor(bases, levels = c("A", "G", "C", "T")))
    
    barplot(conteo_bases, 
            col = c("skyblue", "salmon", "lightgreen", "orange"), 
            main = paste("Composición de bases -", nombre_variante),
            xlab = "Base", ylab = "Frecuencia")
    

    resultados[[nombre_variante]] <- list(
      fasta = fasta,
      longitud = total_bases,
      gc_porcentaje = porcentaje_gc_pares
    )
  }
  
  return(resultados)
}
resultados <- descargar_y_analizar_varias_secuencias(variantes)
```



##Conclusión 
## Con este código podemos comprender de manera más completa la bilogía molecular con base a 10 variantes de SARS-CoV_2. En todas las variantes anñaizadas se mantienen un rango cercano de 29 mil bases, lo cual confirma que la arquitectura genomica del este virus se mantiene conservada en terminos de tamaño con sus variantes. En cuanto a los pares GC, su porporción fue baja al compararse entre todas la variantes. Las graficas de composición de bases, podemos confirmas que la timina y la adenina son basastentes y esto tambien coincide con los estudios previos sobre el COVID-19.
## La estructura basica se mantiene estable auqnue se presenten mutaciones individuales, los cuales afectan su transmisión. Este tipo analisis son importantes ya que nos ayudan a ver posibles cambio evolutivos en nuestros virus. 
