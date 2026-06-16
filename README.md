#´---
#´Title:"Análisis de expresión diferencial génica GSE150299"
#´Author: "Valeria Contreras Saldivar"
#´Date:"Cómputo científico, 26P, UAM-C"
#´---

# GSE160299 — Bulk RNA-seq DESeq2 Analysis

> **Dataset:** [GSE160299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160299)  
> **Comparación:** NC (Normal Control) vs PD (Parkinson's Disease)  
> **Tejido:** Plasma / suero humano  
> **Pipeline:** GEO2R-style + DESeq2 con anotación de símbolos génicos

---

## Estructura del repositorio

```
.
├── run_GSE160299.R               ← Script maestro — corre el análisis completo
├── R/
│   ├── 00_config.R               Parámetros globales, directorios, seed
│   ├── 01_packages.R             Instalación y carga de paquetes (CRAN + Bioconductor)
│   ├── 02_plot_helpers.R         Tema ggplot, save wrappers, paletas de color
│   ├── 03_geo_download.R         Descarga de metadata y archivos suplementarios de GEO
│   ├── 04_metadata_parsing.R     Inferencia de condición (NC/PD) y detección de batch
│   ├── 05_count_matrix_parsing.R Parseo de la matriz de conteos y alineación a metadata
│   ├── 06_metadata_library_qc.R  QC de librería y metadata — figuras 01–04
│   ├── 07_deseq2_model.R         Construcción del modelo DESeq2, filtrado y normalización
│   ├── 08_geo2r_qc_plots.R       Plots QC estilo GEO2R — figuras 05–13
│   ├── 09_differential_expression.R  Resultados DE, anotación génica — figuras 14–19
│   ├── 10_heatmaps_gene_profiles.R   Heatmap y perfiles por gen — figuras 20–22
│   └── 11_run_manifest_session.R     Manifiesto de outputs, run_summary, sessionInfo
└── docs/
    └── README.md
```

---

## Inicio rápido

### RStudio

```r
# 1. Clona el repositorio y abre el proyecto en RStudio
# 2. Abre run_GSE160299.R
# 3. Haz clic en Source (o Ctrl+Shift+Enter / Cmd+Shift+Enter)
```

### Terminal

```bash
git clone https://github.com/<tu-usuario>/GSE160299.git
cd GSE160299
Rscript run_GSE160299.R
```

El script maestro detecta su ubicación automáticamente y hace `source()` de cada sección en orden.

---

## Dependencias

### CRAN
`tidyverse` · `ggrepel` · `pheatmap` · `uwot` · `matrixStats` · `scales` · `janitor` · `ggvenn`

### Bioconductor
`GEOquery` · `DESeq2` · `limma` · `AnnotationDbi` · `org.Hs.eg.db` · `apeglm`

> Los paquetes faltantes se instalan automáticamente al correr el análisis (`INSTALL_MISSING <- TRUE` en `00_config.R`).

---

## Outputs generados

Los resultados se guardan en `GSE160299_DESeq2_analysis/` en el directorio de trabajo:

```
GSE160299_DESeq2_analysis/
├── data_raw/                    Archivos suplementarios crudos descargados de GEO
├── results/
│   ├── GSE160299_DESeq2_PD_vs_NC_all_genes.csv
│   ├── GSE160299_DESeq2_PD_vs_NC_significant_FDR0.05_absLFC1.csv
│   ├── GSE160299_DESeq2_normalized_counts.csv
│   ├── GSE160299_DEG_summary.csv
│   ├── GSE160299_metadata_raw_parsed.csv
│   ├── GSE160299_sample_metadata_final.csv
│   ├── GSE160299_gene_annotation_from_matrix.csv
│   ├── GSE160299_raw_counts_parsed.csv
│   ├── output_manifest.csv
│   ├── run_summary.csv
│   └── sessionInfo.txt
└── figures/                     22 figuras en PDF + PNG
```

### Figuras producidas

| # | Nombre | Descripción |
|---|--------|-------------|
| 01 | `metadata_group_counts` | Composición de muestras por grupo |
| 02 | `metadata_overview_tiles` | Mapa de metadata y campo de batch detectado |
| 03 | `qc_library_sizes` | Tamaño de librería por muestra |
| 04 | `qc_detected_genes` | Genes detectados por muestra (conteo ≥ 10) |
| 05 | `geo2r_style_boxplot` | Boxplot de distribución VST |
| 06 | `geo2r_style_density` | Curvas de densidad por muestra |
| 07 | `batch_qc_pca_by_condition` | PCA coloreado por condición |
| 08 | `batch_qc_pca_by_detected_batch` | PCA con forma por batch *(si detectado)* |
| 09 | `batch_qc_pca_after_visual_batch_removal` | PCA post-remoción visual de batch |
| 10 | `geo2r_style_umap` | UMAP de muestras |
| 11 | `sample_distance_heatmap` | Heatmap de distancias muestra-muestra |
| 12 | `deseq2_dispersion_estimates` | Estimados de dispersión DESeq2 |
| 13 | `geo2r_style_mean_variance_trend` | Tendencia media-varianza |
| 14 | `deg_counts_summary` | Genes Up / Down / No significativos |
| 15 | `geo2r_style_volcano` | Volcano plot |
| 16 | `geo2r_style_MD_MA_plot` | MA plot (mean-difference) |
| 17 | `raw_pvalue_histogram` | Histograma de p-valores crudos |
| 18 | `geo2r_style_adjusted_pvalue_histogram` | Histograma de FDR (p ajustados) |
| 19 | `geo2r_style_qq_plot_pvalues` | Q-Q plot de p-valores |
| 20 | `top_de_genes_heatmap` | Heatmap top 50 genes DE (z-score VST) |
| 21 | `geo2r_style_gene_profile_top_genes` | Perfiles de expresión por gen y muestra |
| 22 | `geo2r_style_venn_diagram` | Venn *(solo con ≥ 2 contrastes)* |

---

## Parámetros configurables (`R/00_config.R`)

| Parámetro | Default | Descripción |
|-----------|---------|-------------|
| `GSE_ID` | `"GSE160299"` | ID del dataset en GEO |
| `TEST_LEVEL` | `"PD"` | Grupo de prueba |
| `REF_LEVEL` | `"NC"` | Grupo de referencia |
| `ALPHA` | `0.05` | Umbral de FDR |
| `LFC_CUTOFF` | `1.0` | \|log2FC\| mínimo para significancia |
| `MIN_COUNT` | `10` | Conteos mínimos para filtro de genes |
| `MIN_SAMPLES` | `2` | Muestras mínimas con conteo ≥ `MIN_COUNT` |
| `TOP_HEATMAP_GENES` | `50` | Genes mostrados en el heatmap |
| `TOP_PROFILE_GENES` | `12` | Genes mostrados en perfiles individuales |
| `ADJUST_FOR_BATCH_IF_POSSIBLE` | `TRUE` | Incluir batch en diseño si es factible |

---

## Notas

- El análisis usa la **matriz de conteos crudos** provista en los archivos suplementarios de GEO.
- El batch solo se incluye en el modelo si se detecta una columna de metadata no confundida con la condición y el modelo resultante es de rango completo.
- Para GSE160299 no se espera columna de batch explícita; la evaluación de batch es solo para QC.
- Los símbolos génicos se anotan desde la matriz de conteos (si están disponibles) o mediante `org.Hs.eg.db` (ENSEMBL / Entrez ID como fallback).

---

## Licencia

MIT
