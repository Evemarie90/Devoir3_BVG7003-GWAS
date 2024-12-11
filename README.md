# GUIDE COMPLET POUR ANALYSE GWAS

## Description
Ce guide détaille les étapes nécessaires pour effectuer une analyse GWAS à l'aide d'un jeu de données génotypiques et phénotypiques. Il utilise le package **rMVP** dans R pour l'exécution de l'analyse, et inclut des étapes pour le prétraitement des données, la génération de graphiques (Manhattan et QQ), et l'identification des SNPs significatifs.

## Table des matières
1. [Prérequis](#prérequis)
2. [Configuration des outils et données](#configuration-des-outils-et-données)
3. [Prétraitement des données](#prétraitement-des-données)
4. [Exécution de l'analyse GWAS](#exécution-de-lanalyse-gwas)
5. [Résultats attendus](#résultats-attendus)
6. [Exemple d'interprétation](#exemple-dinterprétation)

## Prérequis
- R (version 4.0+)
- RStudio (version 1.4+)
- Connexion Internet pour télécharger les données et packages nécessaires.

## Configuration des outils et données
1. Installez R et RStudio en suivant les instructions disponibles sur [CRAN](https://cran.r-project.org/).
2. Installez les packages nécessaires dans R :
   ```r
   install.packages(c("rMVP","ggplot2", "data.table" "readxl", "dplyr"))

## Prétraitement des données

1. Chargez vos fichiers de génotype et de phénotype.
2. Filtrez les SNPs avec une fréquence allélique mineure (MAF) inférieure à 0,05 :
```r
genotype <- genotype[genotype$MAF >= 0.05, ]
```
3. Générez des covariables via une analyse en composantes principales (ACP) pour corriger la structure de la population.

## Exécution de l'analyse GWAS

1. Chargez les données dans rMVP :
```r
MVP.Data(fileGenotype = "genotype.txt", filePhenotype = "phenotype.txt")
```
2. Exécutez l'analyse GWAS
```r
MVP.GWAS(file = "data", covariates = cov_matrix)
```

## Résultats attendus

Un graphique de Manhattan montrant les SNPs significatifs.
Un graphique QQ montrant l'absence ou la présence d'inflation génomique.

## Exemple d'interprétation

Les SNPs au-dessus du seuil dans le graphique de Manhattan indiquent des associations significatives.
Une distribution linéaire dans le graphique QQ suggère une bonne correction des valeurs P.

## Dépannage ##
Advenant que le téléchargement du fichier hmp entraine la modification de l'alignement des colonnes de ce dernier vous pouvez convertir votre fichier hmp en fichier excel (copier-coller les données). Ensuite, vous pouvez utiliser le script ci-dessous pour re-convertir vos données en hmp :

#### Lire le fichier excel dans R ####
```r
library(readxl)

excel_data <- read_excel("excel_dataset.xlsx")
```

#### Écrire les données dans un fichier hmp ####
```r
Rstudio write.table(excel_data, file = "output_file.hmp.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
```


