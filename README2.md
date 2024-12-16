## **README - Analyse GWAS avec R et rMVP**

---

## **Table des Matières**

1. [Introduction](#1-introduction)  
2. [Logiciel Requis](#2-logiciel-requis)  
3. [Chargement des Bibliothèques](#3-chargement-des-bibliothèques)  
4. [Chargement et Préparation du Fichier Hapmap](#4-chargement-et-préparation-du-fichier-hapmap)  
5. [Sélection d'un Sous-Ensemble de 100 Échantillons](#5-sélection-dun-sous-ensemble-de-100-échantillons)  
6. [Conversion des Données pour rMVP](#6-conversion-des-données-pour-rmvp)  
   - 6.1 Données de Test  
   - 6.2 Données Complètes  
7. [Chargement des Données Transformées](#7-chargement-des-données-transformées)  
8. [Analyse en Composantes Principales (ACP)](#8-analyse-en-composantes-principales-acp)  
9. [Analyse GWAS avec le Modèle MLM](#9-analyse-gwas-avec-le-modèle-mlm)  
10. [Interprétation des Figures et Résultats](#10-interprétation-des-figures-et-résultats)  
    - 10.1 Scree Plot (ACP)  
    - 10.2 Manhattan Plot  
    - 10.3 Diagrammes Circulaires  
    - 10.4 Fichier des SNPs Significatifs  
11. [Conclusion](#11-conclusion)

---

## **1. Introduction**

Ce script permet de réaliser une **analyse GWAS (Genome-Wide Association Study)** sur des données de génotype et de phénotype. Il est implémenté dans **R** et utilise la bibliothèque **`rMVP`** pour effectuer les calculs statistiques et produire des visualisations comme le **Scree Plot** et les **Manhattan Plots**. Téléchargez les fichiers African_SNPs.hmp.txt et Phenotype_African.txt présent dans la section Data pour réaliser l'analyse. 

---

## **2. Logiciel Requis**

- **R** (version ≥ 4.0.0)  
- **RStudio** : Environnement de développement intégré (IDE) recommandé pour R.  
- **Packages nécessaires** :  
   - `data.table`  
   - `rMVP`  
   - `bigmemory`  

### **Installation des Packages**

Pour installer les bibliothèques dans RStudio, exécutez la commande suivante :

```r
install.packages(c("data.table", "bigmemory"))
if (!requireNamespace("rMVP", quietly = TRUE)) {
    install.packages("rMVP", repos = "http://cran.us.r-project.org")
}
```

---

## **3. Chargement des Bibliothèques**

```r
library(data.table)
library(rMVP)
library(bigmemory)
```

- **Objectif** : Charger les bibliothèques pour manipuler les données et exécuter les analyses GWAS.

---

## **4. Chargement et Préparation du Fichier Hapmap**

```r
setwd("./")  # Définit le répertoire de travail actuel
hmp_data <- fread("African_SNPs.hmp.txt", header = TRUE)
```

- **Explication** :  
   - **`African_SNPs.hmp.txt`** : Fichier contenant les génotypes des SNPs.  
   - **Colonnes 1 à 11** : Informations des SNPs (chromosome, position, etc.).  
   - **Colonnes suivantes** : Génotype des individus.  

---

## **5. Sélection d'un Sous-Ensemble de 100 Échantillons**

```r
set.seed(42)
remaining_columns <- colnames(hmp_data)[12:ncol(hmp_data)]
random_columns <- sample(remaining_columns, 100)
selected_columns <- c(colnames(hmp_data)[1:11], random_columns)
test <- hmp_data[, ..selected_columns]
write.table(test, file = "test_file.hmp.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```

- **Objectif** : Créer un sous-ensemble réduit de 100 échantillons pour des tests rapides.

---

## **6. Conversion des Données pour rMVP**

### **6.1 Données de Test**

```r
MVP.Data(fileHMP = "test_file.hmp.txt", filePhe = "Phenotype_African.txt", sep.phe = "\t", out = "mvp.hmp.test")
```

### **6.2 Données Complètes**

```r
MVP.Data(fileHMP = "African_SNPs.hmp.txt", filePhe = "Phenotype_African.txt", sep.phe = "\t", out = "mvp.hmp")
```

- **Explication** : Convertit les fichiers Hapmap pour les rendre compatibles avec rMVP.  
- **Sortie** :  
   - Matrice de génotype.  
   - Données de phénotype.  
   - Carte génétique.

---

## **7. Chargement des Données Transformées**

### **Données de Test**

```r
genotype_test <- attach.big.matrix("mvp.hmp.test.geno.desc")
phenotype_test <- read.table("mvp.hmp.test.phe", header = TRUE)
map_test <- read.table("mvp.hmp.test.geno.map", header = TRUE)
```

### **Données Complètes**

```r
genotype <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype <- read.table("mvp.hmp.phe", header = TRUE)
map <- read.table("mvp.hmp.geno.map", header = TRUE)
```

---

## **8. Analyse en Composantes Principales (ACP)**

```r
genotype_matrix <- as.matrix(genotype)
pca_results <- prcomp(genotype_matrix, scale. = TRUE)
eigenvalues <- pca_results$sdev^2
plot(eigenvalues, type = "b", pch = 19, xlab = "Composante principale", ylab = "Valeur propre", main = "Scree Plot", xlim = c(1, 30))
```

- **Objectif** : Déterminer le nombre de composantes principales pour corriger la structure de population.  
- **Interprétation** : Les **5 premières composantes** expliquent le plus de variance.

---

## **9. Analyse GWAS avec le Modèle MLM**

```r
for(i in 2:ncol(phenotype_test)) {
  imMVP <- MVP(phe = phenotype_test[, c(1, i)],
               geno = genotype_test,
               map = map_test,
               nPC.MLM = 5,
               maxLine = 10000,
               vc.method = "BRENT",
               threshold = 0.05,
               method = "MLM",
               file.output = c("pmap", "pmap.signal", "plot", "log"))
  gc()
}
```

- **Objectif** : Réaliser une analyse GWAS pour chaque phénotype en utilisant le modèle MLM.  
- **Paramètres Clés** :  
   - **`nPC.MLM = 5`** : Correction de structure basée sur 5 PC.  
   - **`threshold = 0.05`** : Seuil de significativité pour les SNPs.
- **Pour obtenir les analyses avec les données complètes** : Remplacez les termes "phenotype_test" et "genotype_test" par "phenotype" et "genotype"

---

## **10. Interprétation des Figures et Résultats**

### **10.1 Scree Plot (ACP)**  
- **Objectif** : Identifier le nombre optimal de PC.  
- **Interprétation** : Le **coude** indique le nombre de PC (5 ici).  

### **10.2 Manhattan Plot**  
- **Axe X** : Position des SNPs.  
- **Axe Y** : -log10(p-value).  
- **Interprétation** : Les pics dépassant la ligne seuil représentent des SNPs associés significativement au phénotype.
La localisation des pics permet d'identifier des régions génomiques d'intérêt.  

### **10.3 Diagrammes Circulaires**  
- Visualisent la distribution des SNPs sur les chromosomes.  
- **Interprétation** : Les SNPs significatifs se trouvent éloignés du centre.

### **10.4 Fichier des SNPs Significatifs**  
- **`pmap.signal`** : Liste des SNPs significatifs avec ID, position et p-value.  
- **Interprétation** : Ces SNPs peuvent être explorés pour identifier des gènes candidats.

---

## **11. Conclusion**

Ce script, exécuté sous **RStudio**, offre un workflow complet pour l'analyse GWAS. Les résultats obtenus permettent d'identifier des **SNPs significatifs** et des **régions génomiques** d'intérêt pour des études ultérieures.
