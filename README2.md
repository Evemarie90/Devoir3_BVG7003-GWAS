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
11. [Dépannage](#11-dépannage)
12. [Conclusion](#11-conclusion)

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

1### **10. Interprétation des Figures et Résultats**

L'interprétation des figures et des résultats dans une analyse GWAS est cruciale pour comprendre les associations entre les **SNPs** (polymorphismes nucléotidiques simples) et le **phénotype** d'intérêt. Voici un examen détaillé des différentes figures utilisées dans ce type d'analyse :

---

### **10.1 Scree Plot (Analyse en Composantes Principales - ACP)**

#### **Objectif :**
Le Scree Plot est un graphique essentiel pour évaluer l'importance des différentes composantes principales (PC) dans une analyse en composantes principales (ACP). L'objectif principal de l'ACP est de réduire la dimensionnalité des données tout en conservant un maximum de variance.

#### **Interprétation :**
- **Axe X** : Représente les **composantes principales** (PC). Chaque composante principale représente une direction dans l'espace de données qui maximise la variance.
- **Axe Y** : Représente la **valeur propre** de chaque composante principale, c'est-à-dire la quantité de variance expliquée par cette composante.

Le **Scree Plot** permet de déterminer combien de composantes principales doivent être retenues dans l'analyse. La **règle du coude** (ou *elbow rule*) est utilisée pour choisir le nombre de composantes principales à conserver.

- Le **coude** du graphique est généralement visible lorsque la pente des valeurs propres diminue de manière significative. Avant le coude, chaque composante principale apporte une réduction substantielle de la variance, mais après ce point, l'ajout de nouvelles composantes principales contribue beaucoup moins à expliquer la variance.
  
Dans notre exemple :
- **5 composantes principales** sont retenues après l'observation du coude dans le graphique. Cela signifie que les 5 premières composantes expliquent suffisamment la variance des données, et il n'est pas nécessaire d'inclure davantage de composantes.

---

### **10.2 Manhattan Plot**

#### **Objectif :**
Le **Manhattan Plot** est une figure utilisée pour visualiser les résultats d'une analyse GWAS. Elle permet d'identifier les **SNPs associés** au trait d'intérêt en affichant les p-values des tests statistiques.

#### **Structure du Graphique :**

- **Axe X :** Représente la **position génomique** des SNPs sur les chromosomes. Chaque point sur l'axe X correspond à un SNP particulier. Les SNPs sont regroupés par chromosome, et chaque chromosome est représenté par un secteur sur l'axe X. Cela permet de voir la répartition des SNPs sur le génome.
  
  - **Position des SNPs** : Les SNPs sont ordonnés selon leur position sur le chromosome. Cela permet d'identifier les régions génomiques où se concentrent les SNPs analysés.
  
- **Axe Y :** Représente les **-log10(p-values)** des tests d'association réalisés pour chaque SNP. Les **p-values** sont transformées en **valeurs -log10** pour faciliter l'interprétation. Une **valeur élevée sur l'axe Y** indique une **faible p-value**, et donc une **association significative** entre le SNP et le phénotype.

#### **Interprétation :**
- Les **points élevés** sur l'axe Y correspondent aux SNPs ayant des **p-values faibles**, c'est-à-dire des SNPs qui sont fortement associés au trait d'intérêt.
- La **ligne de seuil** est souvent tracée sur le graphique (par exemple à **-log10(0.05)** ou une valeur similaire), indiquant un seuil au-delà duquel les SNPs sont considérés comme **significativement associés** au trait. Les SNPs au-dessus de cette ligne sont ceux que l'on retient pour des analyses plus poussées.
  
**Exemple d'interprétation :**
- Un pic de SNPs au-dessus du seuil peut suggérer une **région génomique d'intérêt**, potentiellement impliquée dans le trait étudié.
- L'intensité et la localisation des pics peuvent également fournir des indices sur des **gènes candidats** ou des **régions génétiques** importantes.

---

### **10.3 Diagrammes Circulaires**

#### **Objectif :**
Les **diagrammes circulaires** sont utilisés pour visualiser de manière compacte et esthétique la **répartition des SNPs** significatifs sur le génome, tout en maintenant une structure chromosomique ordonnée. C’est une variante du **Manhattan Plot** classique, mais avec un agencement circulaire.

#### **Structure du Graphique :**

- **Disposition Circulaire des Chromosomes :** Le diagramme est constitué de plusieurs **secteurs** représentant chacun un **chromosome**. Chaque secteur est une section du cercle où la position des SNPs est représentée en fonction de leur **position physique** sur le chromosome.
  
  - **Positions des SNPs** : Sur chaque chromosome, les SNPs sont représentés le long de l'arc du cercle, suivant leur position sur le chromosome. Les SNPs sont espacés de manière uniforme le long de chaque secteur.
  
- **Position des SNPs sur l'axe Y** : Les **points** sur l'axe vertical du diagramme indiquent la valeur de **-log10(p-value)** de chaque SNP. Comme dans le Manhattan Plot traditionnel, plus un point est élevé, plus l’association du SNP avec le trait est significative.

#### **Interprétation :**
- Les **SNPs significatifs** apparaissent comme des **points éloignés du centre du cercle**. Plus un SNP est significatif, plus il sera situé près de la périphérie du cercle.
- Les **points situés près du centre** du diagramme correspondent à des SNPs ayant des **p-values plus élevées**, indiquant qu'ils ne sont pas associés de manière significative au trait d'intérêt.
- La **répartition des SNPs** autour du cercle permet d'identifier rapidement les régions chromosomiques où se concentrent les SNPs significatifs.

**Exemple d'interprétation :**
- Si plusieurs **pics significatifs** sont observés dans une région chromosomique particulière, cela pourrait suggérer une **zone génomique d'intérêt** liée au trait étudié.

---

### **10.4 Fichier des SNPs Significatifs**

#### **Objectif :**
Le fichier **`pmap.signal`** contient les **détails des SNPs significatifs** après l'analyse GWAS. Ce fichier est crucial pour une exploration plus approfondie des SNPs associés et pour la recherche de gènes candidats.

#### **Structure du fichier :**
Le fichier contient généralement les informations suivantes :
- **ID du SNP** : Un identifiant unique pour chaque SNP analysé.
- **Position** : La position physique du SNP sur le chromosome.
- **P-value** : La p-value associée au SNP, représentant l'intensité de l'association avec le phénotype.
  
#### **Interprétation :**
- **SNPs significatifs** : Les SNPs dont la **p-value** est inférieure au seuil de significativité défini sont considérés comme **significativement associés** au trait.
- **Identification de gènes candidats** : Ces SNPs peuvent être situés près de gènes connus ou dans des régions du génome qui méritent une exploration approfondie pour comprendre leurs rôles potentiels dans le trait d'intérêt.

**Exemple d'interprétation :**
- Un SNP situé près d'un gène spécifique dans le fichier **`pmap.signal`** pourrait être une **cible potentielle** pour une étude fonctionnelle, permettant de vérifier si ce gène joue un rôle dans l'expression du phénotype étudié.

---

## **11. Dépannage ##
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

## **12. Conclusion**

Ce script, exécuté sous **RStudio**, offre un workflow complet pour l'analyse GWAS. Les résultats obtenus permettent d'identifier des **SNPs significatifs** et des **régions génomiques** d'intérêt pour des études ultérieures.
