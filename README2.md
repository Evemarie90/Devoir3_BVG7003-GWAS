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

Ce script permet de réaliser une **analyse GWAS (Genome-Wide Association Study)** sur des données de génotype et de phénotype. Il est implémenté dans **R** et utilise la bibliothèque **`rMVP`** pour effectuer les calculs statistiques et produire des visualisations comme le **Scree Plot** et les **Manhattan Plots**. 
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

## **2.1 Fichiers d'entrée requis**

Les fichiers suivants sont présents dans la section Data de ce Github:
- **Fichier Génotype** : African_SNPs.hmp.txt
- **Fichier Phénotype** : Phenotype_African.txt

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
   - **`African_SNPs.hmp.txt`** : Fichier contenant les information relatives aux SNPs.  
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

- **Objectif** : Créer un sous-ensemble réduit de 100 échantillons pour des tests rapides afin de valider que le jeu de données est conforme.

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

### **Chargement des Données de Test**

```r
genotype_test <- attach.big.matrix("mvp.hmp.test.geno.desc")
phenotype_test <- read.table("mvp.hmp.test.phe", header = TRUE)
map_test <- read.table("mvp.hmp.test.geno.map", header = TRUE)
```

### **Chargement des Données Complètes**

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
- **IMPORTANT: Pour obtenir les analyses avec les données complètes** : Remplacez les termes "phenotype_test" et "genotype_test" par "phenotype" et "genotype"

---

### **10. Interprétation des Figures et Résultats**

L'interprétation des figures et des résultats dans une analyse GWAS est cruciale pour comprendre les associations entre les **SNPs** (polymorphismes nucléotidiques simples) et le **phénotype** d'intérêt. Voici un examen détaillé des différentes figures utilisées dans ce type d'analyse :

---

### **10.1 Scree Plot (Analyse en Composantes Principales - ACP)**

#### **Objectif :**
Le Scree Plot est un graphique essentiel pour évaluer l'importance des différentes composantes principales (PC) dans une analyse en composantes principales (ACP). L'objectif principal de l'ACP est de réduire la dimensionnalité des données tout en conservant un maximum de variance.

#### **Interprétation :**
- **Axe X** : Représente les composantes principales (PC). Chaque composante principale représente une direction dans l'espace de données qui maximise la variance.
- **Axe Y** : Représente la valeur propre de chaque composante principale, c'est-à-dire la quantité de variance expliquée par cette composante.

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
- La **ligne de seuil** tracée en rouge sur le graphique indique le seuil au-delà duquel les SNPs sont considérés comme **significativement associés** au trait. Les SNPs au-dessus de cette ligne sont ceux que l'on retient pour des analyses plus poussées. Dans notre jeux de données, aucun SNPs est significativement associé à un phénotype quelconque. 
  
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
  
- **Position des SNPs sur l'axe Y** : Les **points** sur l'axe vertical du diagramme indiquent la valeur de **-log10(p-value)** de chaque SNP. Comme dans le Manhattan Plot traditionnel, plus un point est élevé (près du centre), plus l’association du SNP avec le trait est significative.

#### **Interprétation :**
- Les **SNPs significatifs** apparaissent comme des **points se rapprochant du centre du cercle**. Plus un SNP est significatif, plus il sera situé au centre du cercle.
- Les **points éloignés du centre** du diagramme correspondent à des SNPs ayant des **p-values plus élevées**, indiquant qu'ils ne sont pas associés de manière significative au trait d'intérêt.
- La **répartition des SNPs** autour du cercle permet d'identifier rapidement les régions chromosomiques où se concentrent les SNPs significatifs.

**Exemple d'interprétation :**
- Si plusieurs **pics significatifs** sont observés dans une région chromosomique particulière, cela pourrait suggérer une **zone génomique d'intérêt** liée au trait étudié.

---

### **10.4 Graphiques QQplot**
Le QQplot représente la distribution observée des valeurs de \(-\log_{10}(p)\) obtenues à partir des statistiques MLM (Mixed Linear Model) comparée à la distribution théorique attendue sous l'hypothèse nulle.

### Axes :
- **Axe X** : Valeurs attendues de \(-\log_{10}(p)\) sous l'hypothèse nulle.
- **Axe Y** : Valeurs observées de \(-\log_{10}(p)\).

### Éléments du graphique :
- **Ligne rouge** : Distribution théorique sous l'hypothèse nulle.
- **Points bleus** : Valeurs observées de \(-\log_{10}(p)\).
- **Bande bleue** : Intervalle de confiance autour de la distribution théorique.

---

## Interprétation du graphique
1. **Position des points bleus** :
   - Les **points bleus** sont situés **en dessous** de la ligne rouge.

2. **Signification** :
   - Les valeurs observées de \(-\log_{10}(p)\) sont **inférieures** aux valeurs attendues.
   - Cela signifie que les p-valeurs observées sont plus **grandes** que prévu sous l'hypothèse nulle.

3. **Conclusion** :
   - Il n'y a **pas de signal fort** dans les données analysées. Les résultats suggèrent que les associations testées ne sont pas significatives.
   - Le modèle MLM semble **conservateur**, ce qui peut réduire le risque de faux positifs mais augmenter celui de **faux négatifs** (signal réel manqué).

---

### **10.5 Figure de densité des SNPs**

## Objectif  
Cette figure présente la **densité des SNPs** sur les chromosomes, analysée dans des fenêtres de **1 Mb**.

---

## Description des Axes  
- **Axe vertical** : Chromosomes (Chr1 à Chr20).  
- **Axe horizontal** : Positions génomiques le long des chromosomes (en mégabases, Mb).  

---

## Échelle des Couleurs  
La densité des SNPs est représentée par une échelle de couleurs :  
- **Vert foncé à clair** : Faible densité de SNPs (1 à 100 SNPs).  
- **Jaune à orange** : Densité modérée (100 à 232 SNPs).  
- **Rouge** : Haute densité de SNPs (265 à >298 SNPs).  

---

## Interprétation  
- Les **SNPs** sont répartis de manière **hétérogène** sur les chromosomes.  
- Certaines régions présentent une **densité élevée** de SNPs (zones rouges), notamment sur les chromosomes :  
  - **Chr7**, **Chr13**, **Chr16** et **Chr18**.  
- À l'inverse, des régions apparaissent avec une faible densité de SNPs (zones vert foncé).  

---

## Conclusion  
Cette répartition inégale des SNPs pourrait refléter des particularités génomiques telles que :  
- **Régions de sélection** ou soumises à une pression évolutive.  
- **Variabilité génétique accrue** dans certaines zones.  

Les régions de **forte densité** de SNPs constituent des candidats potentiels pour des études approfondies sur la variabilité génétique et la biologie sous-jacente.  

---


## **11. Conclusion**

Ce script, exécuté sous **RStudio**, offre un workflow complet pour l'analyse GWAS. Les résultats obtenus permettent d'identifier des **SNPs significatifs** et des **régions génomiques** d'intérêt pour des études ultérieures.