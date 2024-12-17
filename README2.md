---

## **README - Analyse GWAS avec R et rMVP**

---

## **Table des Matières**

1. [Introduction](#1-introduction)  
2. [Logiciel Requis](#2-logiciel-requis)  
3. [Chargement des Bibliothèques](#3-chargement-des-bibliothèques)  
4. [Définition du Répertoire de Travail](#4-définition-du-répertoire-de-travail)  
5. [Chargement et Préparation du Fichier Hapmap](#5-chargement-et-préparation-du-fichier-hapmap)  
6. [Sélection d'un Sous-Ensemble de 100 Échantillons](#6-sélection-dun-sous-ensemble-de-100-échantillons)  
7. [Conversion des Données pour rMVP](#7-conversion-des-données-pour-rmvp)  
   - 7.1 Données de Test  
   - 7.2 Données Complètes  
8. [Chargement des Données Transformées](#8-chargement-des-données-transformées)  
9. [Analyse en Composantes Principales (ACP)](#9-analyse-en-composantes-principales-acp)  
10. [Filtrage des SNPs avec une Fréquence Allélique Mineure (MAF)](#10-filtrage-des-snps-avec-une-fréquence-allélique-mineure-maf)  
11. [Analyse GWAS avec le Modèle MLM](#11-analyse-gwas-avec-le-modèle-mlm)  
12. [Interprétation des Figures et Résultats](#12-interprétation-des-figures-et-résultats)  
   - 12.1 Scree Plot (ACP)  
   - 12.2 Manhattan Plot  
   - 12.3 Diagrammes Circulaires  
   - 12.4 QQ Plot  
   - 12.5 Densité des SNPs  
13. [Conclusion](#13-conclusion)

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

## **4. Définition du Répertoire de Travail**

```r
setwd("./")  # Définit le répertoire de travail sur le dossier actuel (par exemple, après avoir cloné le dépôt Git)
```

- **Explication** :  
   Cette ligne de code définit le répertoire de travail pour le script R. Le **répertoire de travail** est l'endroit où R va chercher les fichiers d'entrée et où il enregistrera les fichiers de sortie.
   - **`setwd("./")`** : Cette commande définit le répertoire de travail comme étant le **dossier actuel**, c'est-à-dire le dossier où se trouve le script R.  
   - Cela est particulièrement utile lorsque vous avez cloné un dépôt Git contenant le script et les fichiers de données. Il est important de vous assurer que le répertoire de travail est bien défini avant de charger des fichiers ou de générer des résultats.  
   - Cette ligne garantit que le script fonctionne correctement sans avoir à spécifier de chemins absolus pour les fichiers.

---

## **5. Chargement et Préparation du Fichier Hapmap**

```r
hmp_data <- fread("African_SNPs.hmp.txt", header = TRUE)
```

- **Explication** :  
   - **`African_SNPs.hmp.txt`** : Fichier contenant les informations relatives aux SNPs.  
   - **Colonnes 1 à 11** : Informations des SNPs (chromosome, position, etc.).  
   - **Colonnes suivantes** : Génotype des individus.  

---

## **6. Sélection d'un Sous-Ensemble de 100 Échantillons**

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

## **7. Conversion des Données pour rMVP**

### **7.1 Données de Test**

```r
MVP.Data(fileHMP = "test_file.hmp.txt", filePhe = "Phenotype_African.txt", sep.phe = "\t", out = "mvp.hmp.test")
```

### **7.2 Données Complètes**

```r
MVP.Data(fileHMP = "African_SNPs.hmp.txt", filePhe = "Phenotype_African.txt", sep.phe = "\t", out = "mvp.hmp")
```

- **Explication** : Convertit les fichiers Hapmap pour les rendre compatibles avec rMVP.  
- **Sortie** :  
   - Matrice de génotype.  
   - Données de phénotype.  
   - Carte génétique.

Voici la version sans sous-sections pour le README :

---

## 8. Chargement des Données Transformées

Une fois que vous avez préparé vos données et effectué les conversions nécessaires pour rMVP, l'étape suivante consiste à charger les données transformées pour pouvoir les utiliser dans vos analyses. 

Pour charger les données de test, vous devez lire les fichiers suivants :

- **Genotype de Test** : Ce fichier contient les données génétiques sous forme de matrice.
  ```r
  genotype_test <- attach.big.matrix("mvp.hmp.test.geno.desc")
  ```

- **Phénotype de Test** : Ce fichier contient les informations phénotypiques associées aux échantillons de test.
  ```r
  phenotype_test <- read.table("mvp.hmp.test.phe", header = TRUE)
  ```

- **Map de Test** : Ce fichier contient la carte des positions des SNPs pour les échantillons de test.
  ```r
  map_test <- read.table("mvp.hmp.test.geno.map", header = TRUE)
  ```

Pour charger les données complètes, procédez de la manière suivante :

- **Genotype Complet** : Ce fichier contient les données génétiques pour l'ensemble de l'échantillon.
  ```r
  genotype <- attach.big.matrix("mvp.hmp.geno.desc")
  ```

- **Phénotype Complet** : Ce fichier contient les informations phénotypiques pour l'ensemble des échantillons.
  ```r
  phenotype <- read.table("mvp.hmp.phe", header = TRUE)
  ```

- **Map Complet** : Ce fichier contient la carte des positions des SNPs pour l'ensemble des échantillons.
  ```r
  map <- read.table("mvp.hmp.geno.map", header = TRUE)
  ```

Les fonctions `attach.big.matrix()` et `read.table()` sont utilisées pour charger les données depuis les fichiers de votre répertoire de travail.

---

## **9. Analyse en Composantes Principales (ACP)**

```r
genotype_matrix <- as.matrix(genotype)
pca_results <- prcomp(genotype_matrix, scale. = TRUE)
eigenvalues <- pca_results$sdev^2
plot(eigenvalues, type = "b", pch = 19, xlab = "Composante principale", ylab = "Valeur propre", main = "Scree Plot", xlim = c(1, 30))
```

- **Objectif** : Déterminer le nombre de composantes principales pour corriger la structure de population.  
- **Interprétation** : Les **5 premières composantes** expliquent le plus de variance.

---

## **10. Filtrage des SNPs avec une Fréquence Allélique Mineure (MAF)**

Avant d'effectuer l'analyse GWAS, il est important de filtrer les SNPs selon leur fréquence allélique mineure (MAF). Cela permet de supprimer les SNPs rares, qui peuvent être difficiles à interpréter et peuvent introduire du bruit dans les analyses.

Filtrage des SNPs avec MAF < 0.05
Ce filtrage peut être effectué avec le logiciel Tassel 5.0 en appliquant un seuil de 0,05 pour la fréquence allélique mineure. Cela signifie que tous les SNPs dont la MAF est inférieure ou égale à 0,05 seront exclus de l'analyse.

Note importante : Les données partagées dans ce dépôt GitHub sont déjà filtrées pour inclure uniquement les SNPs avec une MAF supérieure à 0,05.


Voici une version détaillée de l'étape "Analyse GWAS avec le Modèle MLM" en format README :

---

## 11. Analyse GWAS avec le Modèle MLM

L'objectif de cette étape est de réaliser une analyse GWAS pour chaque phénotype en utilisant le modèle MLM (Modèle Linéaire Mixte), qui permet de prendre en compte la structure de population et d'effectuer des corrections via des composantes principales (PC). Ce modèle est particulièrement utile pour identifier les associations entre les SNPs et les traits phénotypiques dans un jeu de données.

### Étapes de l'analyse

 **Boucle pour chaque phénotype** :
   La boucle `for(i in 2:ncol(phenotype_test))` permet d'effectuer l'analyse pour chaque colonne de données phénotypiques. On commence à partir de la colonne 2 pour éviter d'inclure la première colonne qui est généralement utilisée pour les identifiants d'échantillons (ou autres informations non phénotypiques).

   ```r
   for(i in 2:ncol(phenotype_test)) {
   ```

 **Appel de la fonction `MVP`** :
   À chaque itération de la boucle, la fonction `MVP` est utilisée pour effectuer l'analyse GWAS en utilisant un modèle MLM. Voici les différents paramètres utilisés dans cette fonction :

   - **`phe = phenotype_test[, c(1, i)]`** : Ce paramètre définit le phénotype à analyser. Il prend la première colonne du fichier (généralement les identifiants des individus) et la colonne i (le phénotype en cours d'analyse).
   - **`geno = genotype_test`** : Ce paramètre définit les données génétiques pour l'analyse. C'est la matrice des génotypes (génome) à utiliser pour identifier les associations avec le phénotype.
   - **`map = map_test`** : Ce paramètre contient la carte des positions des SNPs, nécessaire pour l'analyse afin de connaître la position des SNPs dans le génome.
   - **`nPC.MLM = 5`** : Ce paramètre indique le nombre de composantes principales (PC) à utiliser pour corriger la structure de population dans l'analyse. Ici, 5 PC sont utilisées.
   - **`maxLine = 10000`** : Ce paramètre définit le nombre maximal de lignes à traiter dans l'analyse. Il est utilisé pour optimiser la vitesse du traitement.
   - **`vc.method = "BRENT"`** : Cette option définit la méthode d'estimation de la variance. La méthode "BRENT" est utilisée ici pour maximiser l'efficacité du calcul.
   - **`threshold = 0.05`** : Ce paramètre définit le seuil de significativité pour les SNPs. Un SNP sera considéré comme significatif si sa valeur p est inférieure à 0.05.
   - **`method = "MLM"`** : Ce paramètre spécifie que le modèle utilisé pour l'analyse est un modèle linéaire mixte (MLM).
   - **`file.output = c("pmap", "pmap.signal", "plot", "log")`** : Ces options définissent les types de fichiers de sortie que l'analyse générera, notamment une carte des associations (`pmap`), les signaux significatifs (`pmap.signal`), des graphiques (`plot`) et un fichier log (`log`) pour enregistrer les détails de l'exécution.

   ```r
   imMVP <- MVP(phe = phenotype_test[, c(1, i)],
                geno = genotype_test,
                map = map_test,
                nPC.MLM = 5,
                maxLine = 10000,
                vc.method = "BRENT",
                threshold = 0.05,
                method = "MLM",
                file.output = c("pmap", "pmap.signal", "plot", "log"))
   ```

 **Libération de la mémoire** :
   Après chaque analyse, la fonction `gc()` est appelée pour effectuer un nettoyage de la mémoire et libérer les ressources, ce qui permet d'éviter une surcharge mémoire lorsque vous travaillez avec de grands jeux de données.

   ```r
   gc()
   ```
**Analyse des Données Complètes** :
   Si vous souhaitez utiliser les données complètes pour l'analyse (plutôt que les données de test), il vous suffit de remplacer les objets `phenotype_test` et `genotype_test` par `phenotype` et `genotype` dans le code. Cela permet d'analyser l'ensemble complet des données plutôt que seulement l'échantillon de test.

   - Remplacez :
     ```r
     phenotype_test
     genotype_test
     map_test
     ```
     par :
     ```r
     phenotype
     genotype
     map
     ```

---


### **12. Interprétation des Figures et Résultats**

L'interprétation des figures et des résultats dans une analyse GWAS est cruciale pour comprendre les associations entre les **SNPs** et le **phénotype** d'intérêt. Voici un examen détaillé des différentes figures utilisées dans ce type d'analyse :

---

### **12.1 Scree Plot (Analyse en Composantes Principales - ACP)**

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

## **12.2 Manhattan Plot**

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

### **12.3 Diagrammes Circulaires**

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

### **12.4 Graphiques QQplot**
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

### **12.5 Figure de densité des SNPs**

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


## **13. Conclusion**

Ce script, exécuté sous RStudio, offre un workflow complet pour l'analyse GWAS. Les résultats obtenus permettent d'identifier des **SNP
