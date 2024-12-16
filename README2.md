# **README - Explication Détaillée du Script GWAS**

---

## **1. Introduction**

Ce script permet de réaliser une **analyse GWAS** (Genome-Wide Association Study) sur des données de génotype et de phénotype. Il utilise des bibliothèques R pour charger, traiter et analyser les données, générer des visualisations (ACP, Manhattan Plot, etc.), et identifier les SNPs significatifs. Ce document décrit **chaque partie du script** en détail et explique comment interpréter les résultats.

---

## **2. Description Détaillée du Script**

### **2.1 Chargement des bibliothèques**

```r
library(data.table)
library(rMVP)
library(bigmemory)
```

- **`data.table`** : Utilisée pour charger rapidement de grands fichiers texte comme les fichiers Hapmap.  
- **`rMVP`** : Bibliothèque dédiée pour les analyses GWAS avec des méthodes avancées comme **MLM (Mixed Linear Model)**.  
- **`bigmemory`** : Facilite la manipulation de matrices volumineuses, essentielles pour les analyses de génotype.

---

### **2.2 Chargement et Préparation du Fichier Hapmap**

```r
setwd("./")  # Définit le répertoire de travail actuel
hmp_data <- fread("African_SNPs.hmp.txt", header = TRUE)
```

- **Explication** :  
   Le fichier **`African_SNPs.hmp.txt`** est chargé. Ce fichier Hapmap contient :  
   - **Colonnes 1 à 11** : Métadonnées sur les SNPs (ID, position chromosomique, etc.).  
   - **Colonnes suivantes** : Génotype des individus.  
- **Objectif** : Préparer les données pour l'analyse GWAS.

---

### **2.3 Sélection d'un Sous-Ensemble de 100 Échantillons**

```r
set.seed(42)
remaining_columns <- colnames(hmp_data)[12:ncol(hmp_data)]
random_columns <- sample(remaining_columns, 100)
selected_columns <- c(colnames(hmp_data)[1:11], random_columns)
test <- hmp_data[, ..selected_columns]
write.table(test, file = "test_file.hmp.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```

- **Explication** :  
   - **`set.seed(42)`** garantit la reproductibilité des résultats.  
   - 100 individus sont sélectionnés aléatoirement pour créer un sous-ensemble de test, ce qui accélère les analyses initiales.  
- **Objectif** : Permet de tester rapidement le workflow GWAS avant d'utiliser l'ensemble des données.

---

### **2.4 Conversion des Données pour MVP**

#### **Données de test**

```r
MVP.Data(fileHMP = "test_file.hmp.txt", filePhe = "Phenotype_African.txt", sep.phe = "\t", out = "mvp.hmp.test")
```

#### **Données complètes**

```r
MVP.Data(fileHMP = "African_SNPs.hmp.txt", filePhe = "Phenotype_African.txt", sep.phe = "\t", out = "mvp.hmp")
```

- **Explication** :  
   - **`fileHMP`** : Fichier de génotype au format Hapmap.  
   - **`filePhe`** : Fichier contenant les phénotypes (données associées aux individus).  
   - **`out`** : Préfixe des fichiers de sortie (génotype, phénotype, carte génétique).  
- **Objectif** : Transformer les données brutes pour les rendre compatibles avec **`rMVP`**.

---

### **2.5 Chargement des Données Transformées**

#### **Données de test**

```r
genotype_test <- attach.big.matrix("mvp.hmp.test.geno.desc")
phenotype_test <- read.table("mvp.hmp.test.phe", header = TRUE)
map_test <- read.table("mvp.hmp.test.geno.map", header = TRUE)
```

#### **Données complètes**

```r
genotype <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype <- read.table("mvp.hmp.phe", header = TRUE)
map <- read.table("mvp.hmp.geno.map", header = TRUE)
```

- **Explication** :  
   - **`genotype`** : Matrice de génotype stockée en mémoire pour une manipulation rapide.  
   - **`phenotype`** : Données de phénotype pour chaque individu (valeurs mesurées).  
   - **`map`** : Carte génétique indiquant la position des SNPs sur les chromosomes.  
- **Objectif** : Préparer les données pour l'analyse GWAS.

---

### **2.6 Analyse en Composantes Principales (ACP)**

```r
genotype_matrix <- as.matrix(genotype)
pca_results <- prcomp(genotype_matrix, scale. = TRUE)
eigenvalues <- pca_results$sdev^2
plot(eigenvalues, type = "b", pch = 19, xlab = "Composante principale", ylab = "Valeur propre", main = "Scree Plot", xlim = c(1, 30))
```

- **Objectif** : Identifier le **nombre optimal de composantes principales** pour corriger la structure de la population.  
- **Interprétation du Scree Plot** :  
   - **X** : Composantes principales (PC).  
   - **Y** : Valeurs propres (**eigenvalues**).  
   - Le **"coude"** du graphique indique le nombre minimal de PC à inclure (ici, 5 composantes sont choisies).  

---

### **2.7 Analyse GWAS avec le Modèle MLM**

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

- **Paramètres** :  
   - **`phe`** : Données de phénotype (une colonne à la fois).  
   - **`geno`** : Matrice de génotype.  
   - **`map`** : Carte génétique.  
   - **`nPC.MLM = 5`** : Nombre de composantes principales (ACP).  
   - **`method = "MLM"`** : Utilise le modèle mixte linéaire pour corriger la structure de la population.  
   - **`threshold = 0.05`** : Seuil de significativité pour les SNPs.

---

## **3. Interprétation des Figures et Résultats**

### **3.1 Scree Plot (ACP)**

- **Objectif** : Déterminer le nombre de PC pour corriger la population.  
- **Interprétation** : Le coude de la courbe indique **5 PC** comme optimal.

---

### **3.2 Manhattan Plot**

- **Description** : Graphique des p-values des SNPs.  
   - **Axe X** : Position des SNPs sur les chromosomes.  
   - **Axe Y** : -log10(p-value).  
- **Interprétation** :  
   - Les **pics dépassant la ligne seuil** représentent des SNPs associés significativement au phénotype.  
   - La localisation des pics permet d'identifier des régions génomiques d'intérêt.

---

### **3.3 Diagrammes Circulaires**

- **Objectif** : Montrer la distribution des SNPs significatifs sur les chromosomes.  
- **Interprétation** :  
   - Les SNPs significatifs apparaissent éloignés du centre.  
   - Une **concentration de SNPs** dans une région spécifique indique un locus candidat.

---

### **3.4 Fichier des SNPs Significatifs**

- **Fichier `pmap.signal`** : Contient les SNPs significatifs avec :  
   - **ID du SNP**.  
   - **Position génomique**.  
   - **p-value**.  
- **Interprétation** : Ces SNPs peuvent être analysés dans des bases de données pour identifier des gènes associés.

---

## **4. Conclusion**

Ce script réalise une **analyse complète GWAS**. Les figures et fichiers générés permettent d'identifier des **SNPs candidats** associés à un phénotype.  
- **Scree Plot** : Nombre de PC pour corriger la population.  
- **Manhattan Plot et diagrammes circulaires** : Visualisation des SNPs significatifs.  
- **Fichier signal** : Liste des SNPs à analyser plus en détail.  

**Utilisez ces résultats pour explorer les loci et gènes potentiellement impliqués.**
