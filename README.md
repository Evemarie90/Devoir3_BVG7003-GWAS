***Guide complet pour analyse GWAS ***

**DESCRIPTION**
Ce guide permet d'effectuer l'analyse GWAS en utilisant une jeu de données relatif au génotype et au phénotype soumis au logiciel rMVP. Il inclut des étapes de prétraitement des données, d'exécution d'analyse GWAS et de préparation de graphiques (Manhanttan et QQ). Il permet d'identifier les SNPS significatifs via l'interprétration des graphiques. 

**FONCTIONNALITÉS**
- Explications relatives à la configuration de l'ensemble des données
- Téléchargement des outils R et Rstudio ainsi que du paquet rMVP et des dépendances supplémentaires (ggplot, ggplot2)
- Prétraitement des données
- Exécution de l'analyse GWAS

# Configuration de l'ensemble de données pour l'analyse GWAS

## 1. Accéder aux données

### Projet Figshare
- Télécharger les données relatives au génotype et au phénotype depuis le projet Figshare

## 2. Outils et configuration de l'environnement

### Installation des outils nécessaires
- **R et RStudio** : Suivez les instructions sur CRAN pour installer R et RStudio.
- **rMVP** : Installez le paquet rMVP pour l'analyse GWAS en utilisant la commande suivante dans R :
  ```R
  install.packages("rMVP")
  
**3. Prétraitement des données**
- Formatage des fichiers
- Assurez-vous que les fichiers de génotype et de phénotype sont correctement formatés pour l'analyse.
- Utilisez des scripts R pour convertir les formats si nécessaire.
- Contrôle de qualité
- Filtrez les SNP avec une fréquence allélique mineure (MAF) < 0,05.
- Traitez les données manquantes en utilisant des méthodes appropriées (par exemple, imputation).
- Génération des covariables
- Utilisez l'analyse en composantes principales (ACP) pour la structure de la population ou générez une matrice de parenté.

