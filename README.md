# **Laboratoire de Contrôle PID avec Feedforward (TCLab)**

## **Objectif du projet**
Ce laboratoire vise à appliquer la théorie du contrôle PID avec feedforward sur un système de contrôle de température appelé **TCLab**. L'objectif est de :
- Comprendre le comportement dynamique du système,
- Identifier les modèles mathématiques de la réponse thermique,
- Implémenter et tester un contrôleur PID en boucle ouverte et fermée,
- Comparer les performances avec et sans feedforward.

---

## **Prérequis avant la première séance**
### **1. Installation logicielle**
- Installer **Anaconda** et les packages **tclab** et **Tkinter** via l'invite de commande Anaconda.
- Tester les fichiers **DBR.py** et **DBR.ipynb** fournis sur la plateforme CLACO.
- Installer **Arduino Software**.

### **2. Familiarisation avec TCLab**
- Lire la documentation.
- Comprendre le branchement et l'utilisation du système.

---

## **Travail expérimental**
### **1. Identification des dynamiques du système**
- Expérimentation sur l’effet du chauffage sur la température.
- Utilisation de modèles mathématiques **FOPDT (Broida)** et **SOPDT (Van der Grinten)**.
- Validation des modèles via des simulations et comparaison des réponses.

### **2. Programmation des contrôleurs**
- Développer un package Python **LAB.py** incluant :
  - Une fonction de correction **Lead-Lag**.
  - Un **PID discret** avec feedforward.
- Tester les performances en simulation (JupyterLab).

### **3. Expériences en boucle ouverte et fermée**
- Effectuer plusieurs essais en enregistrant les réponses du système.
- Comparer les performances du PID seul et avec feedforward.

---

## **Livrables et Évaluation**
### **1. Un rapport PDF contenant :**
- Explications des choix de conception.
- Analyse des résultats expérimentaux avec au moins **7 réponses expérimentales** (boucle ouverte, boucle fermée, avec/sans feedforward).
- Développements mathématiques manuscrits en annexe.

### **2. Un package Python fonctionnel**
- **LAB.py** et fichiers Jupyter associés (**Simulation OLP.ipynb, Simulation CLP PID FF.ipynb**).

### **3. Données expérimentales (.csv)**
- Fichiers correctement nommés et documentés.

### **4. Dépôt des fichiers**
- Sur OneDrive avant **le 2 avril 2025 à 16h**.
- Retards pénalisés.

### **5. Respect du matériel (TCLab) et assiduité obligatoire**

---

## **Méthodes de modélisation utilisées**
### **1. Modèle de Broida (FOPDT)**
- Approximation du système par un modèle du premier ordre avec retard.

### **2. Modèle de Van der Grinten (SOPDT)**
- Modélisation par un système du second ordre avec retard.

### **3. Modèle de Strejc**
- Approximation du système par un modèle d’ordre supérieur.

---

## **Conclusion**
Ce projet demande des compétences en **programmation Python, modélisation mathématique et contrôle des systèmes**. Il est crucial de bien préparer l’environnement de travail, de réaliser les expérimentations correctement et de rédiger un rapport clair et structuré.
