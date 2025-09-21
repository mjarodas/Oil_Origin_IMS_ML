# Discrimination of olive oil samples accoding to their geographical origin: a machine learning approach using HS-GC-IMS data

---

## 📌 Overview
This repository contains the code and workflows for the discrimination of olive oil samples according to their geographical origin, considering four origins (Spain, Morocco, Portugal and Italy). The analytical technique employed to obtain the data was headspace-gas chromatography-ion mobility spectrometry (HS-GC-IMS), whose conditions had previously been optimised and with the approximation of the use of IMSSS.

The objective of this research was to develop a machine learning-based algorithm for the classification of olive oil samples according to their geographical origin. This would allow automatic learning tools to be implemented in quality control analysis to ensure the origin, and consequently the quality, of olive oils.


## 📂 Project Structure

The repository is structured as follows:

```
├── figures/                               # Generated figures from data analysis
├── scores/                                # Generated scores from data analysis
├── scripts/                               # Scripts for data analysis
│   ├── EDA_script.R                       # Exploratory Data Analysis (EDA)
│   ├── Spectra_visualization_script.R     # Ion Mobility Sum Spectra plotting
│   ├── Unsupervised_algorithms_script.R   # Unsupervised ML (HCA and PCA)
│   ├── RF_CV_script.R                     # Random Forest with 5-fold CV for different pretreatments
│   ├── RF_Train_Test_script.R             # Random Forest with Train/Test split for different pretreatments
├── requirements.txt                       # Required R packages
├── README.md                              # Project documentation
├── LICENSE                                # License file

```

## 🔄 Workflow
The data analysis workflow follows these main steps:

- **Exploratory Data Analysis (EDA)**: Detection of missing values and outliers.

- **Spectra Visualization**: Based on origin of olive oil simples.

- **Unsupervised techniques**: Exploratory assessment of the dataset using Hierarchical Clustering Analysis (HCA) and Principal Component Analysis (PCA)

- **Supervised Machine Learning** Random Forest 5-fold Cross-Validation using different pre-treatments: raw data, maximum normalisation, Savitzky-Golay filter, first derivative with Savitzky-Golay filter, second derivative with Savitzky-Golay filter, Standard Normal Variation, Multiplicative Scatter Correction pretreatments.

- **Supervised Machine Learning** Random Forest with Train (70% of data) and Test (30% of data) using different pre-treatments: raw data, maximum normalisation, Savitzky-Golay filter, first derivative with Savitzky-Golay filter, second derivative with Savitzky-Golay filter, Standard Normal Variation, Multiplicative Scatter Correction pretreatments.


## 🖥️ Software and Dependencies
The analysis is conducted in R (v4.4.0). The required R packages are specified in requirements.txt.

## 📜 License
This project is licensed under the GNU GENERAL PUBLIC License. See LICENSE for deta
