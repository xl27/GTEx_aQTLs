# Identifying Genetic Regulatory Variants that Affect Transcription Factor Activity
This repository contatins the core code associated with the manuscript "Identifying Genetic Regulatory Variants that Affect Transcription Factor Activity"

![Outline](https://github.com/xl27/GTEx_aQTLs/blob/main/images/Fig1A.png)

## Inferring individual-specific transcription factor activity

1) Process RNA-Seq data from the GTEx project:
Process_GTEx_Count_Matrix.R	
2) Select gene pairs:
Select_Gene_Pairs.Rmd
3) Generate TF perturbation response signature using ENCODE CRISPRi data:
Perturbation_Signature.Rmd	
4) Infer TF activity based on gene-pair model:
Pair_Level_Inference_Model.py	

## Mapping genetic determinants of TF activity (aQTLs)
