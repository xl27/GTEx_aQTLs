# Identifying Genetic Regulatory Variants that Affect Transcription Factor Activity
This repository contatins the core code associated with the manuscript "Identifying Genetic Regulatory Variants that Affect Transcription Factor Activity."

![Outline](https://github.com/xl27/GTEx_aQTLs/blob/main/images/Fig1A.png)

## Inferring individual-specific transcription factor activity
Scripts used for the first step can be found in the folder `TF_activity_inference`
1) Properation of RNA-Seq data from the GTEx project: `Process_GTEx_Count_Matrix.R`	
2) Selection of gene pairs: `Select_Gene_Pairs.Rmd`
3) Generating TF perturbation response signature using ENCODE CRISPRi data: `Perturbation_Signature.Rmd`	
4) Inferring TF activity based on gene-pair model: `Pair_Level_Inference_Model.py`	

### Dependencies:
- R v4.0.1   
- featureCounts v2.0.0    
- Python v3   
- TensorFlow v2.9.1    
- pandas/numpy/sys/sklearn   
 

## Mapping genetic determinants of TF activity (aQTLs)
Scripts used for the second step can be found in the folder `aQTL_mapping`
1) Genome-wide associate studies: `GWAS_mapping.sh`
2) Finemapping: `finemapping.sh`

### Dependencies:
- plink v2.0   
- gcta v1.93.1   
- ldstore v1.1   
- finemap v1.3.1   



