# Guide: Inferring individual-specific transcription factor activity

## Example Datasets:
We use GTEx skeletal muscle RNA-seq data and NR4A1 perturbation data  as examples:
GTEx v8 data: [link](https://gtexportal.org/home/datasets)
ENCODE CRIPSRi data: [link](https://www.encodeproject.org/experiments/ENCSR357LVC/)

## Steps:

1) Preparation of RNA-Seq data from the GTEx project: `Process_GTEx_Count_Matrix.R`	
2) Selection of gene pairs: `Select_Gene_Pairs.Rmd`
3) Generating TF perturbation response signature using ENCODE CRISPRi data: `Perturbation_Signature.Rmd`	
4) Inferring TF activity based on gene-pair model: `Pair_Level_Inference_Model.py`	


## Dependencies:
- R v4.0.1   
- featureCounts v2.0.0    
- Python v3   
- TensorFlow v2.9.1    
- pandas/numpy/sys/sklearn   