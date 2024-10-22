# CKD microbiome
The two main chunk of codes are in CKD-prevalence1 an CKD-incident, followed by CKD-incidenceDA and CKD-ancombc, and functional and functional-cross. Subsequent analysis were mostly done at different taxonmic levels, and were done in separate markdown files to avoid recalculating certain computationally heavy sections. 

```
cd existing_repo
git remote add origin https://gitlab.com/li-fangyeo/ckd-microbiome.git
git branch -M main
git push -uf origin main
```

## Name
Gut microbiome of chronic kidney disease

## Description
- Cross-sectional anlaysis of the gut microbiome and baseline kidney function using linear regression. 
- Prospective analysis of incident CKD and the gut microbiome using Cox proportional hazards model. 
- Differential abundance anlaysis using ANCOM-BC2. 
- Functional pathway abundance anlaysis from MetaCyc using linear regression (baseline kidney functions) and Cox proportional hazards model (incident CKD).

