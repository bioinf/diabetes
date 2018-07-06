#!/bin/bash

pseq ./final.pseq v-assoc --phenotype diab > T2D_VS_OB+C_AllVars.tab
pseq ./final.pseq v-assoc --phenotype obesity > OBESITY+T2D_VS_C_AllVars.tab
pseq ./final.pseq v-assoc --phenotype obvc > OBESITY_VS_C_AllVars.tab
pseq ./final.pseq v-assoc --phenotype t2dvc > T2D_VS_C_AllVars.tab
pseq ./final.pseq glm --phenotype bmi --mask maf=0.05:1 > BMI_AllVars.tab
pseq ./final.pseq glm --phenotype glucose --mask maf=0.05:1 > GLUCOSE_AllVars.tab
pseq ./final.pseq glm --phenotype whr --mask maf=0.05:1 > WHR_AllVars.tab
pseq ./final.pseq glm --phenotype triglyceride --mask maf=0.05:1 > TRIGLYCERIDE_AllVars.tab

pseq ./final.pseq assoc --phenotype diab --mask loc.group=gene_name maf=0:0.07 --tests uniq > T2D_VS_OB+CAllLoci_UNIQ.tab
pseq ./final.pseq assoc --phenotype obesity --mask loc.group=gene_name maf=0:0.07 --tests uniq > OBESITY+T2D_VS_C_AllLoci_UNIQ.tab
pseq ./final.pseq assoc --phenotype t2dvc --mask loc.group=gene_name maf=0:0.07 --tests uniq > T2D_VS_C_AllLoci_UNIQ.tab
pseq ./final.pseq assoc --phenotype obvc --mask loc.group=gene_name maf=0:0.07 --tests uniq > OBESITY_VS_C_AllLoci_UNIQ.tab

for i in *UNIQ.tab
do
	grep -P "\tDESC|\tUNIQ" $i > ${i%%.tab}.Final.tsv
done

mkdir raws
mv *UNIQ.tab raws/
