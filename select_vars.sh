#!/bin/bash

grep -f T2D_Patterns.list T2D_VS_C_AllVars.tab | grep -P '\|MODERATE.*_pred=[PD]|\|HIGH' | grep -P 'BIOBANK.*?=0.0' > T2D_VS_C_SelectedVars.tab

grep -f T2D_Patterns.list T2D_VS_OB+C_AllVars.tab | grep -P '\|MODERATE.*_pred=[PD]|\|HIGH' | grep -P 'BIOBANK.*?=0.0' > T2D_VS_OB+C_SelectedVars.tab

grep -f OBESITY_Patterns.list OBESITY_VS_C_AllVars.tab | grep -P '\|MODERATE.*_pred=[PD]|\|HIGH' | grep -P 'BIOBANK.*?=0.0' > OBESITY_VS_C_SelectedVars.tab

grep -f OBESITY_Patterns.list OBESITY+T2D_VS_C_AllVars.tab | grep -P '\|MODERATE.*_pred=[PD]|\|HIGH' | grep -P 'BIOBANK.*?=0.0' > OBESITY+T2D_VS_C_SelectedVars.tab

for i in *SelectedVars.tab
do
	cat binary_header.tab $i > tmp
	mv tmp $i
done
