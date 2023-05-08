#!/bin/bash
# (C) 2023 Patrick Cahan   

# script to fetch webApps repo, and associated data from S3, and put it in the right place
# N.B.: Execute this from within /srv/shiny-server/


# remove git-related files 
rm -rf .git
rm -f .DS_Store
rm -f README.md


# Fetch Agnostic CellNet Files
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/Hs_expTrain_Jun-20-2017.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/Hs_stTrain_Jun-20-2017.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/heart_GRN_statusQuery.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/heart_broadClassifier100.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/heart_engineeredRef_normalized_expDat_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/heart_engineeredRef_sampTab_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/heart_grnAll.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/heart_studies_iGenes.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/heart_top_TF_display_order.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/heart_trainNormParam.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/hspc_GRN_statusQuery.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/hspc_broadClassifier100.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/hspc_engineeredRef_normalized_expDat_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/hspc_engineeredRef_sampTab_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/hspc_grnAll.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/hspc_studies_iGenes.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/hspc_top_TF_display_order.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/hspc_trainNormParam.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/intestine_colon_GRN_statusQuery.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/intestine_colon_broadClassifier100.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/intestine_colon_engineeredRef_normalized_expDat_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/intestine_colon_engineeredRef_sampTab_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/intestine_colon_grnAll.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/intestine_colon_studies_iGenes.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/intestine_colon_top_TF_display_order.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/intestine_colon_trainNormParam.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/liver_GRN_statusQuery.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/liver_broadClassifier100.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/liver_engineeredRef_normalized_expDat_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/liver_engineeredRef_sampTab_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/liver_grnAll.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/liver_studies_iGenes.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/liver_top_TF_display_order.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/liver_trainNormParam.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/lung_GRN_statusQuery.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/lung_broadClassifier100.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/lung_engineeredRef_normalized_expDat_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/lung_engineeredRef_sampTab_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/lung_grnAll.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/lung_studies_iGenes.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/lung_top_TF_display_order.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/lung_trainNormParam.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/neuron_GRN_statusQuery.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/neuron_broadClassifier100.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/neuron_engineeredRef_normalized_expDat_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/neuron_engineeredRef_sampTab_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/neuron_grnAll.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/neuron_studies_iGenes.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/neuron_top_TF_display_order.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/neuron_trainNormParam.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/skeletal_muscle_GRN_statusQuery.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/skeletal_muscle_broadClassifier100.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/skeletal_muscle_engineeredRef_normalized_expDat_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/skeletal_muscle_engineeredRef_sampTab_all.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/skeletal_muscle_grnAll.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/skeletal_muscle_studies_iGenes.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/skeletal_muscle_top_TF_display_order.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/skeletal_muscle_trainNormParam.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/none_of_the_above_broadClassifier100.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/none_of_the_above_grnAll.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/none_of_the_above_studies_iGenes.rda
wget -P data/ https://cnobjects.s3.amazonaws.com/webApps/agnosticCellNet_web/none_of_the_above_trainNormParam.rda