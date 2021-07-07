import pandas as pd
import numpy as np
import pdb

filter_intensity=0
data_type='KD'
#########################INTENSITY_ANALYSIS#########################

data_cfim_seq = pd.read_csv("data/intensities_cfim_phospho_sequence_" + data_type+ ".txt", sep="\t")
data_cfim_seq_orig=data_cfim_seq
data_cfim_seq = data_cfim_seq.groupby(level=0).sum()
data_cfim_seq.columns = [str(col) + '_cfim' for col in data_cfim_seq.columns]
control=data_cfim_seq.iloc[:, 0:3].median(axis=1)
m25=data_cfim_seq.iloc[:, 3:6].median(axis=1)
m68=data_cfim_seq.iloc[:, 6:9].median(axis=1)


all_conditions_medians_cfim=pd.concat([control, m25, m68], axis=1, sort=True)
all_conditions_medians_cfim=all_conditions_medians_cfim.rename(columns={0: "median_control", 1: "median_m25", 2: "median_m68"})

all_conditions_medians_cfim['log2FC_m25_vs_control'] = np.log2(all_conditions_medians_cfim["median_m25"]/all_conditions_medians_cfim["median_control"])
all_conditions_medians_cfim['log2FC_m68_vs_control'] = np.log2(all_conditions_medians_cfim["median_m68"]/all_conditions_medians_cfim["median_control"])
all_conditions_medians_cfim.to_csv('cfim_median_intensities.txt', sep="\t")
#########################ALL_CONDITIONS#########################
all_peptides_muscles_conditions=data_cfim_seq
w=open("CFIm_" + data_type + "_Phospho_Peptides.xlsx_all_predictions_from_PWM_only_unique_peptide_001_UNIQUE.txt", "r")
w1=open("input_for_ksea_cfim_m25_control_" + data_type + ".txt", 'wt')
w2=open("input_for_ksea_cfim_m68_control_" + data_type + ".txt", 'wt')
w1.write("Peptide_ID" + "\t" + "Log2FC" + "\t" + "Predicted_Kinase" +  "\t" + "Proteins"  +  "\t" + 'median_control' +  "\n")
w2.write("Peptide_ID" + "\t" + "Log2FC" + "\t" + "Predicted_Kinase" +  "\t" + "Proteins" +  "\t" + 'median_control' +  "\n")

lines=w.readlines()

count_stats=0
unique_peptides_all=[]
all_kinase_names=[]
for line in lines:
    header=line.strip("\n").replace('[', '').replace(']', '').replace("'", '').split("\t")
    orig_pep=header[2]
    best_predictions=header[3].replace(' ', '').split(',')
    protein=header[4]
    if best_predictions[0]=='BACKGROUND' and len(best_predictions)>1:
        best_pred=best_predictions[1]
        best_pred=best_pred[:len(best_pred)-2]
        all_kinase_names.append(best_pred)
    else:
        best_pred=best_predictions[0]
        best_pred=best_pred[:len(best_pred)-2]
        all_kinase_names.append(best_pred)
    w1.write(orig_pep + "\t" + str(all_conditions_medians_cfim.loc[orig_pep, 'log2FC_m25_vs_control']) + "\t" + best_pred + "\t" + protein + "\t" + str(all_conditions_medians_cfim.loc[orig_pep, 'median_control']) +  "\n")
    w2.write(orig_pep + "\t" + str(all_conditions_medians_cfim.loc[orig_pep, 'log2FC_m68_vs_control']) + "\t" + best_pred + "\t" + protein + "\t" + str(all_conditions_medians_cfim.loc[orig_pep, 'median_control']) +  "\n")
    count_stats=count_stats+1        
w.close()
w1.close()
w2.close()