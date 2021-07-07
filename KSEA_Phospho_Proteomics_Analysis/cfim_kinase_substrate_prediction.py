import os
import pdb
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from Bio import SeqIO
import numpy as np
import math
import pickle

#######FUNCTION FOR LONGEST COMMON SUBSTRING BETWEEN 2 STRINGS#######
def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])

    return list(lcs_set)

#######FUNCTION FOR NORMALIZING A DICTIONARY BY THE SUM OF ITS VALUES#######
def normalize(d, target=1.0):
   raw = sum(d.values())
   factor = target/raw
   return {key:value*factor for key,value in d.items()}

all_uniprot_sequences = {}
for entry in SeqIO.parse("uniprot_sprot.fasta", 'fasta'): #####DOWNLOADED UNIPROT DATABASE IN FASTA FILE FORMAT#####
    seq = str(entry.seq)
    na=entry.id.split("|")
    if na[1] in all_uniprot_sequences:
        disp("problem! Line 43")
        pdb.set_trace()
    all_uniprot_sequences[na[1]] = seq
 
 
check_predictions=1 ###check the predictions for known kinases
#####GO THROUGH EACH KINASE IN THE DATABASES######
kinases_db=[]
motifs_db=[]
organisms_db=[]
organisms_target_db=[]
proteins_db=[]
proteins_uniprot_db=[]
uniprot_ids_db=[]
kin_with_site_all_db=[]
ind_phosho_database_constant = 7 # number of amino acids around phosphosite in the database, do not change this!
extra_aa_around_phospho = 0 #this is how many extra amino acids to add around phosphosite

with open('Kinase_Substrate_Dataset') as f: #####DOWNLOADED PHOSPHOSITEPLUS DATABASE FOR KINASE-SUBSTRATE ASSOCIATION#####
    lines=f.readlines()

kin_sites_for_all_kinases={}
for entry in lines:
    headers = entry.strip('\n').split("\t")
    if headers[6] not in all_uniprot_sequences:
        continue
    kinases_db.append(headers[1])
    uniprot_ids_db.append(headers[6])
    index_motif=int(headers[9][1:])
    kin_site = headers[9][0]
    kin_with_site=headers[1] + "_" + kin_site
    try:
        seq_prot = all_uniprot_sequences[headers[6]]
        ind_motif = seq_prot.find(headers[11].upper())
        new_mot = seq_prot[ind_motif-extra_aa_around_phospho:ind_motif] + headers[11][:7].upper() + headers[11][7] + headers[11][8:].upper() + seq_prot[ind_motif+len(headers[11]): ind_motif+len(headers[11]) + extra_aa_around_phospho]
        motifs_db.append(new_mot)
    except:
        motifs_db.append(headers[11][:7].upper()+headers[11][7] + headers[11][8:].upper())
    organisms_db.append(headers[3])
    organisms_target_db.append(headers[8])
    proteins_db.append(headers[4])
    proteins_uniprot_db.append(headers[2])
    kin_with_site_all_db.append(kin_with_site)
    
    if headers[1] not in kin_sites_for_all_kinases:
        kin_sites_for_all_kinases[headers[1]]=kin_site
    else:
        kin_sites_for_all_kinases[headers[1]]=kin_sites_for_all_kinases[headers[1]]+kin_site

##########GENERATE POSITION WEIGHT MATRICES USING A MULTI-DIMENSIONAL ARRAY THAT KEEPS IN THE INFORMATION FOR THE KINASE, PHOSPHORYLATED AMINO ACID, FREQUENCY OF EACH AMINO ACID FOR A KINASE AMONG THE WHOLE DATABASE##########

consensus_length = ind_phosho_database_constant + extra_aa_around_phospho
aminocode = "GALMFWKQESPVICYHRNDTstyhdrk_"
possible_sites = "styhdrk"
#aminocode = "GALMFWKQESPVICYHRNDTsty_"
#possible_sites = "sty"
ind_phosho_database = ind_phosho_database_constant + extra_aa_around_phospho
unique_kinases=[]
count = 0
count2 = 0
all_motifs_aa_count ={}
for kin in kinases_db:
    mot = motifs_db[count2]
    try:
        if uniprot_ids_db[count2]=='P62158':
            uniprot_ids_db[count2]='P0DP23'
        if uniprot_ids_db[count2]=='P01233':
            uniprot_ids_db[count2]='P0DN86'
        if uniprot_ids_db[count2]=='P30443':
            uniprot_ids_db[count2]='P04439'            
        sequence = all_uniprot_sequences[uniprot_ids_db[count2]]
    except:    
        id_hyp=uniprot_ids_db[count2].find("-")
        uni=uniprot_ids_db[count2][0:id_hyp]
        try:
            sequence = all_uniprot_sequences[uni]
        except:
            pdb.set_trace()
    a=sequence.replace(mot.upper(), mot)
    sequence = list(a)
    phospho_amino=mot[ind_phosho_database]
    count2=count2+1
    kin_list = kin
    kin_list.split()
    if kin_list not in unique_kinases:   ########IF IT IS A KINASE THAT DOES NOT HAVE A FREQUENCY MATRIX GENERATED FOR IT########
        if phospho_amino.islower() and len(mot)>=consensus_length*2+1:
            site = possible_sites.find(phospho_amino)
            if site == -1:
                print("Something is wrong here, check this phosphorylation site!")
                pdb.set_trace()
            before_after_AA=[np.zeros((consensus_length*2+1,len(aminocode))), np.zeros((consensus_length*2+1,len(aminocode))), np.zeros((consensus_length*2+1,len(aminocode))), np.zeros((consensus_length*2+1,len(aminocode))), np.zeros((consensus_length*2+1,len(aminocode))), np.zeros((consensus_length*2+1,len(aminocode))), np.zeros((consensus_length*2+1,len(aminocode)))]
            unique_kinases.append(kin)
            #pdb.set_trace()
            before_after_phospho = mot[ind_phosho_database-(consensus_length):ind_phosho_database+consensus_length+1]
            for aa in before_after_phospho:
                ind_aa = aminocode.find(aa)
                if ind_aa==-1:
                    pdb.set_trace()
                try:
                    before_after_AA[site][count][ind_aa] = before_after_AA[site][count][ind_aa] + 1
                except:
                    print(" '_' as amino acid!")
                count=count+1
            all_motifs_aa_count[kin] = before_after_AA
        count = 0
    else:
        ########IF IT IS A KINASE THAT HAS ALREADY A FREQUENCY MATRIX GENERATED FOR IT########
        before_after_AA = all_motifs_aa_count[kin] 
        phospho_amino = mot[ind_phosho_database]
        if phospho_amino.islower() and len(mot)>=consensus_length*2+1:
            site = possible_sites.find(phospho_amino)
            if site == -1:
                print("Something is wrong here, check this phosphorylation site!")
                pdb.set_trace()
            before_after_phospho = mot[ind_phosho_database-(consensus_length):ind_phosho_database+consensus_length+1]
            for aa in before_after_phospho:
                ind_aa = aminocode.find(aa)
                if ind_aa==-1:
                    pdb.set_trace()
                try:
                    before_after_AA[site][count][ind_aa] = before_after_AA[site][count][ind_aa] + 1
                except:
                    print(" '_' as amino acid!")
                count=count+1
            all_motifs_aa_count[kin] = before_after_AA
        count = 0
 
names = ["data/CFIm_KD_Phospho_Peptides.xlsx"]
all_proteins_across_conditions = {}
all_motifs_across_conditions = {}
unique_proteins_all=[]
#pdb.set_trace()
ind_phosho_database_constant = 6 # the number of amino acids around phosposite in our dataset, do not change this!
extra_aa_around_phospho = 1 #this is how many extra amino acids to add around phosphosite
aminocode_upper = "GALMFWKQESPVICYHRNDT_"
####FOR EACH DATA FILE GET ALL PEPTIDES, MOTIFS AND EXTEND THE PEPTIDE IF IT DOES NOT INCLUDE THE WHOLE PHOSPHO SITE#####
#flag_cdk=0
for name in names:

    df = pd.read_excel(name)
    peptides = df["peptide"]
    proteins = df["ac"]
    ptm = df["ptm"]
    updated_orig_peptides=[]
    updated_peptides=[]
    updated_motifs=[]
    updated_proteins=[]
    for i in range(len(peptides)):
        pep=peptides[i]
        pep_orig=peptides[i]
        ptm_pep=ptm[i]
        prot=proteins[i]
        flag_ptm=0
        if ptm_pep.find('|')>-1:
            flag_ptm=1
            all_ptm_all=ptm_pep.split('|')
            for ptm_ptm in all_ptm_all:
                if ptm_ptm.find('Phospho')>-1:
                    ptm_pep=ptm_ptm
                    flag_ptm=2
                    continue
        if flag_ptm==1:
            ptm_pep=all_ptm_all[0]
        all_ptm=ptm_pep.split("]")
        mot_ind=all_ptm[0].strip('[')
        mot_ind=int(mot_ind)
        if True:  #########EXTEND THE MOTIF LENGTH IF REQUIRED BY THE VARIABLE "extra_aa_around_phospho" ##########                                     
            ind_hyp=prot.find("-")
            if ind_hyp>-1:
                prot=prot[:ind_hyp]
            if prot == "Q8R311":
                prot="Q91ZV0"
            try:
                seq_prot = all_uniprot_sequences[prot]
            except:
                print('problematic protein is: ' + prot)
                updated_orig_peptides.append(pep_orig)
                updated_proteins.append(proteins[i])
                updated_peptides.append(pep)
                updated_motifs.append(pep)
                continue

            ind_pep = seq_prot.find(pep.upper())
            ind_mot_from_excel=ind_pep+mot_ind
            if ind_mot_from_excel<7:    
                mo=seq_prot[0:ind_mot_from_excel+6]
            else:
                mo=seq_prot[ind_mot_from_excel-7:ind_mot_from_excel+6]
            if mo=='':
                pdb.set_trace()
            ind_motif = seq_prot.find(mo.upper())
            try:
                new_mot = seq_prot[ind_motif-extra_aa_around_phospho:ind_motif] + mo[:ind_phosho_database_constant].upper() + mo[ind_phosho_database_constant] + mo[ind_phosho_database_constant+1:].upper() + seq_prot[ind_motif+len(mo): ind_motif+len(mo) + extra_aa_around_phospho]
            except:
                updated_orig_peptides.append(pep_orig)
                updated_proteins.append(proteins[i])
                updated_peptides.append(pep)
                updated_motifs.append(mo)
                print(pep)
                continue
            mo=new_mot
            updated_orig_peptides.append(pep_orig)
            updated_proteins.append(proteins[i])
            #pdb.set_trace()
            common=lcs(pep, mo)   #########CHECK IF THE PEPTIDE INCLUDES THE WHOLE PHOPHORYLATION SEQUENCE (IN TERMS OF LENGTH) USED TO CALCULATE POSITION WEIGHT MATRICES, IF NOT EXTEND THE PEPTIDE#####
            if len(''.join(common))<len(mo):
                try:
                    mo.index(''.join(common))
                except:
                    print('problematic protein is: ' + prot)
                    updated_peptides.append(pep)
                    updated_motifs.append(mo)
                    continue
                if mo.index(''.join(common))==0:
                    pep=pep + mo[len(''.join(common)):]
                    updated_peptides.append(pep)
                    updated_motifs.append(mo)
                elif mo.index(''.join(common))>0:
                    pep=mo[0:mo.index(''.join(common))] + pep
                    common=lcs(pep, mo)
                    if len(''.join(common))<len(mo):
                        pep=pep + mo[len(''.join(common)):]
                        updated_peptides.append(pep)
                        updated_motifs.append(mo)
                    else:
                        updated_peptides.append(pep)
                        updated_motifs.append(mo)
                else:
                    print("line 271")
                    pdb.set_trace() 
            else:
                updated_peptides.append(pep)
                updated_motifs.append(mo)                           

    unique_peptides=[]
    motifs_for_unique_peptides=[]
    unique_updated_peptides=[]
    unique_uniprot_ids_db=[]
    #pdb.set_trace()
    for ii in range(len(updated_orig_peptides)):
        pep_test=updated_orig_peptides[ii]
        if pep_test not in unique_peptides:
            unique_peptides.append(pep_test)
            try:
                unique_updated_peptides.append(updated_peptides[ii])
            except:
                pdb.set_trace()
            motifs_for_unique_peptides.append(updated_motifs[ii])
            unique_uniprot_ids_db.append(updated_proteins[ii])
        '''else:
            #pdb.set_trace()
            indd=unique_peptides.index(pep_test)
            if updated_motifs[ii]==motifs_for_unique_peptides[indd]:
                continue
            else:
                #pdb.set_trace()
                unique_peptides.append(pep_test)
                motifs_for_unique_peptides.append(updated_motifs[ii])
                unique_updated_peptides.append(updated_peptides[ii])
                unique_uniprot_ids_db.append(updated_proteins[ii])'''
    #pdb.set_trace() 
    
    updated_orig_peptides = unique_peptides
    updated_peptides = unique_updated_peptides
    updated_motifs = motifs_for_unique_peptides   
    no_of_peptides=len(updated_peptides)
    with open(name+'_updated_orig_peptides_with_site_UNIQUE.p', 'wb') as ffffffkk:
        pickle.dump(updated_orig_peptides, ffffffkk)
    with open(name+'_unique_uniprot_ids_db_with_site_UNIQUE.p', 'wb') as ffffffkkk:
        pickle.dump(unique_uniprot_ids_db, ffffffkkk)

    w_test=open("test_peptides_motifs.txt", "wt")
    for i in range(len(updated_peptides)):
        w_test.write(updated_orig_peptides[i] + "\t" +  updated_peptides[i] + "\t" + updated_motifs[i] + "\n")
    w_test.close()
    wm_background={}              
    for peptide in updated_peptides:
        for aa in peptide:
            if aa not in wm_background:
                wm_background[aa]=1
            else:
                wm_background[aa]=wm_background[aa]+1
    #pdb.set_trace()
    wm_background=normalize(wm_background) #######NORMALIZE THE AMINO ACID FREQUENCIES TO CALCULATE WEIGHT MATRICE FOR BACKGROUND MODEL#######
    #pdb.set_trace()
    ##############GO THROUGH EACH PEPTIDE OF PHOSPHO-PROTEOME AND CALCULATE THE PROBABILITY OF A PEPTIDE TO BE PHOPHORYLATED BY A KINASE THAT HAS PWM, SCAN THROUGH THE WHOLE PEPTIDE WITH A WINDOW SIZE OF THE PWM SIZE. ANY PART OF THE PEPTIDE THAT IS OUT OF THIS WINDOW IS EXPLAINED BY THE BACKGROUND MODEL. ALSO FINALLY ADD THE FULL BACKGROUND MODEL IN CASE THERE ARE NO KINASE THAT CAN PHOSPHORYLATE THE PEPTIDE BASED ON THE PWM##############
    window_sizes_for_peptides=[]
    all_probabilities_for_peptides_kin={}
    all_probabilities_for_peptides_back_cont={}
    all_probabilities_for_peptides_back_model={}
    all_kinase_names=[]
    count_test=0
    for pep_for_em in updated_peptides[:no_of_peptides]:
        print(str(count_test))
        count_test=count_test+1
        pep_len_diff=len(pep_for_em) - 15  #########FIND THE LENGTH DIFFERENCE OF PEPTIDE AND PWM SIZE. THIS DETERMINES HOW MANY WINDOW WILL BE SCANNED FOR THE PEPTIDE#######
        window_sizes_for_peptides.append(pep_len_diff)
        all_scores_for_kinases_all_peptides=0
        background_model=[wm_background[x] for x in pep_for_em] ###BACKGROUND MODEL###
        size_of_pep=[*range(len(pep_for_em))]
        if pep_len_diff==0:   ####IF THE PEPTIDE LENGTH IS EQUAL TO THE PWM SIZE, THEN USE THE FULL BACKGORUND MODEL AND PROBABILITIES FOR EACH KINASES#######
            background_contribution=np.prod(background_model)
            site_count=0
            mot_from_nitish_leo=pep_for_em
            mot_list = list(mot_from_nitish_leo)  
            index_of_motif_nitish_lio = [aminocode_upper.index(aai) for li,aai in enumerate(mot_list)]            
            all_scores_for_kinases=[]
            for key, value in all_motifs_aa_count.items():       
                phos_site_per_kinase = value
                for site_i in range(len(phos_site_per_kinase)):
                    if np.amax(phos_site_per_kinase[site_i])>0 and len(mot_list) == consensus_length*2+1: #####SKIP THE KINASES THAT DOES NOT PHOSPHORYLATE THE AMINO ACID ("styhdrk") BASED ON PWM######
                        pwm_values=np.zeros((consensus_length*2+1,1))
                        kinase_name=key + "_" + possible_sites[site_i]
                        if kinase_name not in all_kinase_names:
                            all_kinase_names.append(kinase_name)
                        for ind_mot in range(len(phos_site_per_kinase[site_i])):  #####GO THROUGH EACH AMINO ACID POSITION#####                          
                            if ind_mot is not consensus_length:
                                #pdb.set_trace()
                                phospho_sites = phos_site_per_kinase[site_i][ind_mot]
                                flag_phospho=0
                                if np.sum(phospho_sites)>0:
                                    #pdb.set_trace()
                                    flag_phospho=1
                                    phospho_sites[0:20]=phospho_sites[0:20]+0.01
                                try:
                                    pwm_values[ind_mot]=phospho_sites[index_of_motif_nitish_lio[ind_mot]]/np.sum(phospho_sites) #####CONVERT FREQUENCIES TO WEIGHTS#####
                                    if flag_phospho==1:
                                        phospho_sites[0:20]=phospho_sites[0:20]-0.01
                                        flag_phospho=0
                                    site_count=site_count+1
                                except:
                                    print("Something is wrong with indexing! line 374")
                                    pdb.set_trace()
                            else:
                                phospho_sites = phos_site_per_kinase[site_i][ind_mot]
                                try:
                                    if mot_list[ind_mot] in kin_sites_for_all_kinases[key]:#mot_list[ind_mot]==possible_sites[site_i].upper():
                                        pwm_values[ind_mot]=1
                                        site_count=site_count+1
                                    else:
                                        pwm_values[ind_mot]=0
                                        site_count=site_count+1
                                except:
                                    print("Something is wrong with indexing! line 386")
                                    pdb.set_trace()          
                        #pdb.set_trace()
                        #pwm_values[consensus_length]=1
                        score_for_kinase=pwm_values.prod()
                        all_scores_for_kinases.append(score_for_kinase)
            ########UPDATE THE DICTIONARIES THAT KEEP TRACK OF THE PROBABILITIES USING THE PWM - PARTIAL BACKGROUND MODEL, FULL BACKGROUND MODEL######
            all_probabilities_for_peptides_kin[pep_for_em] = all_scores_for_kinases
            all_probabilities_for_peptides_back_model[pep_for_em] = background_contribution
            all_probabilities_for_peptides_back_cont[pep_for_em]=1

            continue ######SWITCH TO THE NEXT PEPTIDE######
        
        ######DICTIONARIES THAT KEEP THE PROBABILITIES FOR EACH WINDOW FOR THE PEPTIDE USING PWM - PARTIAL BACKGROUND MODEL, FULL BACKGROUND MODEL##### 
        pep_for_em_likelihood_per_window_kin={}
        pep_for_em_likelihood_per_window_back_cont={}
        pep_for_em_likelihood_per_window_back_model={}      
        for pld in range(pep_len_diff+1):
            pep_to_test=pep_for_em[pld: pld+15]
            ######CALCULATE ALL BACKGROUND PROBABILITIES#####
            background_up_mot=size_of_pep[0:pld]
            background_down_mot=size_of_pep[pld+len(pep_to_test):len(size_of_pep)]
            background_up_mot.extend(background_down_mot)
            background_contribution=[background_model[x] for x in background_up_mot]
            background_contribution=np.prod(background_contribution)
            ######CALCULATE PROBABILITIES USING PWM#####           
            site_count = 0
            mot_from_nitish_leo=pep_to_test
            mot_list = list(mot_from_nitish_leo)  
            index_of_motif_nitish_lio = [aminocode_upper.index(aai) for li,aai in enumerate(mot_list)]
            all_scores_for_kinases=[]
            for key, value in all_motifs_aa_count.items():     
                phos_site_per_kinase = value
                for site_i in range(len(phos_site_per_kinase)):
                    if np.amax(phos_site_per_kinase[site_i])>0 and len(mot_list) == consensus_length*2+1:
                        pwm_values=np.zeros((consensus_length*2+1,1))
                        kinase_name=key + "_" + possible_sites[site_i]
                        if kinase_name not in all_kinase_names:
                            all_kinase_names.append(kinase_name)
                        for ind_mot in range(len(phos_site_per_kinase[site_i])):                            
                            if ind_mot is not consensus_length:
                                phospho_sites = phos_site_per_kinase[site_i][ind_mot]
                                flag_phospho=0
                                if np.sum(phospho_sites)>0:
                                    #pdb.set_trace()
                                    flag_phospho=1
                                    phospho_sites[0:20]=phospho_sites[0:20]+0.01
                                try:
                                    pwm_values[ind_mot]=phospho_sites[index_of_motif_nitish_lio[ind_mot]]/np.sum(phospho_sites)
                                    if flag_phospho==1:
                                        phospho_sites[0:20]=phospho_sites[0:20]-0.01
                                        flag_phospho=0
                                    site_count=site_count+1
                                except:
                                    print("Something is wrong with indexing!")
                                    pdb.set_trace()
                            else:
                                phospho_sites = phos_site_per_kinase[site_i][ind_mot]
                                try:
                                    #pdb.set_trace()
                                    if mot_list[ind_mot] in kin_sites_for_all_kinases[key]:#mot_list[ind_mot]==possible_sites[site_i].upper():
                                        pwm_values[ind_mot]=1
                                        site_count=site_count+1
                                    else:
                                        pwm_values[ind_mot]=0
                                        site_count=site_count+1
                                except:
                                    print("Something is wrong with indexing!")
                                    pdb.set_trace()          
                        #pdb.set_trace()
                        #pwm_values[consensus_length]=1
                        score_for_kinase=pwm_values.prod()
                        all_scores_for_kinases.append(score_for_kinase)
            pep_for_em_likelihood_per_window_kin[pld]=all_scores_for_kinases
            pep_for_em_likelihood_per_window_back_cont[pld]=background_contribution
            #pdb.set_trace()
        ######UPDATE THE DICTIONARIES FOR EACH PEPTIDE######
        all_probabilities_for_peptides_kin[pep_for_em]=pep_for_em_likelihood_per_window_kin
        all_probabilities_for_peptides_back_cont[pep_for_em]=pep_for_em_likelihood_per_window_back_cont
        #pdb.set_trace()
        all_probabilities_for_peptides_back_model[pep_for_em]=np.prod(background_model)
        
    all_kinase_names.append("BACKGROUND")  
    overall_probability_all_kin=0
    overall_probability_per_kin=np.zeros((len(all_kinase_names),1))
    probability_per_kinase_per_peptide=np.zeros((len(updated_peptides[:no_of_peptides]), len(all_kinase_names)))
    #pdb.set_trace()
    for i in range(len(all_kinase_names)): ########GO THROUGH EACH KINASE AND EACH PEPTIDE FOR EACH KINASE TO CALCULATE OVERALL PROBABILITY OF KINASE FOR ALL PEPTIDES#######
        for jj in range(len(updated_peptides[:no_of_peptides])):
            if i == len(all_kinase_names)-1:
                back_model=all_probabilities_for_peptides_back_model[updated_peptides[jj]]
                overall_probability_per_kin[i]=back_model
                probability_per_kinase_per_peptide[jj][i]=back_model
                continue
            
            pep_for_em=updated_peptides[jj]
            kin_probs=all_probabilities_for_peptides_kin[pep_for_em] ####ALL PROBABILITIES FOR PWM MODEL#####
            back_cont=all_probabilities_for_peptides_back_cont[pep_for_em]#####ALL PROBABILITIES FOR PARTIAL BACKGROUND MODEL####
            #back_model=all_probabilities_for_peptides_back_model[pep_for_em]####ALL PROBABILITIES FOR FULL BACKGROUND MODEL####
            #overall_probability_per_kin[i]=back_model#####IN CASE NO KINASE FOR A PEPTIDE####
            if len(kin_probs)==len(all_kinase_names)-1:#####IF THE LENGTH OF THE PEPTIDE IS EQUAL TO THE LENGTH OF PWM####
                overall_probability_per_kin[i]=overall_probability_per_kin[i] + kin_probs[i]*back_cont
                probability_per_kinase_per_peptide[jj][i]=kin_probs[i]*back_cont
            else:#####IF THE LENGTH OF THE PEPTIDE IS LONGER THAN THE LENGTH OF PWM#####
                per_peptide=0
                #pdb.set_trace()
                for pld in range(len(kin_probs)):#####GOT THROUGH EACH WINDOW FOR PROBABILITY CALCULATIONS ALONG THE PEPTIDE#####
                    kin_prob_individual=kin_probs[pld][i]
                    back_cont_individual=back_cont[pld]                    
                    per_peptide=per_peptide+kin_prob_individual*back_cont_individual
                    overall_probability_per_kin[i]=overall_probability_per_kin[i]+kin_prob_individual*back_cont_individual
                probability_per_kinase_per_peptide[jj][i]=per_peptide
    with open(name+'_probability_per_kinase_per_peptide_with_site_unique_peptide_001_UNIQUE.p', 'wb') as f:
        pickle.dump(probability_per_kinase_per_peptide, f)
        
    with open(name+'_all_kinase_names_with_site_unique_peptide_001_UNIQUE.p', 'wb') as ff:
        pickle.dump(all_kinase_names, ff)
        
    with open(name+'_updated_peptides_with_site_unique_peptide_001_UNIQUE.p', 'wb') as fff:
        pickle.dump(updated_peptides, fff)
        
    with open(name+'_updated_motifs_with_site_unique_peptide_001_UNIQUE.p', 'wb') as ffff:
        pickle.dump(updated_motifs, ffff)

    with open(name+'_motifs_db_with_site_unique_peptide_001_UNIQUE.p', 'wb') as fffff:
        pickle.dump(motifs_db, fffff)

    with open(name+'_kinases_db_with_site_unique_peptide_001_UNIQUE.p', 'wb') as ffffff:
        pickle.dump(kinases_db, ffffff)

######GET THE EXPERIMENTAL KINASES FOR EACH PEPTIDE/MOTIF FROM THE DATABASE FOR PREDICTION QUALITY FOR THE EM ALGORITHM########
    updated_posterior_probability=probability_per_kinase_per_peptide
    known_peptides=[]
    indices_known_motifs=[]
    all_known_kinases_in_data={}
    w_em=open(name+"_best_kinases_from_PWM_only_unique_peptide_001_UNIQUE.txt", 'wt')
    w_all=open(name+"_all_predictions_from_PWM_only_unique_peptide_001_UNIQUE.txt", 'wt')
    updated_posterior_probability=probability_per_kinase_per_peptide
    if check_predictions:
        known_peptides=[]
        indices_known_motifs=[]
        all_known_kinases_in_data={}
        motifs_db_upper=[mottt.upper() for mottt in motifs_db]
        for counter in range(len(updated_motifs[:no_of_peptides])): #######CHECK IF THE MOTIFS IN PHOSPHO-PROTEOME IS A KNOWN PHOSPHORYLATION SITE WITH A KNOWN KINASE#####
            try:
                indices = [i for i, x in enumerate(motifs_db_upper) if x == updated_motifs[counter]]
                if len(indices)>0:
                    known_kinases=[kinases_db[x] for x in indices]
                    known_kinases=list(set(known_kinases))
                    posterior_pep=updated_posterior_probability[counter]
                    sorted_indices=sorted(range(len(posterior_pep)), key=posterior_pep.__getitem__)[-50:]
                    sorted_indices.reverse()
                    best_kinases=[all_kinase_names[x] for x, x in enumerate(sorted_indices)]
                    w_em.write(updated_peptides[counter] + "\t" + updated_motifs[counter] + '\t' + str(best_kinases) + "\t" + str(known_kinases) + "\n")
                    #pdb.set_trace()
            except:
                continue       

    w_em.close()       #######FIND THE TOP PREDICTIONS FROM POSTERIOR PROBABILITIES FOR EACH PEPTIDE AND WRITE IT TO A FILE FOR ANALYSIS#######
    for pep_ind in range(len(updated_peptides[:no_of_peptides])):
        pep=updated_peptides[:no_of_peptides][pep_ind]
        posterior_pep=updated_posterior_probability[pep_ind]
        sorted_indices=sorted(range(len(posterior_pep)), key=posterior_pep.__getitem__)[-50:]
        sorted_indices.reverse()
        best_kinases=[all_kinase_names[x] for x, x in enumerate(sorted_indices)]
        w_all.write(pep + "\t" + updated_motifs[pep_ind] + '\t' + updated_orig_peptides[pep_ind] +  "\t" + str(best_kinases) +   "\t" + unique_uniprot_ids_db[pep_ind] + "\n")
    w_all.close()
