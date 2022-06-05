# Pulse-labeling-reveals-the-tail-end-of-protein-folding-by-proteome-profiling
This is a code repository for "Pulse labeling reveals the tail end of protein folding by proteome profiling" manuscript.  
It contains the following information:  
  
1. python script used to calculate all statistical data involving self-entanglements.  
2. input files in the inpfiles/ directory.  
3. output files in the outpfiles/ directory.  

## Usage  
'''python stats_v4.5.py [1]     
[1] = tag for output files    
'''  
  
## Input files
  
inpfiles/NEW_clustered_mapped_GE_rep_chain_summary_v2.csv  
| contains the following columns  
0. geneidx - arbitraty counter  
1. Uniprot gene ID  
2. number of unique entanglements in the representative structure  
3. average number of crossings in the representative structure  
4. length of canonical sequence (amino acids) |  

  
inpfiles/genes_with_representative_proteins.txt  
inpfiles/CTRLProteinsidentifiedbyMSthatarenotntSPorntCP_UniProt_genes.txt  
inpfiles/ntSPnewlytranslatedandthermo-sensitiveproteins_UniProt_genes.txt  
inpfiles/ntCPmoreinpelletwhennewtranslated_UniProt_genes.txt  
inpfiles/ntSP_minus_LiP_genes.txt  
inpfiles/ntSPthathaveLiPpeptides_UniProt_genes.txt  
inpfiles/AllLiP_genes.txt  
  
## Output files  
  
AVG_Num_of_Crossings_final.png  
AVG_Num_of_Crossings_final.svg  
AVG_Num_of_Crossings_normL_final.png  
AVG_Num_of_Crossings_normL_final.svg  
AVG_Num_of_Entanglements_final.png  
AVG_Num_of_Entanglements_final.svg  
AVG_Num_of_Entanglements_normL_final.png  
AVG_Num_of_Entanglements_normL_final.svg  
Frac_Entanglement_Presenet_final.png  
Frac_Entanglement_Presenet_final.svg  
bootstrap_results_AllLiP_X_MSCntrl_final.txt  
bootstrap_results_ntCP_X_MSCntrl_final.txt  
bootstrap_results_ntSP_X_MSCntrl_final.txt  
bootstrap_results_ntSP_minus_LiP_X_MSCntrl_final.txt  
bootstrap_results_ntSPwLiP_X_MSCntrl_final.txt  
