# Pulse-labeling-reveals-the-tail-end-of-protein-folding-by-proteome-profiling
This is a code repository for "Pulse labeling reveals the tail end of protein folding by proteome profiling" manuscript.  
It contains the following information:  
  
1. python script used to calculate all statistical data involving self-entanglements.  
2. input files in the inpfiles/ directory.  
3. output files in the outpfiles/ directory.  

## Usage  
  

    python stats_v4.5.py [1]   
    [1] = tag for output files    

    Required packages: 
    sys, os, numpy, intertools, scipy, matplotlib, statsmodels  

  
## Input files
  
inpfiles/NEW_clustered_mapped_GE_rep_chain_summary_v2.csv  
contains entanglement information regarding all genes in yeast with a high quality representative PDB  
contains the following columns:   

    0. geneidx - arbitraty counter  
    1. Uniprot gene ID  
    2. number of unique entanglements in the representative structure  
    3. average number of crossings in the representative structure  
    4. length of canonical sequence (amino acids)   

  
inpfiles/genes_with_representative_proteins.txt  
contains information on the genes in the yeast proteome that have representative structures  
the file has the following columns:   
    
    0. Uniprot gene ID
    1. representative PDB
    3. representative chain

inpfiles/CTRLProteinsidentifiedbyMSthatarenotntSPorntCP_UniProt_genes.txt  
contains control genes in yeast identified by MS that are not part of ntSP or ntCP sets  
contains the following columns:   
    
    0. Uniprot gene ID
        

inpfiles/ntSPnewlytranslatedandthermo-sensitiveproteins_UniProt_genes.txt  
contains genes in yeast that are newly translated and temperature sensative ntSP  
contains the following columns:   
    
    0. Uniprot gene ID
        

inpfiles/ntCPmoreinpelletwhennewtranslated_UniProt_genes.txt  
contains genes in yeast with more in pellet when new translated ntCP  
contains the following columns:   
    
    0. Uniprot gene ID
        

inpfiles/ntSP_minus_LiP_genes.txt  
contains genes in yeast ntSP set minus those with significant LiP peptides  
contains the following columns:   
    
    0. Uniprot gene ID
        

inpfiles/ntSPthathaveLiPpeptides_UniProt_genes.txt  
contains genes in yeast ntSP set that have significant LiP peptides  
contains the following columns:   
    
    0. Uniprot gene ID
        

inpfiles/AllLiP_genes.txt  
contains all genes in yeast that had significant LiP peptides  
contains the following columns:   
    
    0. Uniprot gene ID
        

  
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

