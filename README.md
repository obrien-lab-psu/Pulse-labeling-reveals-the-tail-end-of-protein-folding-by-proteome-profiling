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
  
### inpfiles/NEW_clustered_mapped_GE_rep_chain_summary_v2.csv  
contains entanglement information regarding all genes in yeast with a high quality representative PDB  

    contains the following columns:   

    0. geneidx - arbitraty counter  
    1. Uniprot gene ID  
    2. number of unique entanglements in the representative structure  
    3. average number of crossings in the representative structure  
    4. length of canonical sequence (amino acids)   

  
### inpfiles/genes_with_representative_proteins.txt  
contains information on the genes in the yeast proteome that have representative structures  
    
    containst the following columns:   
    
    0. Uniprot gene ID
    1. representative PDB
    3. representative chain

### inpfiles/CTRLProteinsidentifiedbyMSthatarenotntSPorntCP_UniProt_genes.txt  
contains control genes in yeast identified by MS that are not part of ntSP or ntCP sets  
    
    containst the following columns:   
    
    0. Uniprot gene ID
        

### inpfiles/ntSPnewlytranslatedandthermo-sensitiveproteins_UniProt_genes.txt  
contains genes in yeast that are newly translated and temperature sensative ntSP  
    
    containst the following columns:   
    
    0. Uniprot gene ID
        

### inpfiles/ntCPmoreinpelletwhennewtranslated_UniProt_genes.txt  
contains genes in yeast with more in pellet when new translated ntCP  
    
    containst the following columns:   
    
    0. Uniprot gene ID
        

### inpfiles/ntSP_minus_LiP_genes.txt  
contains genes in yeast ntSP set minus those with significant LiP peptides  
    
    containst the following columns:   
    
    0. Uniprot gene ID
        

### inpfiles/ntSPthathaveLiPpeptides_UniProt_genes.txt  
contains genes in yeast ntSP set that have significant LiP peptides  
    
    containst the following columns:   
    
    0. Uniprot gene ID
        

### inpfiles/AllLiP_genes.txt  
contains all genes in yeast that had significant LiP peptides  
    
    containst the following columns:   
    
    0. Uniprot gene ID
        

  
## Output files  
  
A PNG and SVG file are created for each of the five metrics analyzed showing comparison of datasets in input files  
  
### fraction of genes in a dataset that have an entanglement present  
    
    Frac_Entanglement_Presenet_final.png  
    Frac_Entanglement_Presenet_final.svg  
    
  
### average number of unique entanglements identified for each gene in a particular dataset  
    
    AVG_Num_of_Entanglements_final.png  
    AVG_Num_of_Entanglements_final.svg  
    
  
### average number of crossings amongts the unique entanglements idetified for each gene in a dataset  
    
    AVG_Num_of_Crossings_final.png  
    AVG_Num_of_Crossings_final.svg  
    
  
### average number of unique entanglements identified for each gene in a particular dataset normalized by the gene Length  
    
    AVG_Num_of_Entanglements_normL_final.png  
    AVG_Num_of_Entanglements_normL_final.svg  
    
  
### average number of crossings amongts the unique entanglements idetified for each gene in a dataset normalized by the gene Length  
    
    AVG_Num_of_Crossings_normL_final.png  
    AVG_Num_of_Crossings_normL_final.svg  
    
  
### raw statistical data files resulting from bootstrapping and permutation testing  
For each data compared to the control a file containing the bootstraped 95% confidence intervals and p-values determined from permutation testing the differences between sample metric means  
100,000 bootstrapped iterations and permutations are used   
    
    bootstrap_results_AllLiP_X_PME_final.txt  
    bootstrap_results_ntCP_X_PME_final.txt  
    bootstrap_results_ntSP_X_PME_final.txt  
    bootstrap_results_ntSP_minus_LiP_X_PME_final.txt  
    bootstrap_results_ntSPwLiP_X_PME_final.txt  
    
  
