# afdb-regions-network

Scripts to analyse protein protein similarities based on foldseek output, as described in methods  

Barrio-Hernandez I, Yeo J, Jänes J, Wein T, Varadi M, Velankar S, Beltrao P, Steinegger M. Clustering predicted structures at the scale of the known protein universe. bioRxiv, doi.org:10.1101/2023.03.09.531927 (2023)

SCRIPT 1: partition of Foldseek output for paralelization
SCRIPT 2: filtering edges for evalue<=0.001 plus protein files per partition (loop)
SCRIPT 3: getting ready for hierarchical clustering of regions per protein
SCRIPT 4: hierarchical clustering of regions per protein (loop)
SCRIPT 5: assembling clustering results
SCRIPT 6: recoding the edges (based on clustering results)
SCRIPT 7: assembling the recoded tables, selection of connected components
SCRIPT 8: PFAM annotation of regions, first part, extracting information from database and cuting of files
SCRIPT 9: annotation of regions using pfam (loop)
SCRIPT 10: trimming of the network, clustering of regions from connected components and connecting modules
