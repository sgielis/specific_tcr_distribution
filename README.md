# Epitope-specific TCR repertoires and prediction model performances

This github repository provides the code for the manuscript 'Revealing the hidden sequence distribution of epitope-specific TCR repertoires and its influence on machine learning model performance'. The code is structured according to the subtitles in the results section of the paper. 


### Overview of the collected data  

#### (1) Summarize the number of epitopes and training data
- notebook: data_description/data_description.ipynb  
- environment: description_env 
- input: data/parsed/tcrex_models.csv  
- output: results/data_description  

#### (2) Collecting all epitope-specific TCRs in one dataframe  
- notebook: data_collection/assemble_tcrex_data.ipynb  
- environment: raptcr_env  
- input: data/parsed/TCRex_data  
- output:  
  - data/final/all_tcrs.tsv: contains all TCRs for every epitope after TCRex parsing  
  - data/final/unique_CDR3s.tsv: contains every CDR3 beta once with a list of the epitopes and V/Jgenes  


### Epitope-specific TCR repertoires are an ensemble of different TCR clusters  

#### (1) Cluster the TCRex training data
- notebook: TCRex_space/tcrex_clustering/cluster_tcrex_data.ipynb  
- input: data/final/unique_CDR3s.tsv
- environment: raptcr_env
- output:
  - results/tcrex_clustering/tcrex_clusters.tsv
  - results/tcrex_clustering/clusters_sizes.tsv
  - results/tcrex_clustering/pure_clusters.tsv

#### (2) Create a UMAP of the TCRex training data
- notebook: tcrex_space/umap/umap_tcrex.ipynb  
- input: data/final/unique_CDR3s.tsv  
- environment: raptcr_env  
- output: UMAP of the TCRex space in notebook  


### Epitope-specific TCR clusters are spread out over TCR space 

#### (1) Cluster epitope-specific TCR repertoires separately 
- notebook: epitope_specific_clustering/cluster_epitope_data.ipynb  
- environment: raptcr_env  
- input: /data/final/all_tcrs.tsv  
- output: results/epitope_specific_clustering/epitope_specific_clusters.tsv  

#### (2) Determine the number of unique motifs and the number of overlapping motifs
- notebook: epitope_specific_clustering/unique_motifs.ipynb  
- environment: raptcr_env  
- input: epitope_specific_clustering/epitope_specific_clusters.tsv  
- Output:  
  - results/epitope_specific_clustering/unique_motifs.tsv: one row per motif, multiple epitopes separated by comma  
  - results/epitope_specific_clustering/overlapping_motifs.tsv: table with only the motifs shared between epitopes  

#### (3) Cluster the epitope-specific TCR motifs
- notebook: epitope_specific_clustering/cluster_motifs.ipynb  
- environment: raptcr_env    
- input:epitope_specific_clustering/unique_motifs.tsv   
- output: /results/epitope_specific_clustering/cluster_motifs/motif_clusters.tsv    


### Negative TCRs can differ in only one amino acid with epitope-specific TCRs 

#### (1) Cluster positive and negative data together
- task: For every epitope, cluster positive(1) and negative(0) data together and count how many clusters contain both positive and negative TCRs. Shared CDR3 beta sequences were removed before clustering  
- notebook: background/cluster_together.ipynb   
- environment: raptcr_env  
- input: data/parsed/tcrex_data 
- output: results/background/shared_clusters.tsv: Table with the number of shared clusters and total nr of clusters

#### (2) Cluster positive and negative data separatly
- task: For every epitope, the positive(1) and negative(0) data was clustered separately. Also here, shared CDR3 beta sequences were removed before clustering. For every cluster the motifs were defined. Finally, motifs shared with positive and negative clusters were identified	
- notebook: background/cluster_separately.ipynb
- environment: raptcr_env
- input: data/parsed/TCRex_data 
- output: results/background/shared_motifs.tsv : Table with a list of the shared motifs between positive and negative clusters 


### TCRs originating from different LIgO seeds can have similar sequences

#### (1) Perform LIgO simulations
- environment: ligo_env
- output: results/ligo_simulations: this folder also contains the specs.yaml files that were used to simulate the repertoires

#### (2) Cluster LIgO repertoires and their motifs
- notebook: ligo/clustering/cluster_ligo_results.ipynb 
- environment: raptcr_env
- input:
  - seeds.csv
  - results/ligo_simulations
- output
  - results/ligo_results/clustering_statistics.tsv: table of repertoire clustering
  - results/ligo_results/list_of_unique_motifs.tsv: list of all motifs:
  - results/ligo_results/impure_motif_clusters.tsv: list of impure clusters

#### (3) UMAP of LIgO simulated repertoires
- notebook: ligo/umap/umap_ligo.ipynb 
- environment: raptcr_env
- input: results/ligo_simulations
- output: figures are shown directly in the notebook


### TCRs in between clusters have low generation probability 

#### (1) Simulate large repertoires using only one seed per repertoire
- Task: From the 8 simulated repertoires in results/ligo_simulations, keep one seed per yaml file and generate 3000 TCRs instead of 300 TCRs 
   - Simulation 1: WTGEKHE  
   - Simulation 2: KGTGLYNE  
   - Simulation 3: SLAVGGYE  
   - Simulation 4: PIPPYNE  
   - Simulation 5: PFAGADT  
   - Simulation 6: TPWGGSYE  
   - Simulation 7: VVGRGAYNE  
   - Simulation 8: SPLLAGGPYE  
- environment: ligo_env  
- output: results/ligo_simulations_one_seed: this folder also contains the specs.yaml files that were used to simulate the repertoires  

#### (2) Collect all simulated data into one df and calculate pgen based on junction_aa
- notebook: generation_prob/simulated_repertoires/Parse_simulated_data.ipynb  
- environment: gen_prob  
- input: results/ligo_simulations_one_seed  
- output: results/generation_prob/simulated_repertoires/all_data.csv  

#### (3) Calculate generation probability of TCRex training data
- Task: Select data for the 8 epitopes with the most data and calculate pgen (remove those with pgen of 0)  
- notebook: notebooks/paper/generation_prob/true_repertoires/tcrex_pgen.ipynb  
- environment: gen_prob  
- input: data/final/all_tcrs.tsv  
- output: results/generation_prob/true_repertoires/tcrex_pgen_data.csv  

#### (4) Compare pgen between clustered TCRs and TCRs that link clusters together
- notebook: Link_clusters_together.ipynb  
- environment: gen_prob  
- True repertoires:  
  - input: results/generation_prob/true_repertoires/tcrex_pgen_data.csv  
  - output: results/generation_prob/true_repertoires/link_clusters  

- Simulated repertoires:  
  - input: results/generation_prob/simulated_repertoires/all_data.csv  
  - output: results/generation_prob/simulated_repertoires/link_clusters  

#### (5) Link_singlets.ipynb 
Environment: gen_prob  
- True repertoires:  
  - input: results/generation_prob/true_repertoires/tcrex_pgen_data.csv  
  - output: results/generation_prob/true_repertoires/link_singlets  

- Simulated repertoires:  
  - input: results/generation_prob/simulated_repertoires/all_data.csv  
  - output: results/generation_prob/simulated_repertoires/link_singlets  


### Large overlap between positive and negative TCR sequences influences model performance

#### (1) Fraction  of pairwis edit-distance distribution
- Task: For every TCR in the positive training data set, calculate the min distance with another positive TCR and a background TCR. Finally, plot the median of the fractions for all TCRs.	
- notebook: background/fraction_per_epitope_pairwise_edit_distance.ipynb 
- environment: raptcr_env
- input: data/parsed/TCRex_data 
- output: results/background/min_distances.tsv 


### Clonal epitope-specific training data sets result in more performant models

#### (1) Influence of diversity on model performance  
- notebook: paper/diversity/diversity_performance_bubbles.ipynb  
- Environment: diversity_metrics_env 
- input:   
  - data/parsed/tcrex_models.csv    
  - results/epitope_specific_clustering/epitope_specific_clusters.tsv  
- output: figures are shown directly in the notebook  


#### (2 )Influence of duplicates on model performance	  
- notebook: paper/duplicates/duplicates.ipynb   
- environment: raptcr_env  
- input:  
  - data/parsed/tcrex_models.csv    
  - data/final/all_tcrs.tsv  
- output: figure is shown directly in the notebook  
