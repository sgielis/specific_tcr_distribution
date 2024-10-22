
# Imports
import os

from Bio import pairwise2
import Levenshtein
import matplotlib.pyplot as plt
import pandas as pd
import plotnine as pn
import networkx as nx
import numpy as np

# Load all packages to calculate generation probability
import olga.load_model as load_model
import olga.generation_probability as pgen
import olga.sequence_generation as seq_gen

olga_dir = 'path_to_your_dir/OLGA-master/olga'

# Define the files for loading in generative model/data
params_file_name = os.path.join(olga_dir,'default_models/human_T_beta/model_params.txt')
marginals_file_name = os.path.join(olga_dir,'default_models/human_T_beta/model_marginals.txt')
V_anchor_pos_file =os.path.join(olga_dir,'default_models/human_T_beta/V_gene_CDR3_anchors.csv')
J_anchor_pos_file = os.path.join(olga_dir,'default_models/human_T_beta/J_gene_CDR3_anchors.csv')

# Load data
genomic_data = load_model.GenomicDataVDJ()
genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
# Load model
generative_model = load_model.GenerativeModelVDJ()
generative_model.load_and_process_igor_model(marginals_file_name)

# Process model/data for pgen computation by instantiating GenerationProbabilityVDJ
pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)


# Calculate the generation prob using the CDR3 beta sequence
def calculate_p_gen(cdr3):
    pgen = pgen_model.compute_aa_CDR3_pgen(cdr3)
    return pgen


# Determine the singlets and the clustered TCRs
def in_cluster(clusters,x):
    if x in clusters['junction_aa'].tolist():
        cls = 'Yes'
    else:
        cls='No'
    return cls


# Cluster with a Levenshtein distance of one
def cluster_sequences(sequences):

    # Create a graph
    G = nx.Graph()

    # Add nodes
    G.add_nodes_from(sequences)

    # Add edges if sequences differ by exactly one Levenshtein distance
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            if Levenshtein.distance(sequences[i], sequences[j]) == 1:
                G.add_edge(sequences[i], sequences[j])

    # Find connected components (clusters)
    clusters = list(nx.connected_components(G))
    clusters = [cluster for cluster in clusters if len(cluster) > 1]

    # Create a DataFrame with sequence and cluster ID
    data = []
    for cluster_id, cluster in enumerate(clusters, 1):
        for sequence in cluster:
            data.append((sequence, cluster_id))
    df = pd.DataFrame(data, columns=["junction_aa", "Cluster"])

    return df


# Get the most closest tcr 
def closest_tcr_singlet(x, sequences):
    all_scores = []
    all_seqs = []
    
    for m in sequences:
        if m != x:
            score = Levenshtein.distance(x, m)
            all_scores.append(score)
            all_seqs.append(m)

    # get max similarity score:
    min_value = min(all_scores)
    min_indices = [index for index, value in enumerate(all_scores) if value == min_value]
    final_tcrs = [all_seqs[pos] for pos in min_indices]
    
    return min_value,final_tcrs

# Get the most closest CDR3 beta not in the same cluster
def most_optimal(x, clusters):
    all_scores = []
    all_seqs = []

    # Get cluster id of x
    selection = clusters[clusters['junction_aa'] == x]
    cluster = selection['Cluster'].tolist()[0]

    # Get all TCRs that do not belong to the cluster of x
    selection2 = clusters[clusters['Cluster']!= cluster]
    subset_cdr = selection2['junction_aa'].tolist()

    # Calculate the Levenshtein distance with every TCR in the other clusters
    for m in subset_cdr:
        score = Levenshtein.distance(x, m)
        all_scores.append(score)
        all_seqs.append(m)

    # Get TCR with minimal Levensthein distance
    min_value = min(all_scores)
    min_indices = [index for index, value in enumerate(all_scores) if value == min_value]
    final_tcrs = [all_seqs[pos] for pos in min_indices]
    
    # Use all TCRs, since we also do this for the singlets
    final_selection = clusters[clusters['junction_aa'].isin(final_tcrs)]
    selected_clusters = final_selection['Cluster'].tolist()

    return min_value,final_tcrs,selected_clusters


# Perform global alignment
def align_sequences(tcr1, tcr2):

    alignments = pairwise2.align.globalms(tcr1, tcr2, 2, -1, -1.1, -0.6)

    # Take the alignment with the higest score
    best_alignment = max(alignments, key=lambda x: x[2])
    aligned_tcr1, aligned_tcr2, score, begin, end = best_alignment

    matches = []
    mismatches = []
    gaps = []

    for i, (a, b) in enumerate(zip(aligned_tcr1, aligned_tcr2)):
        if a == b:
            matches.append((i, a))
        elif a == '-' or b == '-':
            gaps.append((i, a, b))
        else:
            mismatches.append((i, a, b))

    return aligned_tcr1, aligned_tcr2, matches, mismatches, gaps


# Insert the amino acid at the specified position
def insert_amino_acid(tcr_sequence, amino_acid, position):
    if position < 0 or position > len(tcr_sequence):
        raise ValueError("Position is out of bounds")
    updated_sequence = tcr_sequence[:position] + amino_acid + tcr_sequence[position:]
    return updated_sequence


# Replace the amino acid at the specified position
def replace_amino_acid(tcr_sequence, new_amino_acid, position):
    updated_sequence = tcr_sequence[:position] + new_amino_acid + tcr_sequence[position + 1:]
    return updated_sequence


# Create new TCRs that link one TCR (seq1) with the TCRs with the lowest Levenshtein distance (every TCR in match)
def create_new_tcrs(seq1, match):
    new_sequences = []

    for seq2 in match:
         # go from short to long tcr
        if len(seq1)<len(seq2):
            tcr1 = seq1
            tcr2 = seq2
        else:
            tcr1 = seq2
            tcr2 = seq1
        aligned_tcr1, aligned_tcr2, matches, mismatches, gaps = align_sequences(tcr1, tcr2)

        # fix all gaps
        while len(gaps) > 0:
            gap = gaps[0]
            if gap[2] =='-':
                tcr2 = insert_amino_acid(tcr2, gap[1], gap[0])
                new_sequences.append(tcr2)
            else:
                tcr1 = insert_amino_acid(tcr1, gap[2], gap[0])
                new_sequences.append(tcr1)
            aligned_tcr1, aligned_tcr2, matches, mismatches, gaps = align_sequences(tcr1, tcr2)

       # fix mismatches
        while len(mismatches) >1:
            mismatch = mismatches[0]
            tcr1 = replace_amino_acid(tcr1, mismatch[2], mismatch[0])
            new_sequences.append(tcr1)
            aligned_tcr1, aligned_tcr2, matches, mismatches, gaps = align_sequences(tcr1, tcr2)
            
        # When seq1 and seq2 are only differing in gaps, the final added sequence might be identical to one of these
        # To only retain a list of TCRs having a levensthein distance of 1 between each other and these end points,
        # remove the TCRs being identical to these end points
        end_points = match +[seq1]
        new_sequences = [seq for seq in new_sequences if seq not in end_points]

    # Create an unique list as duplicates can arise due to presence of gaps in both sequences (tcr1 and tcr2 are both adapted)
    new_sequences = list(set(new_sequences))

    return new_sequences


def get_size(cluster,sizes):
    selection = sizes[sizes['Cluster']==cluster]
    return selection['junction_aa'].tolist()[0]


def get_cluster(clusters,x):
    selection = clusters[clusters['junction_aa']==x]
    return selection['Cluster'].tolist()[0]


def neighbour_analysis(data, folder):
    
    # Get all epitopes for true repertoires and simulation nrs for simulated repertoires
    epitopes = set(data['epitope'].tolist())
    
    # Perform the neighbour analysis for every epitope/simulation separately
    for epitope in epitopes:

        # Cluster training data for selected epitope
        subset = data[data['epitope']==epitope]
        sequences = subset['junction_aa'].tolist()
        clusters = cluster_sequences(sequences)

        # find number of neighbours with score 1
        clustered_tcrs = clusters['junction_aa'].tolist()
        df = pd.DataFrame()

        for tcr1 in clustered_tcrs:
            counter = 0
            for tcr2 in clustered_tcrs:
                score = Levenshtein.distance(tcr1, tcr2)
                if score ==1: 
                    counter +=1

            new = pd.DataFrame({'TCR':[tcr1],'counts':[counter]})
            df = pd.concat([df,new], axis=0)

        df['p_gen'] = df['TCR'].apply(lambda x: calculate_p_gen(x))
        df['log_p'] = df['p_gen'].apply(lambda x: -np.log10(x))

        # Boxplot of log pgen versus number of neighbours
        plot = (pn.ggplot(df, pn.aes(x='counts', y='log_p', group='counts')) +
                pn.geom_boxplot()+
                pn.scale_x_continuous(breaks=[1]+ list(range(5, df['counts'].max() + 2, 5)))+ #+2 instead of +1 since one epitope ends at 4
                pn.ggtitle(epitope)+
                pn.xlab('Number of Neighbours')+
                pn.ylab('-log(generation probability)')+
                pn.theme(figure_size=(10, 8)))

        # Save the plot to a file
        plot.save('./results/pgen/'+ folder +'/nr_neighbours/'+ epitope +'.png', dpi=600)

        # Boxplot of log pgen versus cluster size
        sizes = clusters.groupby('Cluster').count().reset_index()
        df['cluster'] = df['TCR'].apply(lambda x: get_cluster(clusters,x))
        df['size'] = df['cluster'].apply(lambda x: get_size(x,sizes))

        plot = (pn.ggplot(df, pn.aes(x='size', y='log_p', group='size')) +
                pn.geom_boxplot()+
                pn.xlab('Cluster size')+
                pn.ylab('-log(generation probability)')+
                pn.ggtitle(epitope))

        # Save the plot to a file
        plot.save('./results/pgen/'+ folder +'/cluster_size/'+ epitope +'.png', dpi=600)


def connect_clusters(data, folder):

    # Get all epitopes for true repertoires and simulation nrs for simulated repertoires
    epitopes = set(data['epitope'].tolist())

    for epitope in epitopes:

        # Select training data for selected epitope
        subset = data[data['epitope']==epitope]

        # Cluster the data using a levenshtein distance of one
        sequences = subset['junction_aa'].tolist()
        clusters = cluster_sequences(sequences)

        # For every TCR in a cluster, find the closest TCR in another cluster
        clusters[['score','match','matched_clusters']] = clusters['junction_aa'].apply(lambda x: pd.Series(most_optimal(x,clusters)))

        # For every original junction_aa and its match, calculate all needed connecting TCRs
        clusters['connecting_tcrs'] = clusters.apply(lambda row: create_new_tcrs(row['junction_aa'], row['match']), axis=1)

        # Collect all original junction_aa and the new connecting tcrs in one df
        clustered = clusters['junction_aa'].tolist()
        clustered = pd.DataFrame({'junction_aa':clustered ,'origin':['clustered']*len(clustered )})

        new = clusters['connecting_tcrs'].tolist()
        new = [item for sublist in new for item in sublist]
        new = pd.DataFrame({'junction_aa':new,'origin':['connecting_TCRs']*len(new)})

        all_sequences = pd.concat([clustered,new])
        all_sequences = all_sequences.drop_duplicates(subset='junction_aa', keep='first')

        # Calculate p gen for every sequence using only the CDR3 beta sequence
        all_sequences['p_gen'] = all_sequences['junction_aa'].apply(lambda x: calculate_p_gen(x))
        all_sequences['logp'] = all_sequences['p_gen'].apply(lambda x: -np.log10(x))

        # Create the plot
        plot = (pn.ggplot(all_sequences, pn.aes(x='origin', y='logp', group='origin')) +
                pn.geom_boxplot()+
                pn.ggtitle(epitope)+
                pn.xlab('Type of TCR')+
                pn.ylab('-log10(generation probability)')+
                pn.theme(figure_size=(10, 8)))

        # Save the plot to a file
        plot.save('./results/pgen/'+folder+'/link_clusters/'+ epitope +'.png', dpi=600)


def connect_singlets(data,folder):
    
    # Get all epitopes for true repertoires and simulation nrs for simulated repertoires
    epitopes = set(data['epitope'].tolist())

    # Perform the analysis for every epitope/simulation separately
    for epitope in epitopes:

        # Select training data for selected epitope
        print('Performing analysis for epitope: ', epitope)
        subset = data[data['epitope']==epitope]

        # Cluster the data using a Levenshtein distance of one
        sequences = subset['junction_aa'].tolist()
        clusters = cluster_sequences(sequences)

        # Get all singlets
        subset['clustered'] = subset['junction_aa'].apply(lambda x: in_cluster(clusters,x))
        singlets = subset[subset['clustered'] == 'No']

        # Get the TCR with the shortest Levenshtein distance for each singlet
        singlets[['score', 'match']] = singlets['junction_aa'].apply(lambda x: pd.Series(closest_tcr_singlet(x,sequences)))

        # Find new TCRs connecting the singlets with clusters or other singlets
        singlets['connecting_tcrs'] = singlets.apply(lambda row: create_new_tcrs(row['junction_aa'], row['match']), axis=1)

        # Create one df with clustered data, singlets and new connecting tcrs
        clustered = subset[subset['clustered'] == 'Yes']
        clustered = clustered['junction_aa'].tolist()
        clustered = pd.DataFrame({'junction_aa':clustered ,'origin':['clustered']*len(clustered )})

        original = singlets['junction_aa'].tolist()
        original = pd.DataFrame({'junction_aa':original,'origin':['singlets']*len(original)})

        connection = singlets['connecting_tcrs'].tolist()
        connection = [item for sublist in connection for item in sublist]
        connection = list(set(connection))
        connection = pd.DataFrame({'junction_aa':connection,'origin':['connecting TCRs']*len(connection)})

        all_sequences = pd.concat([clustered,original,connection])
        all_sequences = all_sequences.drop_duplicates(subset='junction_aa', keep='first')

        # Compare pgen values between clustered, singlets and connecting TCRs
        all_sequences['p_gen'] = all_sequences['junction_aa'].apply(lambda x: calculate_p_gen(x))
        all_sequences['logp'] = all_sequences['p_gen'].apply(lambda x: -np.log10(x))

        # Create the plot
        plot = (pn.ggplot(all_sequences, pn.aes(x='origin', y='logp', group='origin')) +
                pn.geom_boxplot()+
                pn.ggtitle(epitope)+
                pn.xlab('Type of TCR')+
                pn.ylab('-log10(generation probability)')+
                pn.theme(figure_size=(10, 8)))

         # Save the plot to a file
        plot.save('./results/pgen/'+folder+'/link_singlets/'+ epitope +'.png', dpi=300)


