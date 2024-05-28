#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 28 August 2023 at 15:06
Modified on 16 May 2024

Author: Aleksandra Panfilova
"""

import functools
import os
import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import textwrap
from Bio.Seq import Seq
import re
import statistics
import math
from scipy import stats
from itertools import combinations

"""
Function that takes in the .txt file with sequence counts and filters out
sequences with less counts than 'threshold' variable value, as well as sequences
with unexpected lengths (for example, forward or reverse adapters untrimmed)
"""
def read_DNA_counts(path_to_txt, threshold, min_len, max_len, rep, print_yn=False):
    count = list() # 'count' column of the future dataframe
    sequence = list() # 'sequence' column of the future dataframe
    f = open(path_to_txt, "r")
    lines = f.readlines()
    # Go through the .txt line by line and get a count and a sequence
    for line in lines: 
        # Finds the numbers in the string, returns a list (still in str format)
        numbers = re.findall(r'\d+', line) 
        res = list(map(int, numbers)) # Turns str into int
        seq = " ".join(re.findall("[a-zA-Z]+", line)) # Gets the sequence as a str
        # Write only the reads with sufficient counts and proper lenght
        if (res[0] >= threshold)&(len(seq)>=min_len)&(len(seq)<=max_len): 
            count.append(res[0])
            sequence.append(seq)
    # Create a DataFrame to store the sequences and their counts
    DNA_df = pd.DataFrame() 
    DNA_df['sequence'] = sequence
    DNA_df['count'] = count
    # To keep the information on replicates through the downstream analysis, 
    # put the replicate number in the column name
    DNA_df.rename(columns={"count": ("count_" + rep)}, inplace=True)
    # Throughout the code, if passed True, print_yn returns the detailed 
    # information on the sequencing data
    if print_yn == True: 
        print("\nInitial N of variants: ", len(lines))
        print("After filtering: ", len(DNA_df))
    return DNA_df

"""
Function that takes in sequences, checks them for being the correct
length, in-frame, and without additional mismatches, and classifies 
them as synonymous variant, stop-codon variant, insertion or deletion.
For indels, mismatches on the DNA level are allowed, since
this code is focused on protein-level changes, and they are most
probably the sequencing errors. For stop and synonymous codons, mismatches
are not allowed. For the rare cases when the direction of the reads is not 
forward, user can specify 'direction' variable (from default 'fwd' to 'rev')
in the 'sequence_data_dict'.
Residue position numeration starts with 1.
"""
def var_call(DNA_df, wt_dna, rep, direction='fwd', print_yn = False):
    # Turn the str into the biopython Seq type and make sure it's in uppercase 
    # (important for the aligment)
    wt_dna = Seq(wt_dna.upper()) 
    wt_aa = wt_dna.translate() # Protein sequence from DNA sequence
    old_res = ' ' # Previous residue
    old_codon = ' ' # Previous codon
    sequential_positions = list()
    
    """
    Identify positions with sequential identical amino acids coded by 
    non-identical codons. This is important because most of the analysis is 
    made on the protein level, and deletion of one of the identical residues 
    going one after another will be the same on the protein level, 
    but not on the DNA level. In this case we want to check the sequence
    on the DNA level to assign the correct position to the variant.
    Position variable is set to 0 for the wild type and some of the sequences we 
    want to filter out: synonymous (to wild type) variants with more than one 
    codon mutated, stop-codon variants with additional DNA mismatches, 
    out-of-frame sequences, protein-level mismatches.
    """
    for n_r, res in enumerate(wt_aa):
        # Find the codon corresponding to the residue
        codon = wt_dna[((n_r)*3):(((n_r)*3)+3)] 
        if res == old_res: # Are these two sequential residues?
            if codon != old_codon: # Are they coded by non-identical codons?
                sequential_positions.append(n_r)
        old_res = res
        old_codon = codon
        
    syn_count = 0 # Counter to sum all of the synonymous variants for reporting
    mut_type_list = list() # 'mut_type' column to add to the df
    pos_list = list() # 'pos' column to add to the df
    pos = 0 # Variable for mutation position tracking
    for i, row in DNA_df.iterrows():
        # Get the Seq object from the str variant sequence,
        # read direction considered.
        if direction == 'rev': 
            var_dna = Seq(row['sequence']).reverse_complement()
        else:
            var_dna = Seq(row['sequence'])
        
        # This checks for DNA length and classifies variants as out-of-frame
        if len(var_dna)%3 == 0:
            # If correct length, get protein sequence
            var_aa = var_dna.translate() 
            # In case this is a stop-variant, get the protein sequence which 
            # terminates at the stop codon. By default, Seq().translate() 
            # gives out the full sequence with * in place of the stop-codon.
            var_aa_stop = var_dna.translate(to_stop=True)
        else:
            mut_type = 'out_of_frame' 
            pos = 0
            mut_type_list.append(mut_type)
            pos_list.append(pos)
            continue
        
        # This statement allows to fish out non-indels
        if len(var_aa)==len(wt_aa):
            if var_dna == wt_dna: # We found the wild type
                mut_type = 'wt'
                pos = 0
                # Get the WT count for detailed reporting
                wt_count = row[('count_'+rep)] 
            elif (var_dna != wt_dna)&(var_aa==wt_aa): # Synonymous
                bp_pos = 0 # Position in the DNA seq
                mismatch = 0 # Variable for counting mismatching codons
                res_pos = 0 # Residue position
                for bp_wt, bp_var in zip(wt_dna, var_dna):
                    if bp_wt != bp_var: # syn to wt DNA mismatch
                        if res_pos == 0: # Is this the first mismatch?
                            # Residue position from DNA position
                            # int() cuts off the decimals
                            res_pos = int(bp_pos/3)
                            mismatch += 1
                        else: # next mismach
                            # Is it in the same codon?
                            # If not, it's a multimutant
                            if int(bp_pos/3) != res_pos:
                                mismatch += 1
                    bp_pos += 1
                if mismatch == 1: # Only one codon is replaced
                    mut_type = 'syn'
                    pos = res_pos + 1
                else: # Multiple synonymous substitutions
                    mut_type = 'multi_syn'
                    pos = 0
                # Get the syn count for detailed reporting
                syn_count += row[('count_'+rep)]
            elif '*' in var_aa: # Is this a stop codon?
                pos = var_aa.find('*') + 1 
                # Get the variant sequence with the stop codon cut out
                no_stop_var = var_dna[0:((pos-1)*3)]+var_dna[(pos*3):]
                # WT with mutated codon cut out
                no_stop_wt = wt_dna[0:((pos-1)*3)]+wt_dna[(pos*3):]
                # Is everything apart from the one mutated codon identical?
                if no_stop_var==no_stop_wt:
                    mut_type = 'stop'
                else: # If not, there is a mismatch on top of the stop mutation
                    pos = 0
                    mut_type = 'stop_mismatch' 
            # If none of the above is true, it's a protein-level mismatch variant:
            # same length as WT, but not synonymous or stop
            else:  
                mut_type = 'mismatch'
                pos = 0
        #Is this an indel? +-1 residue, and not a stop variant
        elif (((len(var_aa)==(len(wt_aa)-1))|(len(var_aa)==(len(wt_aa)+1)))
              &('*' not in var_aa)):
            # Align the variant to the WT to get the aligment score and 
            # the gapped output.
            align_var = pairwise2.align.globalms(wt_aa, var_aa, 1, -1, -1, -.5)
            # If gap is in seqA (WT), it's an insertion.
            # The scores in the aligment are giving 1 for a matching residue,
            # -1 for mismatching, -1 for a gap opening. So, a single residue  
            # insertion (protein sequence level) should have a score of
            # L_seq * 1 - 1 (number of matches = number of residues,
            # minus one gap opening).
            if ('-' in align_var[0].seqA)&(align_var[0].score==(len(wt_aa)-1)):
                # Position of the residue after which the insertion was made
                ins_pos = align_var[0].seqA.find('-') 
                ins_dna_start = ins_pos*3 # First bp of the insertion codon
                # Is the inserted residue Gly, and is it the codon we used 
                # for all of the insertions?
                if ((align_var[0].seqB[ins_pos]=='G')&
                    (var_dna[(ins_dna_start):ins_dna_start+3]=='GGT')):
                    mut_type = 'ins'
                else: # If not, it's not the insertion we designed
                    mut_type = 'wrong_ins'
                pos = ins_pos
            # If gap is in seqB (variant), it's a deletion.
            # A single deletion (protein sequence level) should have a score of
            # L_seq * 1 - 2 (number of matches = number of residues - 1,
            # minus one gap opening)
            elif ('-' in align_var[0].seqB)&(align_var[0].score==(len(wt_aa)-2)):
                # Position of the deleted residue 
                # (+1 because Python counts from 0)
                pos = align_var[0].seqB.find('-')+1
                mut_type = 'del'
                if pos in sequential_positions:
                    # If this is one of the positions where the same residues go 
                    # one after another, but are encoded by different codons, 
                    # we want to make sure that the position we assign is 
                    # according to the DNA aligment
                    align_dna = pairwise2.align.globalms(wt_dna, 
                                                        var_dna, 1, -1, -1, -.5)
                    for a in align_dna:
                        if len(a.seqB.split('---')[0])%3 == 0:
                            pos = int(a.seqB.find('-')/3)+1
            # If the score is wrong, there are additional mismatches
            # in the protein sequence
            else: 
                pos = 0
                mut_type = 'indel_mismatch'
        # Everything else we did not design to be in the DNA library.
        else: 
            mut_type = 'something_else'
            pos = 0

        mut_type_list.append(mut_type)
        pos_list.append(pos)
        pos = 0 # Reset the position
        
    if print_yn == True:
        print('WT count: ', wt_count)
        print('Synonymous count: ', syn_count)
    DNA_df['mut_type'] = mut_type_list
    DNA_df['pos'] = pos_list
    DNA_df.sort_values(by=['pos'], ignore_index=True, inplace=True)
    # Order the columns
    DNA_df = DNA_df[['sequence', 'pos', 'mut_type', ('count_'+rep)]]
    return DNA_df

"""
Main analysis function that handles data in the following order: reads data per 
replicate, then merges together replicates from different conditions and 
calculates for every replicate a ratio of variant counts to the wild type or 
synonymous counts depending on the selected normalisation type.
This function outputs data in various stages of analysis to simplify
troubleshooting, but tile_shared_averaged is what is used to produce heatmaps 
and can be considered the main output.

It normalizes and scales the data according to the specified method

sequence_data_dict (one of the inputs) is a dictionary constructed by user,
which contains different data parameters: number and names of reps, conditions
and tiles, as well as wilt type sequences for tiles ans for the full protein,
min and max read lenghts, values for desired mut_types, and read directions.
"""
def analysis(path_to_folder, sequence_data_dict, normalisation_scaling_type, 
             print_yn = False, threshold = 35):
    # A list with unconcateneted data frames for each tile
    tile_list = list() 
    # Dictionary of dataframes with filtered out (not a wt, syn, stop, or indel) 
    # reads for each replicate. Has one dataframe per key, key name structured
    # as 'tile_condition_replicate'
    filtered_out = dict()
    # Data frame with wild type sequence mean counts for all tiles,
    # all conditions
    wt_df = pd.DataFrame()
    # Dictionary with filtered sequence counts for every replica. 
    # Structured as following: rep_counts[tile][condition][N],
    # where N is the index of replicate (0, 1, 2), and both tile and condition
    # are subdictionaries
    rep_counts = dict()
    # Gets keys of the reps_cond_tile dictionary in sequence_data_dict:
    # those keys are tiles
    for tile in sorted(sequence_data_dict['reps_cond_tile'].keys()):
        rep_counts[tile] = dict()
        # Gets a wt sequence for this tile
        wt_tile = sequence_data_dict['wt'][tile] 
        # Local aligment -- we align short tile sequence vs the long 
        # full-length sequence to get the starting residue of the tile.
        # We are then using it to have a continuous numeration through the
        # tiles, as var_call() function numerates according to the single tile
        # sequence, starting from 1.
        # May cause problems if you have tiles with repeats on 
        # both ends of the sequence?
        align = pairwise2.align.localms(sequence_data_dict['wt']['full'], 
                                        wt_tile,  1, -1, -1, -.5)
        # Get a position of the first residue of the tile in the full-
        # length sequence from the first aligned DNA position.
        # This number is going to be added to all position numbers in 
        # the tile dataframe. int() cuts off the decimals.
        tile_modifier = int(align[0].start/3)

        # Build a 'scaffold' initial dataframe, which contains all of the 
        # expected variants: 4 mutations at every position.
        # This allows to easily see what are the missing variants.
        positions_list = list()
        mutation_list = list()
        aa_seq_list = list()
        # Designed mutation classes
        mut_set = [mut for mut in sequence_data_dict['values']
                   if mut != 'wt']
        # Go through the indexed wild type tile protein sequence
        for i, aa in enumerate(list(Seq(wt_tile.upper()).translate())):
            # Put position number and residue type into the column 4 times 
            # (number of mutation classes), so that each mut_type has pos
            # and aa identifiers
            positions_list.extend([i+1+tile_modifier]*4)
            aa_seq_list.extend(aa*4)
            mutation_list.extend(mut_set)
            
        scaffold_df = pd.DataFrame({'tile': [tile]*len(positions_list),
                                   'pos': positions_list,
                                   'aa': aa_seq_list,
                                   'mut_type': mutation_list})

        condition_all_df_list = list()
        for condition in sequence_data_dict['reps_cond_tile'][tile].keys():
            # List to collect replicate dataframes to merge them later
            replicates_all_df_list = list() 
            for rep in sequence_data_dict['reps_cond_tile'][tile][condition]:
                if print_yn == True: # To separate between tile info blocks
                    print('\n'+tile.upper()+' '+condition+' '+rep)
                
                # Counstruct path to the .txt count file
                file_path = os.path.join(path_to_folder,
                    tile, (tile + '-' + condition + '-' + rep + '.txt'))
                df = var_call(read_DNA_counts(file_path, threshold, 
                    sequence_data_dict['min_max_len'][tile][0], 
                    sequence_data_dict['min_max_len'][tile][1], rep, print_yn),
                    wt_tile, rep, sequence_data_dict['direction'][tile], 
                    print_yn)
                # Drop everything that is not a designed mutation
                vars_clean = df.loc[df.mut_type.isin(
                    sequence_data_dict['values']) == True
                                   ].reset_index(drop=True)
                # Keep everything that is not a designed mutation
                filtered_out[(f'{tile}_{condition}_{rep}')] = df.loc[
                    df.mut_type.isin(sequence_data_dict['values']) == False
                                   ].reset_index(drop=True)
                # Replace the DNA sequence with protein sequence to sum together 
                # varriants that are synonymous on the protein level, but 
                # different on the DNA level (only for indels).
                vars_clean['sequence'] = vars_clean['sequence'].apply(
                    lambda x: str(Seq(x).translate()))
                clean_sum = vars_clean.groupby(['sequence', 'pos',
                                    'mut_type'], dropna=False, as_index=False,
                                    group_keys = False).sum().sort_values(
                                    by=['pos'], ignore_index=True)
                
                # Calculate the per replicate ratio of count_var/count_syn
                # using the sum of all synonymous variants in the replicate. 
                # Depending on the chosen normalisation type, ratios are going 
                # to be calculated to wild type (wt), sum of all synonymous 
                # variants in the tile (syn_sum), or for synonymous score of 
                # each position (syn_pos)
                if normalisation_scaling_type == 'syn_sum':
                    clean_sum[f'ratio_{rep}'] = clean_sum[f'count_{rep}'
                            ].values/float(clean_sum.loc[clean_sum['mut_type'] 
                                == 'syn'][f'count_{rep}'].sum())
                elif normalisation_scaling_type == 'wt':
                    clean_sum[f'ratio_{rep}'] = clean_sum[f'count_{rep}'
                            ].values/float(clean_sum.loc[clean_sum['mut_type'] 
                                == 'wt'][f'count_{rep}'])
                #elif normalisation_scaling_type == 'syn_pos':
                # Do we even need it?
                else:
                    print('Wrong normalisation type!')
                    print('It should be one of the 3: syn_sum, syn_pos, or wt')
                    return
                # Add replicate dataframe to the list to merge later
                replicates_all_df_list.append(clean_sum)
                
            # Now we enter the condition level -- we've collected all of the
            # replicate dataframes for condition in one list.
            
            # # Save the list with counts and ratios to the 'raw' counts dictionary
            # rep_counts[tile][condition] = replicates_all_df_list
            # Merge data for all replicates based on position, type of mutation,
            # and *sequence*. This means that variants that are only present 
            # in some but not all replicates will have NaN values for some
            # of the columns. We use this to drop all variants like that.
            all_rep_df = functools.reduce(lambda x, y: x.merge(y, on=['mut_type', 
                           'pos', 'sequence'], how='left'), replicates_all_df_list)
            clean_df = all_rep_df.dropna().reset_index(drop=True)
            # Save the list with counts and ratios to the 'raw' counts dictionary
            ratio_cols = [c for c in clean_df.columns if c.startswith('ratio')]
            rep_counts[tile][condition] = clean_df.drop(columns=ratio_cols).copy()
            # Calculates and prints Spearman's correlation between the replicates.
            if print_yn == True:
                print('\n')
                # Form a list of count columns
                reps = [x for x in clean_df.columns if x.startswith('count')]
                # Get all possible combinations of the values in the list
                # (order is not important)
                combos = list(combinations(reps, 2))
                for i in range(len(combos)):
                    # Get the first (x) and the second (y) value form the 
                    # combination
                    x = combos[i][0].split('_')[1]
                    y = combos[i][1].split('_')[1]
                    print(f'Pearson {x} to {y}: ', 
                          stats.spearmanr(clean_df[combos[i][0]],
                                          clean_df[combos[i][1]])[0])

            #STATISTICAL ANALYSIS STARTS
            # Form a list of ratio columns to calculate mean ratio and a standard
            # error between reps
            ratio_list = [r for r in clean_df.columns if r.startswith('ratio')]
            clean_df['mean_ratio_'+condition] = clean_df.loc[:,ratio_list].mean(
                axis=1).round(6)
            clean_df['std_ratio_'+condition] = clean_df.loc[:, ratio_list].std(
                axis=1).round(6)
            # Drop the count columns and per replicate ratio columns, then save 
            # the resulting condition dataframe in a list for merging with other
            # conditions
            condition_all_df_list.append(clean_df[['sequence', 'pos', 'mut_type',
                'mean_ratio_'+condition, 'std_ratio_'+condition]])
        # Merge dataframes for all the conditions in one full tile dataframe
        tile_df = functools.reduce(lambda x, y: x.merge(y, on=['mut_type', 'pos',
                    'sequence'], how='left'), condition_all_df_list)
        # Remember the tile modifier to keep residue numeration consistent?
        tile_df['pos'] = tile_df['pos']+tile_modifier
        # Save the row with the wild type data in a separate dataframe
        wt_row = tile_df[tile_df.mut_type == 'wt']
        wt_df = pd.concat([wt_df, wt_row.join(pd.DataFrame({'tile': [tile]}),
                            how='left')])
        # Drop the wild type row
        tile_all = tile_df.drop(tile_df[tile_df.pos == 'wt'].index)
        # Finally, merge the tile dataframe with ratios per condition to the
        # scaffold dataframe to get structured ordered data
        tile_all = scaffold_df.merge(tile_all, on=['mut_type', 'pos'], 
                                     how='left')
        
        # Now, we need to calculate log ratio scores. The formula is 
        # ln(condition_ratio/control_ratio)
        
        # The first condition in the list is the control condition:
        # this is specified by the user when constructing the sequence_data_dict()
        control = list(sequence_data_dict['reps_cond_tile'][tile].keys())[0]
        # Here we get the list of names of the columns that contain values of
        # ratios and stds
        ratios_sele_col = [col for col in tile_all if (col.startswith('mean_'))&
                    ~col.endswith(control)]
        error_sele_col = [col for col in tile_all if (col.startswith('std_'))&
                    ~col.endswith(control)]
        # Get the score and error column names from the ratios and stds.
        # These are the names for new columns we are going to assign calculated
        # values to.
        score_col = [f'score_{c.split("_")[-1]}' for c in ratios_sele_col]
        error_col = [f'err_{c.split("_")[-1]}' for c in error_sele_col]
        # Divide ratios for the selection by ratios for the control.
        # pd.DataFrame.div() is used to divide all 4 ratio columns
        # by the same column -- control ratio column.
        # Output is the np.array (.values converts dataframe into the np.array)
        # with the size of (160, 4)
        ratio_of_ratio = tile_all[ratios_sele_col].div(tile_all
                                          ['mean_ratio_'+control], axis=0).values
        # See the formulas in the manuscript. Basically,
        # we separate one long formula into sub-operations for clarity.
        # The full formula for propagation of uncertainty of the standard
        # deviation for the division C = A/B with stds sA and sB is:
        # sC = abs(C)*sqrt((sA/A)^2+(sB/B)^2).
        # First, we calculate the (sB/B)^2, or, in our case, control component.
        # The expected output is the np.array of the size (160, ), and this
        # is an element-wise division
        std_to_mean_ctrl = (tile_all['std_ratio_'+control].values/tile_all
                             ['mean_ratio_'+control].values)**2
        # Next, we get the (160, 4) np.array insteas of (160, ) for the (sB/B)^2.
        # We need this to add these values to all 4 selection conditions. For all 
        # 4 of them, control values are the same, of course. The output is a 
        # np.array where all 4 columns are identical -- (sB/B)^2 operation.
        std_to_mean_ctrlx4 = np.repeat(std_to_mean_ctrl2[:, np.newaxis], 4,
                                         axis=1)
        # Now we calculate the (sA/A)^2, or, in our case, condition component.
        # Same logic as the control, (160, 4) output.
        std_to_mean_sele = (tile_all[error_sele_col].values/tile_all
                             [ratios_sele_col].values)**2
        # Implementation of the full formula
        error_div = abs(ratio_of_ratio) * np.sqrt(std_to_mean_sele
                                                  +std_to_mean_ctrlx4)
        # Getting log ratio score as per formula in the beggining of this
        # section. Natural logarithm of condition to control ratio.
        tile_all[score_col] = np.log(ratio_of_ratio)
        # Another error propagation step, this time using formula:
        # If we have y = ln(x), where sx is the std for x,
        # sy = sx/x.
        tile_all[error_col] = error_div/ratio_of_ratio
        
        """
        The last step of the score calculations, rescaling the scores from 0 to
        1, with 0 at mean stop-codon score per condition, and 1 at mean synonymous
        (which is assumed to be 0 before rescaling, because this is what we 
        normalise to). Synonimous variant count divided by sum of all varaint 
        counts should have arounf the same value across all conditions,
        so ln(1)=0.
        Formula used for rescaling is the standard one:
        ((max_value_new - min_value_new) * (data_point - min_value_original_data))/
        (max_value_original_data - min_value_original_data) + min_value_new. 
        Since our new min and max are 0 and 1, this simplifies into 
        (data_point - min_value_original_data)/(-min_value_original_data)
        """
        min_params = tile_all.loc[tile_all['mut_type'] == 'stop'][score_col].mean()
        tile_all[score_col] = ((tile_all[score_col]-min_params)/
                               (-min_params)).round(3)
        # To keep error's scale consistnt with the score, we divide it with the 
        # max-min difference of the origianl data
        tile_all[error_col] = (tile_all[error_col]/(-min_params)).round(3)
        """
        Next part is for copying the values of the missing variants for mutations 
        that have one sequence for two or more positions. For example, there can only
        be one sequence for insertion of Gly before and after a Gly in the sequence, 
        because this is going to be two Gly in a row in both cases. Same goes for
        deletions when we have two (or more) identical codons one after another.
        So, this part of the code goes through the dataframe, finding these cases
        and copying score, sequence and error values to the positions which were
        missing because of that.
        """
        # Get the variants with missing sequences (remeber, we have a scaffold 
        # dataframe that has all possible variants for the set-up). Most of 
        # those are impossiblesynonymous variants at positions with Met and Trp.
        # Some were really absent from the sequencing data. And a handful are
        # same on the DNA sequence level, as described above.
        nans = tile_all[(tile_all['sequence'].isnull())].copy()
        num_rows = len(nans)
        # Get the joined list of score and error columns, plus sequence. 
        score_err_cols = score_col + error_col
        score_err_cols.insert(0, 'sequence')

        # Going through the dataframe backwards beacause of the 3-in-a-row cases -
        # because of the way pairwise2 works for deletions and the fact that we 
        # always pick the first aligment from the list, the varinat is classified
        # as the first residue of the set, so we are copying from first to last.
        for index in range(num_rows - 1, -1, -1): 
            row = nans.iloc[index]
            if row['mut_type'] == 'ins':
                # Find the row in the dataframe that is an insertion after Gly
                # at the position after the one of the missing variant we are 
                # looking at. For insertions (again, because of how pairwise2
                # works), the variant in a case like this is classified as the
                # last one of the set.
                score_ins = tile_all.loc[(tile_all.mut_type == 'ins')&(tile_all.aa 
                         == 'G')&(tile_all.pos == (row['pos']+1))
                                        ][score_err_cols].values
                if print_yn == True:
                        print('Found a missing insertion before G')
                        print('Position: ', row['pos']+1)
                        print(row)
                        print(score_ins)
                        
                if len(score_ins)==1:
                    # Get the index of missing sequecne variant for the row 
                    # editing
                    ind = tile_all.loc[(tile_all.mut_type == 'ins')&
                                       (tile_all.pos == row['pos'])].index
                    tile_all.loc[ind, score_err_cols]  = score_ins
                else:
                    print('Do you have long Gly repeats?')

            # Same idea for the deletions, but different direction
            elif row['mut_type'] == 'del':
                score_del = tile_all.loc[(tile_all.mut_type == 'del')&(tile_all.aa 
                        == row['aa'])&(tile_all.pos == (row['pos']-1)
                                      )][score_err_cols].values
                if len(score_del) == 1:
                    ind = tile_all.loc[(tile_all.mut_type == 'del')&(tile_all.pos 
                                    == row['pos'])].index
                    tile_all.loc[ind, score_err_cols] = score_del

        # Finally, add the tile to the list of tiles
        tile_list.append(tile_all.reset_index(drop=True))
        
    # After iterating through all of the tiles, concat tile dfs together
    df_full = pd.concat(tile_list, ignore_index=True)

    # Get a dataframe, where scores for the overlapping tile regions are averaged,
    # and ratio columns dropped, leaving just scores and errors
    tile_shared_averaged = df_full.groupby(['pos', 'mut_type', 'aa'], dropna=False, 
            as_index=False, group_keys = False)[[col for col in df_full if 
                        col.startswith('score_')|col.startswith('err_')]].mean()
    
    return tile_shared_averaged, df_full, tile_list, wt_df, rep_counts, filtered_out