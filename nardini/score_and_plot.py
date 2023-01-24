# -*- coding: utf-8 -*-
import os
import math
import string
import random
import numpy as np
import pandas as pd
from pprint import pprint
from datetime import datetime
from zipfile import ZipFile
from tabulate import tabulate
from nardini.constants import MAPPING_8x8, MAPPING_9x9, LABELS_8x8, LABELS_9x9
from nardini.core import count_residues_in_sequence
from nardini.core import get_org_seq_vals, get_scramble_seqs_vals
from nardini.plotting import plot_zscore_matrix


def export_analysis(to_export):
    random_string = str(''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)))
    zip_filename = f'nardini-data-{random_string}.zip'
    zip_filepath = os.path.join(zip_filename)
    zfile = ZipFile(zip_filepath, 'w')

    # create zip file
    d = dict()
    d['ID'] = list()
    d['original_seq'] = list()
    d['most_similar_seq'] = list()
    d['sum_abs_zscore_original_seq'] = list()
    d['sum_abs_zscore_scrambled_seq'] = list()
    for seq_id in to_export:
        oseq, sseq, index, zplot, splot, zm, sm = to_export[seq_id]
        z_abs_sum = np.sum(abs(zm))
        s_abs_sum = np.sum(abs(sm))

        d['ID'].append(seq_id)
        d['original_seq'].append(oseq)
        d['most_similar_seq'].append(sseq)
        d['sum_abs_zscore_original_seq'].append(z_abs_sum)
        d['sum_abs_zscore_scrambled_seq'].append(s_abs_sum)

        zfile.write(zplot, os.path.basename(zplot))
        zfile.write(splot, os.path.basename(splot))

        labels = LABELS_8x8
        if zm.shape[0] == 9:
            labels = LABELS_9x9
        z_content = tabulate(zm, labels, tablefmt='plain')
        s_content = tabulate(sm, labels, tablefmt='plain')
        zfile.writestr(f'zscore-original-sequence-{seq_id}.tsv', z_content)
        zfile.writestr(f'zscore-scrambled-sequence-{seq_id}.tsv', s_content)

    df = pd.DataFrame(d)
    tsv_content = tabulate(df.values.tolist(), list(df.columns), tablefmt="plain")
    zfile.writestr('sequences.tsv', tsv_content)
    zfile.close()
    print(f'Analysis results saved to: "{zip_filename}"')


def print_amino_acid_groupings_header(typeall):
    mapping = None
    if len(typeall) == len(MAPPING_8x8):
        mapping = MAPPING_8x8

    elif len(typeall) == len(MAPPING_9x9):
        mapping = MAPPING_9x9

    length = max(map(len, mapping))
    padding = 5
    space = length + padding
    for group_name, group in mapping.items():
        name = '\t\t{0:{space}}'.format(group_name, space=space)
        group_ouput = f'{sorted(group)}'
        print(name + group_ouput)
    print()


def calculate_zscore_and_plot(orthseqs, typeall, num_seqs, random_seed):
    olen        = len(orthseqs)
    tlen        = len(typeall)
    zvecdb      = np.zeros((olen, tlen**2))
    zvecdbscr   = np.zeros((olen, tlen**2))
    countseqs   = -1
    start       = datetime.now()

    all_zscore_plots = list()
    all_zscore_scrambled_plots = list()
    to_export = dict()
    print(f'Number of amino acid groupings used for Nardini z-score analysis: {tlen}', end='\n\n')
    print(f'Amino acid groupings used for Nardini z-score analysis: ', end='\n\n')
    print_amino_acid_groupings_header(typeall)

    for seq_record in orthseqs:
        seq_name = str(seq_record.id)
        myseq = str(seq_record.seq)
        fracsall = list()
        countseqs = countseqs + 1

        print('=' * 80, end='\n\n')
        print(f'[ FASTA {seq_name} : {myseq} ]')
        print(f'[ SEQ: {seq_name} | START ] Beginning analysis of FASTA sequence: "{seq_name}"...', end='\n\n')

        for type1 in typeall:
            count = count_residues_in_sequence(myseq, type1)
            fraction = count / len(myseq)
            fracsall.append(fraction)

        print(f'[ SEQ: {seq_name} | 1 / 8 ] Analyzing original sequence...')
        myarr = get_org_seq_vals(myseq, typeall, fracsall)
        print(f'[ SEQ: {seq_name} | 2 / 8 ] Analysis of original sequence complete.')

        # Returns mean of scrambles, std of scrambles, all values in a number of scramble x tlen**2 list, and all scramble sequences
        print(f'[ SEQ: {seq_name} | 3 / 8 ] Performing analysis of {num_seqs} scrambled sequences: mean, stddev, etc...')
        alpha, amean, avar, allscrvals, allscrseqs = get_scramble_seqs_vals(myseq, num_seqs, typeall, fracsall, random_seed)
        print(f'[ SEQ: {seq_name} | 4 / 8 ] Statistical analysis complete.')

        # Get difference of scramble from input sequence
        difffromseq = list()
        for x in range(0, len(allscrvals)):
            difffromseq.append(sum(abs(myarr[0] - allscrvals[x]))) # if care about everything

        # Find most similar scramble
        val, idx = min((val, idx) for (idx, val) in enumerate(difffromseq))

        # Create nxn matrix of original sequence
        for x in range(0, myarr.shape[1]):
            if myarr[0, x] == 0:
                zvecdb[countseqs, x] = 0
            else:
                zvecdb[countseqs, x] = (myarr[0, x] - amean[x]) / math.sqrt(avar[x])

        # Create an nxn matrix of most similar scramble
        for x in range(0, len(allscrvals[idx])):
            if allscrvals[idx, x] == 0:
                zvecdbscr[countseqs, x] = 0
            else:
                zvecdbscr[countseqs, x] = (allscrvals[idx, x] - amean[x]) / math.sqrt(avar[x])

        num_types = len(typeall)
        reshaped_zvecdb = np.array(zvecdb[countseqs, :]).reshape((num_types, num_types))
        reshaped_zvecdbscr = np.array(zvecdbscr[countseqs, :]).reshape((num_types, num_types))

        zscore_savename = f'regular-{seq_name}.png'
        zscore_scrambled_savename = f'scrambled-{seq_name}.png'
        all_zscore_plots.append(zscore_savename)
        all_zscore_scrambled_plots.append(zscore_scrambled_savename)

        print(f'[ SEQ: {seq_name} | 5 / 8 ] Plotting Z-Score Matrix for FASTA: "{seq_name}"...')
        plot_zscore_matrix(seq_name, zvecdb, typeall, countseqs, zscore_savename, is_scrambled=False)
        print(f'[ SEQ: {seq_name} | 6 / 8 ] Plot of Z-Score matrix saved as: "{zscore_savename}".')

        print(f'[ SEQ: {seq_name} | 7 / 8 ] Plotting Scrambled Z-Score Matrix for FASTA: "{seq_name}"...')
        plot_zscore_matrix(seq_name, zvecdbscr, typeall, countseqs, zscore_scrambled_savename, is_scrambled=True)
        print(f'[ SEQ: {seq_name} | 8 / 8 ] Plot of Scrambled Z-Score matrix saved as: "{zscore_scrambled_savename}".')
        print(end='\n\n')

        to_export[seq_name] = (myseq, allscrseqs[idx], countseqs, zscore_savename, zscore_scrambled_savename, reshaped_zvecdb, reshaped_zvecdbscr)

    export_analysis(to_export)

    end = datetime.now()
    print('Total processing time: {}'.format(end - start))


def calculate_zscore(orthseqs, typeall, num_seqs, random_seed):
    olen        = len(orthseqs)
    tlen        = len(typeall)
    zvecdb      = np.zeros((olen, tlen**2))
    zvecdbscr   = np.zeros((olen, tlen**2))
    countseqs   = -1
    start       = datetime.now()

    to_export = dict()
    print(f'Number of amino acid groupings used for Nardini z-score analysis: {tlen}', end='\n\n')
    print(f'Amino acid groupings used for Nardini z-score analysis: ', end='\n\n')
    print_amino_acid_groupings_header(typeall)
    
    for seq_record in orthseqs:
        seq_name = str(seq_record.id)
        myseq = str(seq_record.seq)
        fracsall = list()
        countseqs = countseqs + 1

        print('=' * 80, end='\n\n')
        print(f'[ FASTA {seq_name} : {myseq} ]')
        print(f'[ SEQ: {seq_name} | START ] Beginning analysis of FASTA sequence: "{seq_name}"...', end='\n\n')

        for type1 in typeall:
            count = count_residues_in_sequence(myseq, type1)
            fraction = count / len(myseq)
            fracsall.append(fraction)

        print(f'[ SEQ: {seq_name} | 1 / 4 ] Analyzing original sequence...')
        myarr = get_org_seq_vals(myseq, typeall, fracsall)
        print(f'[ SEQ: {seq_name} | 2 / 4 ] Analysis of original sequence complete.')

        # Returns mean of scrambles, std of scrambles, all values in a number of scramble x tlen**2 list, and all scramble sequences
        print(f'[ SEQ: {seq_name} | 3 / 4 ] Performing analysis of {num_seqs} scrambled sequences: mean, stddev, etc...')
        alpha, amean, avar, allscrvals, allscrseqs = get_scramble_seqs_vals(myseq, num_seqs, typeall, fracsall, random_seed)
        print(f'[ SEQ: {seq_name} | 4 / 4 ] Statistical analysis complete.')

        # Get difference of scramble from input sequence
        difffromseq = list()
        for x in range(0, len(allscrvals)):
            difffromseq.append(sum(abs(myarr[0] - allscrvals[x]))) # if care about everything

        # Find most similar scramble
        val, idx = min((val, idx) for (idx, val) in enumerate(difffromseq))

        # Create an nxn matrix of original sequence
        for x in range(0, myarr.shape[1]):
            if myarr[0, x] == 0:
                zvecdb[countseqs, x] = 0
            else:
                zvecdb[countseqs, x] = (myarr[0, x] - amean[x]) / math.sqrt(avar[x])

        # Create an nxn matrix of most similar scramble
        for x in range(0, len(allscrvals[idx])):
            if allscrvals[idx, x] == 0:
                zvecdbscr[countseqs, x] = 0
            else:
                zvecdbscr[countseqs, x] = (allscrvals[idx, x] - amean[x]) / math.sqrt(avar[x])

        num_types = len(typeall)
        reshaped_zvecdb = np.array(zvecdb[countseqs, :]).reshape((num_types, num_types))
        reshaped_zvecdbscr = np.array(zvecdbscr[countseqs, :]).reshape((num_types, num_types))
        print(end='\n\n')

        to_export[seq_name] = (myseq, allscrseqs[idx], countseqs, reshaped_zvecdb, reshaped_zvecdbscr)

    end = datetime.now()
    print('Total processing time: {}'.format(end - start))
    return to_export
