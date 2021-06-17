# -*- coding: utf-8 -*-
import os
import math
import numpy as np
from datetime import datetime
from nardini.core import count_residues_in_sequence, typeall
from nardini.core import get_org_seq_vals, get_scramble_seqs_vals
from nardini.plotting import plot_zscore_matrix


def calculate_and_plot(orthseqs, typeall, num_seqs, random_seed):
    olen        = len(orthseqs)
    tlen        = len(typeall)
    zvec        = np.zeros((olen, int(tlen + (tlen * (tlen - 1))/2)))
    zvecdb      = np.zeros((olen, tlen**2))
    zvecdbscr   = np.zeros((olen, tlen**2))
    countseqs   = -1
    start       = datetime.now()
    for seq_record in orthseqs:
        seq_name = str(seq_record.id)
        myseq = str(seq_record.seq)
        fracsall = list()
        countseqs = countseqs + 1

        print(f'[ START ] Beginning analysis of FASTA sequence: "{seq_name}"...')

        for type1 in typeall:
            count = count_residues_in_sequence(myseq, type1)
            fraction = count / len(myseq)
            fracsall.append(fraction)

        print('[ 1 / 8 ] Analyzing original sequence...')
        myarr = get_org_seq_vals(myseq, typeall, fracsall)
        print('[ 2 / 8 ] Analysis of original sequence complete.')

        # Returns mean of scrambles, std of scrambles, all values in a number of scramble x 64 list, and all scramble sequences
        print(f'[ 3 / 8 ] Performing analysis of {num_seqs} scrambled sequences: mean, stddev, etc...')
        alpha, amean, avar, allscrvals, allscrseqs = get_scramble_seqs_vals(myseq, num_seqs, typeall, fracsall, random_seed)
        print('[ 4 / 8 ] Statistical analysis complete.')

        # Get difference of scramble from input sequence
        difffromseq = list()
        for x in range(0, len(allscrvals)):
            difffromseq.append(sum(abs(myarr[0] - allscrvals[x]))) # if care about everything

        # Find most similar scramble
        val, idx = min((val, idx) for (idx, val) in enumerate(difffromseq))

        # Create 8x8 matrix of original sequence
        for x in range(0, myarr.shape[1]):
            if myarr[0, x] == 0:
                zvecdb[countseqs, x] = 0
            else:
                zvecdb[countseqs, x] = (myarr[0, x] - amean[x]) / math.sqrt(avar[x])

        # Create 8x8 matrix of most similar scramble
        for x in range(0, len(allscrvals[idx])):
            if allscrvals[idx, x] == 0:
                zvecdbscr[countseqs, x] = 0
            else:
                zvecdbscr[countseqs, x] = (allscrvals[idx, x] - amean[x]) / math.sqrt(avar[x])

        zscore_savename = f'regular-{seq_name}.png'
        zscore_scrambled_savename = f'scrambled-{seq_name}.png'

        print(f'[ 5 / 8 ] Plotting Z-Score Matrix for FASTA: "{seq_name}"...')
        plot_zscore_matrix(seq_name, zvecdb, typeall, countseqs, zscore_savename, is_scrambled=False)
        print(f'[ 6 / 8 ] Plot of Z-Score matrix saved as: "{zscore_savename}".')

        print(f'[ 7 / 8 ] Plotting Scrambled Z-Score Matrix for FASTA: "{seq_name}"...')
        plot_zscore_matrix(seq_name, zvecdbscr, typeall, countseqs, zscore_scrambled_savename, is_scrambled=True)
        print(f'[ 8 / 8 ] Plot of Scrambled Z-Score matrix saved as: "{zscore_scrambled_savename}".')
        print(end='\n\n')
    end = datetime.now()
    print('Total processing time: {}'.format(end - start))
