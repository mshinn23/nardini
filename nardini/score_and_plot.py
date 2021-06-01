# -*- coding: utf-8 -*-
import os
import math
import numpy as np
from nardini.core import count_residues_in_sequence, typeall
from nardini.core import get_org_seq_vals, get_scramble_seqs_vals
from nardini.plotting import plot_zscore_matrix


def calculate_and_plot(orthseqs, typeall, num_seqs):
    olen        = len(orthseqs)
    tlen        = len(typeall)
    zvec        = np.zeros((olen, int(tlen + (tlen * (tlen - 1))/2)))
    zvecdb      = np.zeros((olen, tlen**2))
    zvecdbscr   = np.zeros((olen, tlen**2))
    countseqs   = -1
    for seq_record in orthseqs:
        seq_name = str(seq_record.id)
        myseq = str(seq_record.seq)
        fracsall = list()
        countseqs = countseqs + 1

        for type1 in typeall:
            count = count_residues_in_sequence(myseq, type1)
            fraction = count / len(myseq)
            fracsall.append(fraction)

        myarr = get_org_seq_vals(myseq, typeall, fracsall)

        # Returns mean of scrambles, std of scrambles, all values in a number of scramble x 64 list, and all scramble sequences
        alpha, amean, avar, allscrvals, allscrseqs = get_scramble_seqs_vals(myseq, num_seqs, typeall, fracsall)

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
        plot_zscore_matrix(seq_name, zvecdb, typeall, zscore_savename, is_scrambled=False)
        plot_zscore_matrix(seq_name, zvecdbscr, typeall, zscore_scrambled_savename, is_scrambled=True)