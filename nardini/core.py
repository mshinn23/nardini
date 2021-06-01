# -*- coding: utf-8 -*-
import numpy as np
import random
import scipy.stats as stats


NUM_SEQUENCES = 10000


pol = ['S','T','N','Q','C','H']
hyd = ['I','L','M','V']
pos = ['R','K']
neg = ['E','D']
aro = ['F','W','Y']
ala = ['A']
pro = ['P']
gly = ['G']


# We use a tuple since this should be immutable
typeall = (pol, hyd, pos, neg, aro, ala, pro, gly)


def count_residues_in_sequence(sequence, residue_types):
    residue_count = 0
    for residue in residue_types:
        residue_count += sequence.count(residue)
    return residue_count


def get_kappa(seq, type1, type2):
    blobsz = [5, 6]
    kappab = list()
    for b in blobsz:
        # Get full sequence asymmetry
        count1 = count_residues_in_sequence(seq, type1)
        count2 = count_residues_in_sequence(seq, type2)

        count1_frac = count1/len(seq)
        count2_frac = count2/len(seq)

        sigAll = (count1_frac - count2_frac)**2 / (count1_frac + count2_frac)

        # Get asymmetry for each blob
        sigX = list()
        for x in range(0, len(seq)-b+2):
            subseq = seq[x:x+b]

            ss_count1 = count_residues_in_sequence(subseq, type1)
            ss_count2 = count_residues_in_sequence(subseq, type2)

            ss_count1_frac = ss_count1/b
            ss_count2_frac = ss_count2/b

            if ss_count1 + ss_count2 == 0:
                sigX.append(0)
            else:
                sig = (ss_count1_frac - ss_count2_frac)**2 / (ss_count1_frac + ss_count2_frac)
                sigX.append(sig)

        asym = list()
        for x in range(0, len(sigX)):
            asym.append((sigX[x] - sigAll)**2)

        kappab.append(np.mean(asym))

    kappa = np.mean(kappab)
    return kappa


def get_omega(seq, type1):
    blobsz = [5, 6]
    omegab = []
    for b in blobsz:
        # Get full sequence asymmetry
        count = count_residues_in_sequence(seq, type1)
        count_frac = count/len(seq)

        sigAll = (count_frac - (1 - count_frac))**2

        # Get asymmetry for each blob
        sigX = list()
        for x in range(0, len(seq)-b+2):
            subseq = seq[x:x+b]
            ss_count = count_residues_in_sequence(subseq, type1)

            ss_count_frac = ss_count / b
            sig = (ss_count_frac - (1 - ss_count_frac))**2
            sigX.append(sig)

        asym = list()
        for x in range(0, len(sigX)):
            asym.append((sigX[x] - sigAll)**2)
        omegab.append(np.mean(asym))

    omega = np.mean(omegab)
    return omega


def get_org_seq_vals(myseq, typeall, fracsall):
    type_all_len = len(typeall)
    org_seq_arr = np.zeros((type_all_len, type_all_len))

    for count1 in range(0, type_all_len):
        type1 = typeall[count1]

        for count2 in range(count1, type_all_len):
            type2 = typeall[count2]

            if type1 == type2 and fracsall[count1] > 0.10:
                org_seq_arr[count1, count2] = get_omega(myseq, type1)

            if type1 != type2 and fracsall[count1] > 0.10 and fracsall[count2] > 0.10:
                org_seq_arr[count1, count2] = get_kappa(myseq,type1,type2)

    org_seq_1d = org_seq_arr.reshape([1, type_all_len**2])
    return org_seq_1d


def get_scramble_seqs_vals(myseq, num_seqs, typeall, fracsall, random_seed=None):
    if random_seed is not None:
        random.seed(random_seed)
        np.random.seed(random_seed)  # for `scipy.stats` (and consistency)

    currseq         = list()
    allseqs         = list()
    scr_vals        = np.zeros((num_seqs, len(typeall)**2))
    type_all_len    = len(typeall)

    for x in range(0, num_seqs):
        currseq = ''.join(random.sample(myseq, len(myseq)))
        scr_seq_arr = np.zeros((type_all_len, type_all_len))

        for count1 in range(0, type_all_len):
            type1 = typeall[count1]

            for count2 in range(count1, type_all_len):
                type2 = typeall[count2]

                if type1 == type2 and fracsall[count1] > 0.10:
                    scr_seq_arr[count1, count2] = get_omega(currseq, type1)

                if type1 != type2 and fracsall[count1] > 0.10 and fracsall[count2] > 0.10:
                    scr_seq_arr[count1, count2] = get_kappa(currseq, type1, type2)

        scr_vals[x, 0:type_all_len**2] = scr_seq_arr.reshape([1, type_all_len**2])
        allseqs.append(currseq)

    # Fit to a gamma distribution and obtain mean and variance
    alpha   = list()
    beta    = list()
    amean   = list()
    avar    = list()

    scr_vals_t = scr_vals.transpose()
    scr_vals_row = scr_vals_t.shape[0]
    for i in range(0, scr_vals_row):
        fit_alpha, fit_loc, fit_beta = stats.gamma.fit(scr_vals_t[i,:])

        cmean   = stats.gamma.mean(fit_alpha, fit_loc, fit_beta)
        cvar    = stats.gamma.var(fit_alpha, fit_loc, fit_beta)

        alpha.append(fit_alpha)
        beta.append(fit_beta)
        amean.append(cmean)
        avar.append(cvar)

    return (alpha, amean, avar, scr_vals, allseqs)
