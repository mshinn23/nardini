# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


def plot_zscore_matrix(seq_name, zvec_db, typeall, index, savename, is_scrambled):
    x_label_list = ['µ', 'h', '+', '-', 'π', 'A', 'P', 'G']
    num_types = len(typeall)
    reshaped_zvec_db = np.array(zvec_db[index, :]).reshape((num_types, num_types))

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    img = ax.imshow(reshaped_zvec_db, vmin=-3, vmax=3, cmap='bwr', aspect='auto')
    fig.colorbar(img)

    if is_scrambled:
        ax.set_title('Z-Score Matrix With the Most Similar Scramble to Sequence {}'.format(seq_name))
    else:
        ax.set_title('Z-Score Matrix for Sequence {}'.format(seq_name))

    ax.set_xticks(list(range(len(x_label_list))))
    ax.set_xticklabels(x_label_list)
    ax.set_xlabel('Amino Acid Type')

    ax.set_yticks(list(range(len(x_label_list))))
    ax.set_yticklabels(x_label_list)
    ax.set_ylabel('Amino Acid Type')

    fig.savefig(savename)
