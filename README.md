# Overview

This repository houses a Python 3 package that implements the Nardini method of uncovering non-random binary patterns within IDRs, as described in the paper: [Uncovering Non-random Binary Patterns Within Sequences of Intrinsically Disordered Proteins](https://doi.org/10.1016/j.jmb.2021.167373).

The Python package provides two interfaces: a command-line script `nardini`, and a programmatic interface for the development of more complex analysis or programs.

# Installation

This package can be installed in two ways: `setup.py` and `pip`. The following commands can also be run in a virtual Python environment (e.g. one created using [venv](https://docs.python.org/3/library/venv.html) and related tools, or with [Anaconda](https://www.anaconda.com/)). Although this package can be installed to the System's Python installation root (e.g. `/usr/bin/python3`) it is recommended to install the package in a virtual environment. If the installation is performed on the System Python, administrator permissions may be required (e.g. prepending `sudo` on a Linux/Mac system, or granting approval via Window's UAC).

## Setup.py

To install this package using `setup.py`, download a zip of this repository or clone it with git (i.e. `git clone https://github.com/mshinn23/nardini.git`). After uncompressing the ZIP archive or after cloning the repository, change to the sub-directory where `setup.py` is present. Then, run: `python setup.py install`. Any requisite packages and dependencies such as BioPython, will be download and installed.

## PIP

To install the package with `pip`, this can be performed via: `pip install git+https://github.com/mshinn23/nardini.git`.

# Command-Line Usage

Upon installation of the Python package, a command-line tool, `nardini` is made available. This can be verified by entering `which nardini` in a terminal. That command should report a filepath to a nardini script. Note that Nardini was installed to a virtual environment created by Conda or Python's Virtual Environment, that environment must be active.

The `nardini` script has a few options, which are enumerated below. These options can be output by entering `nardini -h`:

```
usage: nardini [-h] [-s RANDOM_SEED] [-n NUM_SCRAMBLED_SEQUENCES] sequences

positional arguments:
  sequences             The sequence, or filename containing FASTA sequences to be analyzed.

optional arguments:
  -h, --help            show this help message and exit
  -s RANDOM_SEED, --random-seed RANDOM_SEED
                        The random seed to use for generating the random distribution. (default:
                        None)
  -n NUM_SCRAMBLED_SEQUENCES, --num-scrambled-sequences NUM_SCRAMBLED_SEQUENCES
                        The number of sequences to generate while bootstrapping. (default: 100000)
  -t NUM_AMINO_ACID_TYPES, --num-amino-acid-types NUM_AMINO_ACID_TYPES
                        The number of amino acid type groupings. (default: 8)
```

If one had several IDR fastas saved to a file, e.g. `my_fastas.fsa`, then one can perform the analysis using the default options via: `nardini my_fastas.fsa`. If analysis on a sequence is needed (i.e. no FASTA header), such an analysis can be performed by passing the sequence directly to the script. For e.g.: `nardini PRQEFEVMEDHAGTYGLGDRKDQGGYTMHQ`.

Upon initiating the analysis a report will be generated that includes the number of sequences read as well as the random seed used. Progress on the analysis will be reported to the user. When the analysis is complete, a ZIP file containing both TSV representations of the Nardini analysis and its corresponding plots will be generated in the current directory.

It should be noted that the option to set a `RANDOM_SEED` via `-s` is implemented for generating reproducible runs. If no seed is supplied, a seed is selected and reported to the user, and the analysis proceeds.

The default number of scrambled sequences (1E5) is not ideal for all sequence lengths. For short sequences (<=8), a smaller number should be used as the analysis is likely to converge to the same outcome regardless of the set `RANDOM_SEED`. For longer sequences (>20), more sequence scrambles should be performed for adequate sampling.


# Package Usage Examples

## Simultaneous Analysis and Plots

To perform the Nardini z-score analysis and obtain its corresponding files and plots all at once, the following can be done. This example will produce a ZIP file.

```
import os
import time
from Bio import SeqIO
from io import StringIO
from nardini.core import NUM_SCRAMBLED_SEQUENCES, typeall, typeall_9x9
from nardini.score_and_plot import calculate_zscore_and_plot


# Define a random seed based on the current time.
# If a reproducible run is required, replace `int(time.time())`
# with an integer of one's choice.
RANDOM_SEED = int(time.time())
sequence_or_filename = 'your_sequence_or_path_here'
sequences = list()


# Prepare the sequences by using BioPython to read in FASTA records
# Or, create a fasta record from the sequence.
if os.path.exists(sequence_or_filename):
    with open(sequence_or_filename, 'r') as seqfile:
        sequences = list(SeqIO.parse(seqfile, 'fasta'))

elif type(sequence_or_filename) is str:
    # This means that we have to create a fake record using the sequence content.
    seq = sequence_or_filename[:]
    fasta = f'>fasta-1\n{seq}'
    fasta_record = SeqIO.read(StringIO(fasta), 'fasta')
    sequences.append(fasta_record)

# After preparing the sequences, we can finally perform the calculations and save
# the analysis. Here, `typeall` refers to the coarse-graining of amino acids by a given
# type to improve the search of similar sequences. The default number of types is
# 8. However, a 9x9 grouping is also available where Histidine is moved from the 
# Polar amino acids into its own grouping. If that grouping is desired, the user
# can use `typeall_9x9` instead.
#
# pol = ['S','T','N','Q','C','H']
# hyd = ['I','L','M','V']
# pos = ['R','K']
# neg = ['E','D']
# aro = ['F','W','Y']
# ala = ['A']
# pro = ['P']
# gly = ['G']
calculate_zscore_and_plot(sequences, typeall, NUM_SCRAMBLED_SEQUENCES, RANDOM_SEED)
```

This example will produce a filename with a randomly generated name such as: `nardini-data-SSDF5M91UQ.zip`. The ZIP file will contain the file `sequences.tsv` (the sequences analyzed). And, for every sequence included, 4 files additional will also be included: `regular-<seq_name>.png` (the plot of the Nardini matrix of the original sequence), `scrambled-<seq_name>.png` (the plot of the Nardini matrix of the closest matching scrambled sequence), `zscore-original-sequence-<seq_name>.tsv` (the text file corresponding to the z-score matrix of the original sequence), and `zscore-scrambled-sequence-<seq_name>.tsv` (the text file corresponding to the z-score matrix of the scrambled sequence).

## Isolated Analysis

This scenario is suited for a case where the user desires to perform additional analysis on the Nardini z-score matrices.

```
import os
import time
from Bio import SeqIO
from nardini.core import NUM_SCRAMBLED_SEQUENCES, typeall
from nardini.score_and_plot import calculate_zscore


# Define a random seed based on the current time.
# If a reproducible run is required, replace `int(time.time())`
# with an integer of one's choice.
RANDOM_SEED = int(time.time())
sequence_or_filename = 'your_sequence_or_path_here'
sequences = list()


if os.path.exists(sequence_or_filename):
    with open(sequence_or_filename, 'r') as seqfile:
        sequences = list(SeqIO.parse(seqfile, 'fasta'))

elif type(sequence_or_filename) is str:
    # This means that we have to create a fake record using the sequence content.
    seq = sequence_or_filename[:]
    fasta = f'>fasta-1\n{seq}'
    fasta_record = SeqIO.read(StringIO(fasta), 'fasta')
    sequences.append(fasta_record)

# Analyze the sequence(s) only.
#
# The `data_to_export` is a dictionary where the keys are the names of the FASTA records
# and the keys are tuples of length 5. The tuples are ordered by variables which contain
#
# 1. Original_sequence.
# 2. Scrambled_sequence.
# 3. The sequence number (used for book-keeping for many sequences).
# 4. The `reshaped_zvecdb` corresponding to the original sequence.
# 5. The `reshaped_zvecdbscr` corresponding to the scrambled sequences.
#
# Here, `typeall` refers to the coarse-graining of amino acids by a given
# type to improve the search of similar sequences. The default number of types is
# 8:
#
# pol = ['S','T','N','Q','C','H']
# hyd = ['I','L','M','V']
# pos = ['R','K']
# neg = ['E','D']
# aro = ['F','W','Y']
# ala = ['A']
# pro = ['P']
# gly = ['G']
data_to_export = calculate_zscore(sequences, typeall, NUM_SCRAMBLED_SEQUENCES, RANDOM_SEED)
```

# Citation

If you use the Nardini package for scientific research that will be published, please cite it. For e.g., the BibTeX entry for the Nardini manuscript is:

```
@article{cohan2021uncovering,
  title={Uncovering non-random binary patterns within sequences of intrinsically disordered proteins},
  author={Cohan, Megan C and Shinn, Min Kyung and Lalmansingh, Jared M and Pappu, Rohit V},
  journal={Journal of Molecular Biology},
  pages={167373},
  year={2021},
  publisher={Elsevier}
}
```

# Contributing

If you're interested in contributing to the development of Nardini - be it bug fixes, suggestions, edits, or extensions - please do the following. To submit identified bugs as well as suggestions (including ideas), please open a ticket under [Issues](https://github.com/mshinn23/nardini/issues); and, for direct edits and extensions to the package, please submit a [Pull Request](https://github.com/mshinn23/nardini/pulls).

This is an overview of [Github's Issues](https://docs.github.com/en/issues/tracking-your-work-with-issues/about-issues), as well as [Github's Pull Requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests) and how they can be used to aid development.
