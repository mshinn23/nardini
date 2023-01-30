# Overview

This repository contains a pure Python 3 package that implements the Nardini method of uncovering non-random binary patterns within IDRs, as described in the paper: [Uncovering Non-random Binary Patterns Within Sequences of Intrinsically Disordered Proteins](https://doi.org/10.1016/j.jmb.2021.167373).

Two interfaces are provided: 1) a command-line script `nardini`; and, 2) a programmatic interface for the development of analysis pipelines and programs.

# Installation

This package can be installed in one of two ways: `pip`, or `conda`. UNIX and Unix-like systems (Linux) are supported natively, and can be installed as described in their respective subsections below.

Windows is a special case, as the Windows Subsystem for Linux (WSL) is the recommended means of managing `nardini` installations. To install and configure WSL, one can follow this excellent guide from [How-To Geek](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/).

## PIP

To install the package with `pip`, this can be performed via: `pip install nardini`. Although this package can be installed to the System's Python installation root (e.g. `/usr/bin/python3`) it is recommended to install the package in a virtual environment. One created using [venv](https://docs.python.org/3/library/venv.html) and related tools, or with [Anaconda](https://www.anaconda.com/). 

If the installation is performed on the System Python, administrator permissions may be required (e.g. prepending `sudo` on a Linux/Mac system, or granting approval via Window's UAC).

## PIP - Local Build and Install

To install this package from source, download a zip of the source repository or clone it with `git` (i.e. `git clone https://github.com/mshinn23/nardini.git`). After decompressing the ZIP archive or after cloning the repository, change to the sub-directory where `pyproject.toml` is present. Then, run: `pip install .` to download and install any dependencies, build the package as a wheel, and install it. 

This approach is preferred for development use as it allows the user the capability to make edits to the package as needed. In such a scenario, the package should be installed via `pip install -e .` (i.e. an editable install), which allows changes to the source code to be made immediately available for testing purposes.


## CONDA

Similar to `pip`, this package can be installed via `conda install nardini -c conda-forge`. Note that since Anaconda also supports [interoperability with pip](https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/pip-interoperability.html), `nardini` could be installed into an Anaconda environment using `pip`. However, as this feature is experimental, installation via this method may be bug prone.


# Command-Line Usage

Upon installation of the Python package, a command-line tool, `nardini` is made available. This can be verified by entering `command -v nardini` or `which nardini` in a terminal. Either command should report a filepath to a nardini script. Note that Nardini was installed to a virtual environment created by Conda or Python's Virtual Environment, that environment must be active.

The `nardini` script has a few options, which are enumerated below. These options can be output by entering `nardini -h`:

```
$ nardini -h
usage: nardini [-h] [--sequences SEQUENCES [SEQUENCES ...]]
               [--sequences-filename SEQUENCES_FILENAME] [-r RANDOM_SEED]
               [-n NUM_SCRAMBLED_SEQUENCES] [-t {8,9}]
               [--fasta-name-prefix FASTA_NAME_PREFIX] [--verbose]

optional arguments:
  -h, --help            show this help message and exit
  --sequences SEQUENCES [SEQUENCES ...], -s SEQUENCES [SEQUENCES ...]
                        The sequences to analyze. (default: None)
  --sequences-filename SEQUENCES_FILENAME, -f SEQUENCES_FILENAME
                        The filename containing the sequences for analysis. (default:
                        None)
  -r RANDOM_SEED, --random-seed RANDOM_SEED
                        The random seed to use for generating the random
                        distribution. (default: 0)
  -n NUM_SCRAMBLED_SEQUENCES, --num-scrambled-sequences NUM_SCRAMBLED_SEQUENCES
                        The number of sequences to generate while bootstrapping.
                        (default: 100000)
  -t {8,9}, --num-amino-acid-types {8,9}
                        The number of amino acid type groupings. (default: 8)
  --fasta-name-prefix FASTA_NAME_PREFIX, -p FASTA_NAME_PREFIX
                        The prefix name to use for generated FASTA records (default:
                        fasta)
  --verbose, -v         Output verbose status messages during analysis. (default:
                        False)
```

## Command-Line Examples

If one had several IDR fastas saved to a file, e.g. `my_fastas.fsa`, then one can perform the analysis using the default options via: `nardini -f my_fastas.fsa`. If analysis on several sequences is needed (i.e. no FASTA header), such an analysis can be performed by passing the sequences directly to the script. For e.g.: `nardini -s ASEQUENCE ANOTHERSEQUENCE YETANOTHERSEQUENCE`.

Upon initiating the analysis a report will be generated that includes the number of sequences read as well as the random seed used. Progress on the analysis will be reported to the user. When the analysis is complete, a ZIP file containing both TSV representations of the Nardini analysis and its corresponding plots will be generated in the current directory.

It should be noted that the option to set a `RANDOM_SEED` via `-r` is implemented for generating reproducible runs. If no seed is supplied, a seed is selected using the current timestamp and reported to the user. Analysis then proceeds.

The default number of scrambled sequences (100,000) is not ideal for all sequence lengths. For short sequences (<=8), a smaller number should be used as the analysis is likely to converge to the same outcome regardless of the chosen `RANDOM_SEED`. For longer sequences (>20), more sequence scrambles should be performed for adequate sampling.


# Package Usage Examples

These examples illustrate how to use the package to analyze sequences stored as strings or [FASTA](https://en.wikipedia.org/wiki/FASTA_format).

In each case, a filename with a randomly generated name such as `nardini-data-SSDF5M91UQ.zip` will be produced along with PNGs of the Nardini Z-score matrix. The ZIP file will contain the file `sequences.tsv` (the sequences analyzed). And, for every sequence included, 4 files additional will also be included: `regular-<seq_name>.png` (the plot of the Nardini matrix of the original sequence), `scrambled-<seq_name>.png` (the plot of the Nardini matrix of the closest matching scrambled sequence), `zscore-original-sequence-<seq_name>.tsv` (the text file corresponding to the z-score matrix of the original sequence), and `zscore-scrambled-sequence-<seq_name>.tsv` (the text file corresponding to the z-score matrix of the scrambled sequence).

## Simultaneous Analysis and Plots - Sequence strings

To perform the Nardini z-score analysis and obtain its corresponding files and plots all at once given input sequences only, the following can be done.

```
import os
import time
from Bio import SeqIO
from io import StringIO
from nardini.constants import NUM_SCRAMBLED_SEQUENCES, TYPEALL, TYPEALL_9x9
from nardini.score_and_plot import calculate_zscore_and_plot
from nardini.utils import read_sequences_from_string_list

# Define a random seed based on the current time.
# If a reproducible run is required, replace `int(time.time())`
# with an integer of one's choice.
RANDOM_SEED = int(time.time())
prefix_name = 'fasta'
sequences = ['YOURSEQUENCE', 'ANOTHERSEQUENCE']
fasta_sequences = read_sequences_from_string_list(sequences, prefix_name)

# After preparing the sequences, we can perform the calculations and save
# the analysis. Here, `TYPEALL` refers to the coarse-graining of amino acids by a given
# type to improve the search of similar sequences. The default number of types is
# 8. However, a 9x9 grouping is also available where Histidine is moved from the 
# Polar amino acids into its own grouping. If that grouping is desired, the user
# can use `TYPEALL_9x9` instead.
#
# pol = ['S','T','N','Q','C','H']
# hyd = ['I','L','M','V']
# pos = ['R','K']
# neg = ['E','D']
# aro = ['F','W','Y']
# ala = ['A']
# pro = ['P']
# gly = ['G']
calculate_zscore_and_plot(fasta_sequences, TYPEALL, NUM_SCRAMBLED_SEQUENCES, RANDOM_SEED)
```

## Simultaneous Analysis and Plots - FASTA

In the case where the sequences are in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format), the following illustrates how they can be analyzed:

```
import os
import time
from Bio import SeqIO
from io import StringIO
from nardini.constants import NUM_SCRAMBLED_SEQUENCES, TYPEALL, TYPEALL_9x9
from nardini.score_and_plot import calculate_zscore_and_plot
from nardini.utils import read_sequences_from_filename

# Define a random seed based on the current time.
# If a reproducible run is required, replace `int(time.time())`
# with an integer of one's choice.
RANDOM_SEED = int(time.time())
prefix_name = 'fasta'
filename = 'mysequences.fasta'
fasta_sequences = read_sequences_from_filename(filename, prefix_name)

# After preparing the sequences, we can perform the calculations and save
# the analysis. Here, `TYPEALL` refers to the coarse-graining of amino acids by a given
# type to improve the search of similar sequences. The default number of types is
# 8. However, a 9x9 grouping is also available where Histidine is moved from the 
# Polar amino acids into its own grouping. If that grouping is desired, the user
# can use `TYPEALL_9x9` instead.
#
# pol = ['S','T','N','Q','C','H']
# hyd = ['I','L','M','V']
# pos = ['R','K']
# neg = ['E','D']
# aro = ['F','W','Y']
# ala = ['A']
# pro = ['P']
# gly = ['G']
calculate_zscore_and_plot(fasta_sequences, TYPEALL, NUM_SCRAMBLED_SEQUENCES, RANDOM_SEED)
```

## Isolated Analysis

This scenario is suited for a case where the user desires to perform additional analysis on the Nardini z-score matrices, or where only the analysis is needed - plots are excluded.

```
import os
import time
from pprint import pprint
from Bio import SeqIO
from io import StringIO
from nardini.constants import NUM_SCRAMBLED_SEQUENCES, TYPEALL, TYPEALL_9x9
from nardini.utils import read_sequences_from_filename
from nardini.score_and_plot import calculate_zscore


# Define a random seed based on the current time.
# If a reproducible run is required, replace `int(time.time())`
# with an integer of one's choice.
RANDOM_SEED = int(time.time())
prefix_name = 'fasta'
filename = 'mysequences.fasta'
fasta_sequences = read_sequences_from_filename(filename, prefix_name)

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
# Here, `TYPEALL` refers to the coarse-graining of amino acids by a given
# type to improve the search of similar sequences. The default number of types is
# 8. However, a 9x9 grouping is also available where Histidine is moved from the 
# Polar amino acids into its own grouping. If that grouping is desired, the user
# can use `TYPEALL_9x9` instead.
#
# pol = ['S','T','N','Q','C','H']
# hyd = ['I','L','M','V']
# pos = ['R','K']
# neg = ['E','D']
# aro = ['F','W','Y']
# ala = ['A']
# pro = ['P']
# gly = ['G']
data_to_export = calculate_zscore(fasta_sequences, TYPEALL, NUM_SCRAMBLED_SEQUENCES, RANDOM_SEED)
pprint(data_to_export)
```

# Citing

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
