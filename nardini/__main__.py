#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This script is an entry point that combines several functions from `nardini` into a
# stand-alone script that provides the ability to provide stand-alone analysis.
import os
import sys
import signal
from io import StringIO
from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from nardini.constants import NUM_SCRAMBLED_SEQUENCES, DEFAULT_RANDOM_SEED, DEFAULT_PREFIX_NAME
from nardini.utils import read_sequences_from_filename, read_sequences_from_string_list
from nardini.utils import set_random_seed
from nardini.constants import TYPEALL, TYPEALL_8x8, TYPEALL_9x9
from nardini.validation import validate_arguments
from nardini.score_and_plot import calculate_zscore_and_plot


# -------------------------------------------------------------------------------------------------
# Script-specific utilities


def shutdown(signal_number, frame):
    """A utility function to handle `Ctrl-C` via the command-line."""
    print('Shut down requested. Terminating process.')
    sys.exit(0)


def shutdown_not_suspend(signal_number, frame):
    """A helper function to handle `Ctrl-Z` via the command-line."""
    print('Suspend not supported. Shutting down instead. Terminating process.')
    sys.exit(0)


# -------------------------------------------------------------------------------------------------


def main():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('--sequences',
                        '-s',
                        help='The sequences to analyze.',
                        type=str,
                        nargs='+',
                        default=None)

    parser.add_argument('--sequences-filename',
                        '-f',
                        help='The filename containing the sequences for analysis.',
                        type=str,
                        default=None)

    parser.add_argument('-r', '--random-seed',
                        help='The random seed to use for generating the random distribution.',
                        type=int,
                        default=DEFAULT_RANDOM_SEED)

    parser.add_argument('-n', '--num-scrambled-sequences',
                        help='The number of sequences to generate while bootstrapping.',
                        type=int,
                        default=NUM_SCRAMBLED_SEQUENCES)
    
    parser.add_argument('-t', '--num-amino-acid-types',
                        help='The number of amino acid type groupings.',
                        type=int,
                        choices=[len(TYPEALL_8x8), len(TYPEALL_9x9)],
                        default=len(TYPEALL))

    parser.add_argument('--fasta-name-prefix',
                        '-p',
                        help='The prefix name to use for generated FASTA records',
                        type=str,
                        default=DEFAULT_PREFIX_NAME)

    parser.add_argument('--verbose',
                        '-v',
                        help='Output verbose status messages during analysis.',
                        action='store_true')

    args = parser.parse_args()

    sequences = args.sequences
    sequences_filename = args.sequences_filename
    verbose = args.verbose
    num_scrambled_sequences = args.num_scrambled_sequences
    fasta_name_prefix = args.fasta_name_prefix
    
    amino_acid_groupings = validate_arguments(args)
    random_seed = set_random_seed(args.random_seed)
    
    if sequences is not None:
        parsed_sequences = read_sequences_from_string_list(sequences, fasta_name_prefix, verbose)

    if sequences_filename is not None:
        parsed_sequences = read_sequences_from_filename(sequences_filename, fasta_name_prefix, verbose)

    calculate_zscore_and_plot(parsed_sequences, amino_acid_groupings, num_scrambled_sequences, random_seed)


if __name__ == '__main__':
    # Register the signals that will be handled by this script.
    signal.signal(signal.SIGTSTP,   shutdown_not_suspend)   # Ctrl-Z
    signal.signal(signal.SIGINT,    shutdown)               # Ctrl-C

    # Run as normal otherwise.
    main()
