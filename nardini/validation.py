import os
import sys
from nardini.constants import TYPEALL, TYPEALL_8x8, TYPEALL_9x9


def validate_arguments(args):
    """A helper function to check several of the command line options.

    @param args (argparse.Namespace):       The object from the argument parser
                                            whose options need to be checked.
    
    @returns amino_acid_groupings (list):   The groupings of amino acids to use."""
    sequences = args.sequences
    sequences_filename = args.sequences_filename
    num_amino_acid_types = args.num_amino_acid_types

    amino_acid_groupings = TYPEALL
    num_types_1 = len(TYPEALL_8x8)
    num_types_2 = len(TYPEALL_9x9)

    if sequences is not None and sequences_filename is not None:
        raise RuntimeError('Options for sequences or filename are both set. Please choose one.')

    if sequences_filename is str and not os.path.exists(sequences_filename):
        raise RuntimeError(f'No filename found despite being set ("{sequences_filename}").')

    if num_amino_acid_types <= 0:
        raise RuntimeError('Error: the number of amino acid types cannot be less than or equal to 0.')
    
    elif num_amino_acid_types < num_types_1:
        raise RuntimeError(f'Error: only groupings of {num_types_1} or {num_types_2} amino acids are currently allowed.')
    
    elif num_amino_acid_types == num_types_2:
        amino_acid_groupings = TYPEALL_9x9
    
    elif num_amino_acid_types > num_types_2:
        raise RuntimeError('Error: the number of amino acid types cannot be more than {num_types_2}.')

    if sequences is None and sequences_filename is None:
        raise RuntimeError('No sequences or filename passed for analysis. Please supply an option.')

    return amino_acid_groupings
