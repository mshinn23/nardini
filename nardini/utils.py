# -*- coding: utf-8 -*-
#
# This module contains parallelization functions that are adapted from
# the `damfret_classifier.utils` module from https://github.com/Pappulab/damfret_classifier .
import os
from copy import copy
from io import StringIO
from Bio import SeqIO
import multiprocessing as mp
from collections import OrderedDict
from datetime import datetime
from nardini.constants import DEFAULT_RANDOM_SEED


def start_process():
    """This is just a notifier function to indicate that a worker process has been launched.
    Since this is meant for a pool, the pool reuses these processes."""
    print('Starting :: {}'.format(mp.current_process().name))


def parallelize(func, list_of_func_args, initializer=None, num_processes=None):
    """This function parallelizes calling a function with different argument values.

    @param func (func):                 The function reference which will be called.
    @param func_args (tuple | list):    The arguments to pass to the function call.
    @param initializer (callable):      The function to call when initializing the process
                                        pool.
    @param num_processes (int):         The number of processes to use (default = all).

    @returns results (OrderedDict):     A dictionary of function arguments as keys, and
    a values of 2-tuples comprised of datetime instances and the function result. Note:
    the result is a lazy reference to an `ApplyAsync` result. To extract the values,
    they can be obtained by calling `.get()` on the result.

    Notes: calling the `ApplyAsync` object's get (`result.get()`) here will block, which
    defeats the entire purpose behind this configuration.
    """
    # If no processes are set, create a pool from half the number of available
    # virtual CPUs (the number of virtual CPUs is chipset-specific; consult
    # your processor manufacturer for more).
    num_procs = num_processes
    if num_processes is None:
        num_procs = mp.cpu_count() // 2

    pool = mp.Pool(processes=num_procs, initializer=initializer)

    results = OrderedDict()
    for func_args in list_of_func_args:
        async_result = pool.apply_async(func, args=func_args)
        results[func_args] = (async_result, datetime.now())

    # Stop accepting jobs.
    pool.close()
    pool.join()

    # To obtain the results as they arrive, we need to use `AsyncResult.get()` - i.e.
    # `result.get()`. However, doing this in the processing loop will block and never
    # yield performance increases above those using 1 CPU. So, to parallelize the
    # results call, we do the `get()` after the pool has been closed and joined.
    #
    # The `actual_results` OrderedDict allows us to order the results as they arrive.
    # In Python 3.6+, `dict` supports ordering by default. For backwards compatability,
    # we use OrderedDict.
    #
    # We also return the start and end times along with the result as it provides
    # functionality to the user to sort results by any of those criteria.
    actual_results = OrderedDict()
    for fargs in results:
        async_result, start_time = results[fargs]
        actual_result = async_result.get()
        end_time = datetime.now()
        actual_results[fargs] = (actual_result, start_time, end_time)

    return actual_results


def read_sequences_from_filename(sequence_filename, default_name, verbose=False):
    """This is a helper function to read in sequences from a sequence FASTA file.
    This is a companion function to `read_sequences_from_string_list`.

    @param sequence_filename (str): The filepath of the file containing sequences.
                                    This file can contain FASTA records, or
                                    sequences that are saved on each newline.

    @param default_name (str):      This is the prefix name that's used when the
                                    file is found to only contain raw sequences
                                    on each new line.

    @param verbose (bool):          The extent to which information should be
                                    displayed throughout the analysis. Default:
                                    False.

    @returns seqio_sequences:       A list of sequence strings that were
                                    extracted from the sequence file.
    """
    seqio_sequences = list()
    if sequence_filename is None:
        raise RuntimeError('This parameter cannot be `None`. Exiting.')

    if type(sequence_filename) is not str:
        raise RuntimeError('This parameter can only be a string. Exiting.')

    if os.path.exists(sequence_filename):
        with open(sequence_filename, 'r') as seqfile:
            content = seqfile.read()
            if content.count('>') > 0:
                seqio_sequences = list(SeqIO.parse(sequence_filename, 'fasta'))
            else:
                # This means that there is a list of sequences that are separated
                # by newlines. We split by any whitespace so as to capture
                # carriage returns and line feeds.
                raw_sequences = [s.strip() for s in content.split() if len(s.strip()) > 0]
                seqio_sequences = read_sequences_from_string(raw_sequences, default_name)
    else:
        raise RuntimeError(f'Sequence filename: "{sequence_filename}" not found.')
    return seqio_sequences


def read_sequences_from_string_list(list_of_sequences, default_name, verbose=False):
    """This is a helper function to read in sequences from a list of raw sequences.
    This is a companion function to `read_sequences_from_filename`.

    @param sequence_filename (str): The filepath of the file containing sequences.
                                    This file can contain FASTA records, or
                                    sequences that are saved on each newline.

    @param default_name (str):      This is the prefix name that's used when the
                                    file is found to only contain raw sequences
                                    on each new line.

    @param verbose (bool):          The extent to which information should be
                                    displayed throughout the analysis. Default:
                                    False.

    @returns seqio_sequences:       A list of sequence strings that were
                                    extracted from the sequence file.
    """
    sequences = list()
    # This means that we have to create a fake record using the sequence content.
    for index, sequence in enumerate(list_of_sequences, start=1):
        fasta = f'>{default_name}-{index}\n{sequence}'
        fake_record = SeqIO.read(StringIO(fasta), 'fasta')
        sequences.append(fake_record)

    if verbose:
        print('Number of sequences read: {num}'.format(num=len(sequences)), end='\n\n')
    return sequences


def set_random_seed(seed):
    """This is a utility function to set the random seed. This function should be called
    ahead of any analysis.

    @param seed (int):          The random seed to use. Choosing a good random seed is
                                important for reproducibility and sampling quality.
    
    @return random_seed (int):  The random seed chosen or set.
    """
    random_seed = copy(seed)
    if seed == DEFAULT_RANDOM_SEED:
        random_seed = int(datetime.now().timestamp())
        print(f'No random seed specified. Generating random seed: {random_seed}\n')
    else:
        print(f'Using user-supplied random seed: {random_seed}\n')
    return random_seed
