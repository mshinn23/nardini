# -*- coding: utf-8 -*-
#
# This module contains parallelization functions that are adapted from
# the `damfret_classifier.utils` module from https://github.com/Pappulab/damfret_classifier .
import os
import multiprocessing as mp
from collections import OrderedDict
from datetime import datetime


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
