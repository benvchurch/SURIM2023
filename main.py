#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Main
Litt/Perry Group
Columbia University Mathematics REU 2018
Navtej Singh <singhnav@umich.edu>
'''

# Imports #

import cProfile
import gc
import itertools
import signal

import numpy as np
from variety           import AffineVariety, DefiningEquation, ProjectiveVariety


# Functions #

# Timing

import time

# Varities class usage:

if __name__ == '__main__':
    V = AffineVariety(5, DefiningEquation([1,1,1], [3,3,3]))
    start_time = time.time()
    print(V.is_supersingular())
    print("--- %s seconds ---" % (time.time() - start_time))
