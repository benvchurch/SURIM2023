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

from glob              import glob
from main              import run_affine, run_projective
from newton_polygons   import *
from plot              import plot_lines, plot_multiple_lines
from rational_function import *
from sage.all          import *
from smooth            import my_WP_smooth
from supersingular     import *
from time              import time
from variety           import AffineVariety, DefiningEquation, ProjectiveVariety


# Functions #

# Timing

import time

# Varities class usage:

if __name__ == '__main__':
    V = AffineVariety(19, DefiningEquation([1,1,1], [1000,1000,1000]))
    start_time = time.time()
    print(V.is_supersingular())
    print("--- %s seconds ---" % (time.time() - start_time))
