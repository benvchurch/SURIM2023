#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Explicit Calculation
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

from sage.all import *
from sage.matrix.operation_table import OperationTable
from itertools import product
from utils.characters import UCF, gauss_sum, char, add_char
import numpy as np
import time

# Functions #

def gen_gauss_sums(p, k):
    q = p**k
    els = []

    if k == 1:
        w = GF(p)(primitive_root(p))
    else:
        w = GF(q).gen()

    for i in range(0, q-1):
        els.append(add_char(w**i, p))

    X = np.array(els)

    return (q-1)*np.fft.ifft(X)

# tests

for k in range(1, 10):
    for p in [2,3,5,7,11]:
        q = p**k
        print("p = {}, k = {}, q = {}".format(p, k, q))
        t = time.time()
        sms = gen_gauss_sums(p,k)
        print("FFT in {} s".format(time.time() - t))

        t = time.time()
        for i in range(0, q-1):
            if abs( sms[i] - CC(gauss_sum(Integer(i)/Integer(q-1), p, k)) ) > 0.0001:
                print("OH NO: ", i)
        print("My Method in {} s".format(time.time() - t))
