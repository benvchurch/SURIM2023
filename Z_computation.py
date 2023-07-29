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
from utils.characters import UCF, gauss_sum, char

# Functions #

def gauss_sum_Z(p, k, y, d):
    q = p**k
    m = 0
    z = UCF.gen(q-1)
    # y = w**k for now use w**y
    for i in range(1, q-1):
        alpha = Integer(i)/Integer(q-1)
        m += conjugate(my_char(alpha, y, p, k)) * gauss_sum(alpha, p, k)**d

    m += (q-1)**(d-1) + 1
    return m / Integer(q*(q-1))


def char_prod(x, y, p, k):
    m = 0
    q = p**k
    for i in range(0, q-1):
        alpha = Integer(i)/Integer(q-1)
        m += my_char(alpha, x, p, k) * conjugate(my_char(alpha, y, p, k))

    return m




def calc_m(y, z, p, k, d):
    n = 0
    els = list(GF(p**k))
    all_els = product(list(GF(p**k)), repeat = d)
    for X in all_els:
        X_sum = 0
        X_prod = 1
        for x in X:
            X_sum += x
            X_prod *= x
        if X_sum == z and X_prod == y:
            n += 1
    return Integer(n)

def calc_Z(y, p, k, d):
    return calc_m(y, 0, p, k, d)/Integer(p**k - 1)



def my_char(alpha, y, p, k):
    if alpha == 0:
        return 1
    if y == 0:
        return 0
    if k == 1:
        g = primitive_root(p)
        r = Mod(y, p).log(g)
    else:
        g = GF(p**k).gen()
        r = y.log(g)
    z = UCF.gen(denominator(alpha))
    return z**(numerator(alpha)*r)

def main(p, k):
    q = p**k
    d = (q - 1)
    print("M-TEST")
    for y in GF(p):
        for z in GF(p):
            print(y,z, calc_m(y, z, p, k, d))
    print("DONE")
    print("CHAR TEST")
    for i in GF(p):
        for j in GF(p):
            print(i, j, char_prod(i, j, p, k))
    print("DONE")
    print("G-SUM TEST")
    for i in range(0, q-1):
        alpha = Integer(i)/Integer(q-1)
        m = 0
        for y in GF(p**k):
            if y != 0:
                m += calc_Z(y, p, k, d) * my_char(alpha, y, p, k)
        if i == 0:
            print(i, q*m - ((q-1)**(d-1) + 1), gauss_sum(alpha, p, k)**d)
        else:
            print(i, q*m, gauss_sum(alpha, p, k)**d)
    print("DONE")
    for i in GF(p**k):
        print(i, calc_Z(i, p, k, d),  gauss_sum_Z(p, k, i, d))

# Execution #

if __name__ == '__main__':
    p = 5
    k = 1
    print("M-TEST")
    for y in GF(p**k):
        print(y, calc_m(0, y, p, k, 3))
    print("DONE")
    main(p, k)
