#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Testing Correctness
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

from sage.all import *
import random
from utils.affine_patch_calculations import get_sorted_alpha_tuples
from utils.gauss_sum_factorization import s
from tests import primesfrom2to
from variety import *



def compute_newton_polygon(r, p, powers):
    m = lcm(powers)
    f = Mod(p, m).multiplicative_order()
    q = p**f
    D = dict()
    alphas = get_sorted_alpha_tuples(r, p, f, [1,1,1,1], powers)
    for alpha in alphas:
        S = 0
        for a in alpha[1]:
            j_num = a * m
            S += s((q - 1) * j_num / m, p, q)
        if S not in D:
            D[S] = alpha[0]
        else:
            D[S] += alpha[0]

    point = (0,0)
    points = []
    points.append(copy(point))
    for i in sorted(D):
        slope = Integer(i) / Integer((p - 1) * f) - 1
        if slope == 1:
            D[i] += 1
        point = (point[0]  + D[i], point[1] + D[i] * slope)
        points.append(copy(point))

    return points

print(compute_newton_polygon(3, 11, [5,5,5,5]))

for n in [19]:
    for p in [5,13]:
        V = ProjectiveVariety(p, DefiningEquation([1,1,1], [n,n,n]))
        print(V.newton_polygon())
        s = set()
        for i in range(euler_phi(n)):
            s.add(Mod(p**i, n))
        print(sorted(s))


'''
for alpha in alphas:
    is_root_of_unity = True
    for mu in Integers(m).list_of_elements_of_multiplicative_group():
        S = 0
        S_prime = 0
        for i in range(len(alpha[1])):
            j_num = alpha[1][i] * m
            S += s((q - 1) * mu * j_num / m, p, q)
            S_prime += (p-1) * Integer((mu * j_num) % m) / Integer(m)
        print mu, alpha[1], S, (p-1) * f * (r + 1)/2, S_prime
        if S != (p-1) * f * (r + 1)/2:
            print "RIP"
'''
