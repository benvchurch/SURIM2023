#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Analysis
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

import gc
import itertools

from glob              import glob
from rational_function import *
from sage.all          import gcd, lcm
from tests             import primesfrom2to
from variety           import *

# Functions #

# Utilities

'''
Code below used to express a prime as the sum of two squares
https://math.stackexchange.com/questions/5877/efficiently-finding-two-squares-which-sum-to-a-prime
'''

def mods(a, n):
    if n <= 0:
        raise Exception("negative modulus")
    a = a % n
    if (2 * a > n):
        a -= n
    return a

def powmods(a, r, n):
    out = 1
    while r > 0:
        if (r % 2) == 1:
            r -= 1
            out = mods(out * a, n)
        r /= 2
        a = mods(a * a, n)
    return out

def quos(a, n):
    if n <= 0:
        raise Exception("negative modulus")
    return (a - mods(a, n))/n

def grem(w, z):
    # remainder in Gaussian integers when dividing w by z
    (w0, w1) = w
    (z0, z1) = z
    n = z0 * z0 + z1 * z1
    if n == 0:
        raise Exception("division by zero")
    u0 = quos(w0 * z0 + w1 * z1, n)
    u1 = quos(w1 * z0 - w0 * z1, n)
    return(w0 - z0 * u0 + z1 * u1,
           w1 - z0 * u1 - z1 * u0)

def ggcd(w, z):
    while z != (0,0):
        w, z = z, grem(w, z)
    return w

def root4(p):
    # 4th root of 1 modulo p
    if p <= 1:
        raise Exception("too small")
    if (p % 4) != 1:
        raise Exception("not congruent to 1")
    k = p/4
    j = 2
    while True:
        a = powmods(j, k, p)
        b = mods(a * a, p)
        if b == -1:
            return a
        if b != 1:
            raise Exception("not prime")
        j += 1

def sq2(p):
    a = root4(p)
    return ggcd((p,0),(a,1))

# Analysis

def prime_mod_lcm_exp(prime, exponents):
    '''
    N.B. For Affine Case Only
    For each 0 <= i <= r, let L_i = lcm({n_j}|j =/= i) and let n_i_prime = gcd(n_i, L_i). 
    Then take k = lcm of all n_i_prime and determine if there exists a prime power q such that q ~ -1 mod k.
    Refer Theorem 1.1 in Conjectures/Theorems document.
    '''
    if any([n % prime == 0 for n in exponents]):
        return None
    n_primes = [gcd(n, lcm([n_j for (j, n_j) in enumerate(exponents) if j != i])) for (i, n) in enumerate(exponents)]
    lcm_n_primes = lcm(n_primes)
    power = 1
    while prime ** power % lcm_n_primes != 1:
        if power > lcm_n_primes + 5:
            return None
        if prime ** power % lcm_n_primes - lcm_n_primes == -1:
            return True
        power += 1
    return False

def analyze_prime_mod_lcm_exp(infile_path='data/affine/surface_search_over_primes.*.total.csv',
    outfile_path='data/analyses/prime_mod_lcm_exp_analysis.csv'):
    with open(outfile_path, 'a') as outfile:
        outfile.write('prime,exponents,supersingularity,prime_minus_one_mod_lcm_exponents_condition\n')
        for infile in glob(infile_path):
            for index, inline in enumerate(open(infile, 'r')):
                try:
                    prime, exponents, numerator, denominator, supersingularity, error, speed = inline.split(',')
                    outline = '{},{},{},{}\n'.format(prime, exponents, supersingularity, 
                        prime_mod_lcm_exp(int(prime), list(map(int, exponents.split('_')))))
                    outfile.write(outline)
                    print(outline)
                except ValueError:
                    print('Skipping Line {}: Improper number of columns'.format(index + 1))

def conjecture_surface_exp_4(min_prime=2, max_prime=250):
    total = []
    primes = [int(p) for p in primesfrom2to(max_prime) if p >= min_prime and p % 4 == 1]
    defining_equation = DefiningEquation([1, 1, 1, 1], [4, 4, 4, 4])
    for prime in primes:
        try:
            v = ProjectiveVariety(prime, defining_equation)
            x, y = sq2(prime)
            actual = int(zeta_function(v)[1].factor()[-1][0].coefficients()[1]*(prime**2))
            predicted = 2*(x**2-y**2)
            result = actual == predicted or -actual == predicted
            total.append(result)
            print('Prime {}, {}, {}, {}'.format(prime, 'Pass' if result else 'Fail', actual, predicted))
        except ArithmeticError:
            print('Skipping Prime {}'.format(prime))
    print('All Pass {}'.format(all(total)))

if __name__ == '__main__':
    conjecture_surface_exp_4()
