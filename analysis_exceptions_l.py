#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Analysis
Litt/Perry Group
Columbia University Mathematics REU 2018
'''

# Imports #

import itertools

from glob     import glob
from sage.all import *
from tests    import primesfrom2to

# Functions #

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

def one_exp_coprime_to_all_others(exponents):
    result = True
    for (i,n) in enumerate(exponents):
        if gcd(n, gcd([n_j for (j, n_j) in enumerate(exponents) if j != i])) != 1:
            result = False
    return result

def analyze_prime_mod_lcm_exp(infile_path='supersingular.csv',
    outfile_path='data/analyses/exceptions_2conditions_large.csv'):
    with open(outfile_path, 'a') as outfile:
        outfile.write('prime,exponents,supersingularity,prime_minus_one_mod_lcm_exponents_condition\n')
        for infile in glob(infile_path):
            for index, inline in enumerate(open(infile, 'r')):
                if index == 0:
                    continue
                try:
                    prime, exponents, numerator, denominator, supersingularity, error, speed = inline.split(',')
                    outline = '{},{},{},{}\n'.format(prime, exponents, supersingularity,
                        prime_mod_lcm_exp(int(prime), list(map(int, exponents.split('_')))))
                    if supersingularity.strip().lower() == 'true':
                        if prime_mod_lcm_exp(int(prime), list(map(int, exponents.split('_')))) == False:
                            if one_exp_coprime_to_all_others(list(map(int, exponents.split('_')))) == False:
                                outfile.write(outline)
                                print(outline)
                except ValueError:
                    print(index, inline, len(inline.split(',')))
                    print('Skipping Line {}: Improper number of columns'.format(index + 1))

if __name__ == '__main__':
    analyze_prime_mod_lcm_exp()
