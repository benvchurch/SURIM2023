#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Calculations
Litt/Perry Group
Columbia University Mathematics REU 2018
Navtej Singh <singhnav@umich.edu>
'''

# Imports #

from sage.all import *

import numpy as np

import functools
import itertools

# Functions #

# Utilities

def primesfrom2to(n):
    '''
    https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    Input n>=6, Returns a array of primes, 2 <= p < n
    '''
    sieve = np.ones(n/3 + (n%6==2), dtype=np.bool)
    sieve[0] = False
    for i in range(int(n**0.5)/3+1):
        if sieve[i]:
            k=3*i+1|1
            sieve[((k*k)/3)::2*k] = False
            sieve[(k*k+4*k-2*k*(i&1))/3::2*k] = False
    return list(np.r_[2,3,((3*np.nonzero(sieve)[0]+1)|1)])

def numberToBaseN(number, base):
    base, number = int(base), int(number)
    if number == 0:
        return [0]
    else:
        digits = []
        while number:
            digits.append(int(number % base))
            number /= base
        return digits[::-1]

def getS(v, p):
    return sum(numberToBaseN(v, p))

def factors(n):
    if n in [0, 1]:
        return [n]
    return [f for f in list(set(functools.reduce(list.__add__, 
        ([i, n//i] for i in range(1, int(pow(n, 0.5) + 1)) if n % i == 0)))) if f != 1 and f != n]

def get_largest_factor(n):
    return max(factors(n))

# Main

def generate_gauss_sum_factors(min_prime=3, max_prime=100, max_power=10,
    outfile_path='data/analyses/gauss_sum_factor_calculations.csv'):
    with open(outfile_path, 'a') as outfile:
        outfile.write('q,p,f,mu,r\n')
        for p in [prime for prime in primesfrom2to(max_prime) if prime >= min_prime]:
            for f in range(1, max_power + 1):
                q = p ** f
                mus = [m for m in range(1, q - 1) if gcd(m, q - 1) == 1]
                rs = [r for r in sorted(itertools.product(*[list(range(1, q - 1))] * 4)) if len(set(r)) == len(r)]
                for mu in mus:
                    for r in rs:
                        r_ = [mu * element for element in r]
                        S = [getS(v, p) for v in r_]
                        if S[0] + S[1] == S[2] + S[3]:
                            line = '{},{},{},{},{}\n'.format(q, p, f, mu, '_'.join(map(str, r)))
                            outfile.write(line)
                            print(line)

def generate_counterexamples(prime, power, outfile_path='data/analyses/gauss_sum_factor_counterexamples'):
    with open('{}.{}.{}.csv'.format(outfile_path, prime, power), 'a') as outfile:
        outfile.write('q,p,f,n,largest_factor,mu_failure\n')
        q = prime ** power
        mus = [m for m in range(1, q - 1) if gcd(m, q - 1) == 1]
        rs = [r for r in sorted(itertools.product(*[list(range(1, q - 1))] * 4)) if len(set(r)) == len(r)]
        solutions = {}
        for r in rs:
            for mu in mus:
                count, largest_factor = 0, get_largest_factor(mu)
                r_ = [mu * element for element in r]
                S = [getS(v, prime) for v in r_]
                if S[0] + S[1] == S[2] + S[3]:
                    count += 1
                print('Q {}, Mu {}, Solutions {}\n'.format(q, mu, count))
                if largest_factor not in solutions:
                    solutions[largest_factor] = [(mu, count)]
                else:
                    solutions[largest_factor].append((mu, count))
        for largest_factor, v in solutions.items():
            if len(set(v)) > 1:
                for mu, count in v:
                    line = '{},{},{},{},{},{}\n'.format(q, prime, power, count, largest_factor, mu)
                    outfile.write(line)
                    print(line)

# Execution #

if __name__ == '__main__':
    generate_counterexamples(2, 3)
