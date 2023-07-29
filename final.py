#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Final Calculations
Litt/Perry Group
Columbia University Mathematics REU 2018
Navtej Singh <singhnav@umich.edu>
'''

# Imports #

import itertools
import numpy as np
import os

from sage.all import *
from time     import time
from variety  import AffineVariety, DefiningEquation

# Functions #

# Utilities

def output_path_exists(output_filepath):
    output_directory = '/'.join(output_filepath.split('/')[:-1])
    if not os.path.isdir(output_directory):
        raise Exception('(ERROR) Final - Specified output filepath does not exist')
    return True

def congruence(number, congruent, modulus):
    return Mod(number, modulus) == congruent

def coprime(a, b):
    return gcd(a, b) == 1

def generate_prime_tuples(minimum, maximum, tuple_length):
    primes = generate_primes(minimum, maximum)
    return itertools.combinations(primes, tuple_length)

def generate_primes(minimum, maximum):
    '''
    https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    Input n>=6, Returns a array of primes, 2 <= p < n
    '''
    maximum += 1
    sieve = np.ones(int(maximum/3) + (maximum%6==2), dtype=np.bool)
    sieve[0] = False
    for i in range(int(maximum**0.5) / 3 + 1):
        if sieve[i]:
            k=3*i+1|1
            sieve[((k*k)/3)::2*k] = False
            sieve[(k*k+4*k-2*k*(i&1))/3::2*k] = False
    return [Integer(prime) for prime in list(np.r_[2,3,((3*np.nonzero(sieve)[0]+1)|1)]) if minimum <= prime <= maximum]

def generate_variety(n0, n1, n2, n3, omega):
    return AffineVariety(omega, DefiningEquation([1, 1, 1, 1], [n0, n1, n2, n3]))

def multiplicative_order(number, modulus):
    try:
        return Mod(number, modulus).multiplicative_order()
    except ArithmeticError:
        return None

def primitive_root(number, modulus):
    try:
        return multiplicative_order(number, modulus) == modulus - 1
    except ArithmeticError:
        return None

def Shioda(number, modulus):
    power = 1
    while not congruence(number, 1, modulus):
        if power > modulus + 5:
            return False
        if Mod(number ** power, modulus) - modulus == -1:
            return True
        power += 1
    return False

# Data Collection

def check_conjecture_1(p, q, r, omega):
    '''
    1. Instantiate variety
    2. Return its supersingularity (primarily interested in those cases that are)
    3. Return the values of the following five conditions
        (i) omega generates (Z/pZ)^x, i.e. omega is a primitive root mod p
        (ii) omega generates (Z/qZ)^x, i.e. omega is a primitive root mod q
        (iii) omega is congruent to 1 mod s
        (iv) p is congruent to 1 mod s
        (v) q is congruent to 1 mod s
    '''
    variety = generate_variety(p, r*p, q, r*q, omega)
    supersingularity = variety.is_supersingular()
    condition_1 = primitive_root(omega, p)
    condition_2 = primitive_root(omega, q)
    condition_3 = congruence(omega, 1, r)
    condition_4 = congruence(p, 1, r)
    condition_5 = congruence(q, 1, r)
    return supersingularity, condition_1, condition_2, condition_3, condition_4, condition_5

def check_conjecture_2(p, q, r, omega):
    '''
    1. Instantiate variety
    2. Return its supersingularity
    3. Return the values of the following 3 conditions
        (i) ord_p(omega)
        (ii) ord_q(omega)
        (iii) ord_r(omega)
    '''
    variety = generate_variety(p, q, r, p*q*r, omega)
    supersingularity = variety.is_supersingular()
    condition_1 = multiplicative_order(omega, p)
    condition_2 = multiplicative_order(omega, q)
    condition_3 = multiplicative_order(omega, r)
    return supersingularity, condition_1, condition_2, condition_3

def check_exception_1(p, q, r, omega):
    '''
    Find supersingular varieties given p, q, r (not even necessarily prime) which are exceptions to the Shioda condition, 
    but don't satisfy all conditions from Conjecture 1
    '''
    variety = generate_variety(p, r*p, q, r*q, omega)
    supersingularity = variety.is_supersingular()
    condition_1 = primitive_root(omega, p) or False
    condition_2 = primitive_root(omega, q) or False
    condition_3 = congruence(omega, 1, r)
    condition_4 = congruence(p, 1, r)
    condition_5 = congruence(q, 1, r)
    exception = not all([condition_1, condition_2, condition_3, condition_4, condition_5])
    if supersingularity and exception:
        return p, q, r, omega
    return None

def iterate_conjecture_inputs(minimum_prime, maximum_prime, minimum_omega, maximum_omega, exponents_prime=True):
    triples = generate_prime_tuples(minimum_prime, maximum_prime, 3) if exponents_prime is True else [list(map(Integer, triple)) for triple in itertools.combinations(range(minimum_prime, maximum_prime), 3)]
    for p, q, r in triples:
        omegas = [prime for prime in generate_primes(minimum_omega, maximum_omega) if not Shioda(prime, p*q*r)]
        for omega in omegas:
            yield p, q, r, omega

def range_conjecture_1(minimum_prime, maximum_prime, minimum_omega, maximum_omega, output_filepath=None, r_coprime=True):
    '''
    1. Find all distinct combinations of primes p, q, r such that (pq, r) = 1
    2. Find all omegas that are non-Shioda mod lcm(p, q, s)
    3. For each such combination of p, q, and s, and for each such omega, check_conjecture_1
    4. Record data
    '''
    if output_filepath and output_path_exists(output_filepath):
        with open(output_filepath, 'a') as outfile:
            header = 'p,q,r,omega,supersingularity,condition_1,condition_2,condition_3,condition_4,condition_5\n'
            outfile.write(header)
            print(header.strip())
            for p, q, r, omega in iterate_conjecture_inputs(minimum_prime, maximum_prime, minimum_omega, maximum_omega):
                if r_coprime:
                    for r in [Integer(n) for n in range(1, maximum_prime)]:
                        supersingularity, condition_1, condition_2, condition_3, condition_4, condition_5 = check_conjecture_1(p, q, r, omega)
                        line = '{},{},{},{},{},{},{},{},{},{}\n'.format(p, q, r, omega, supersingularity, condition_1, condition_2, condition_3, condition_4, condition_5)
                        outfile.write(line)
                        print(line.strip())
                else:
                    supersingularity, condition_1, condition_2, condition_3, condition_4, condition_5 = check_conjecture_1(p, q, r, omega)
                    line = '{},{},{},{},{},{},{},{},{},{}\n'.format(p, q, r, omega, supersingularity, condition_1, condition_2, condition_3, condition_4, condition_5)
                    outfile.write(line)
                    print(line.strip())      
    else:
        print('p,q,r,omega,supersingularity,condition_1,condition_2,condition_3,condition_4,condition_5')
        for p, q, r, omega in iterate_conjecture_inputs(minimum_prime, maximum_prime, minimum_omega, maximum_omega):
            supersingularity, condition_1, condition_2, condition_3, condition_4, condition_5 = check_conjecture_1(p, q, r, omega)
            line = '{},{},{},{},{},{},{},{},{},{}'.format(p, q, r, omega, supersingularity, condition_1, condition_2, condition_3, condition_4, condition_5)
            print(line)

def range_conjecture_2(minimum_prime, maximum_prime, minimum_omega, maximum_omega, output_filepath=None):
    '''
    1. Find all distinct combinations of primes p, q, r
    2. Find all omegas that are non-Shioda mod lcm(p, q, r)
    3. For each such combination of p, q, and r, and for each such omega, check_conjecture_2
    4. Record data
    '''
    if output_filepath and output_path_exists(output_filepath):
        with open(output_filepath, 'a') as outfile:
            header = 'p,q,r,omega,supersingularity,condition_1,condition_2,condition_3\n'
            outfile.write(header)
            print(header.strip())
            for p, q, r, omega in iterate_conjecture_inputs(minimum_prime, maximum_prime, minimum_omega, maximum_omega):
                supersingularity, condition_1, condition_2, condition_3 = check_conjecture_2(p, q, r, omega)
                line = '{},{},{},{},{},{},{},{}\n'.format(p, q, r, omega, supersingularity, condition_1, condition_2, condition_3)
                outfile.write(line)
                print(line.strip())
    else:
        print('p,q,r,omega,supersingularity,condition_1,condition_2,condition_3')
        for p, q, r, omega in iterate_conjecture_inputs(minimum_prime, maximum_prime, minimum_omega, maximum_omega):
            supersingularity, condition_1, condition_2, condition_3 = check_conjecture_2(p, q, r, omega)
            line = '{},{},{},{},{},{},{},{}'.format(p, q, r, omega, supersingularity, condition_1, condition_2, condition_3)
            print(line)

def range_exception_1(minimum_exponent, maximum_exponent, minimum_omega, maximum_omega, output_filepath=None):
    if output_filepath and output_path_exists(output_filepath):
        with open(output_filepath, 'a') as outfile:
            header = 'supersingularity,exception,p,q,r,omega\n'
            outfile.write(header)
            print(header.strip())
            for p, q, r, omega in iterate_conjecture_inputs(minimum_exponent, maximum_exponent, minimum_omega, maximum_omega, exponents_prime=False):
                result = check_exception_1(p, q, r, omega)
                if result is not None:
                    valid_p, valid_q, valid_r, valid_omega = result
                    line = 'True,True,{},{},{},{}\n'.format(valid_p, valid_q, valid_r, valid_omega)
                    outfile.write(line)
                    print(line.strip())
    else:
        print('supersingularity,exception,p,q,r,omega')
        for p, q, r, omega in iterate_conjecture_inputs(minimum_exponent, maximum_exponent, minimum_omega, maximum_omega, exponents_prime=False):
            result = check_exception_1(p, q, r, omega)
            if result is not None:
                valid_p, valid_q, valid_r, valid_omega = result
                line = 'True,True,{},{},{},{}'.format(valid_p, valid_q, valid_r, valid_omega)
                print(line)


def check_req_ab_primroot_condition(input_filepath='data/final/conjecture_1.primes_2-1000.omegas_2-1000.1532572380.txt'):
    with open(input_filepath) as infile:
        for index, line in enumerate(infile):
            try:
                if index == 0:
                    continue
                a, b, c, p, supersingular, c1, c2, c3, c4, c5 = line.split(',')
                supersingular = bool(supersingular)
                if supersingular:
                    a, b, c, p = list(map(int, [a, b, c, p]))
                    if a % c == 1 and b % c == 1:
                        if primitive_root(p, a) is False and primitive_root(p, b) is False:
                            print(a, b, c, p)
            except ValueError:
                continue

# Main

def main(conjectures=None, exceptions=None):
    if conjectures is None and exceptions is None:
        raise Exception('(ERROR) Final.py - Function main requires conjectures or exceptions must be defined')
    elif exceptions:
        main_exceptions()
    else:
        main_conjectures()

def main_conjectures():
    conjecture_number = Integer(input('Conjecture Number (1, 2) [1]: ').strip() or 1)
    minimum_prime = Integer(input('Minimum Prime [2]: ').strip() or 2)
    maximum_prime = Integer(input('Maximum Prime [1000]: ').strip() or 1000)
    minimum_omega = Integer(input('Minimum Order [2]: ').strip() or 2)
    maximum_omega = Integer(input('Maximum Order [1000]: ').strip() or 1000)
    default_output_filepath = 'data/final/conjecture_{}.primes_{}-{}.omegas_{}-{}.{}.txt'.format(conjecture_number, minimum_prime, maximum_prime, minimum_omega, maximum_omega, int(time()))
    output_filepath = input('Output Filepath (<filepath>, "None") [{}]: '.format(default_output_filepath))
    output_filepath = output_filepath if output_filepath.strip().lower() not in ['none', ''] else default_output_filepath
    globals()['range_conjecture_{}'.format(conjecture_number)](minimum_prime, maximum_prime, minimum_omega, maximum_omega, output_filepath)

def main_exceptions():
    conjecture_number = Integer(input('Conjecture Number (1, 2) [1]: ').strip() or 1)
    if conjecture_number == 2:
        raise NotImplementedError
    minimum_prime = Integer(input('Minimum Prime [2]: ').strip() or 2)
    maximum_prime = Integer(input('Maximum Prime [1000]: ').strip() or 1000)
    minimum_omega = Integer(input('Minimum Order [2]: ').strip() or 2)
    maximum_omega = Integer(input('Maximum Order [1000]: ').strip() or 1000)
    default_output_filepath = 'data/final/exception_{}.primes_{}-{}.omegas_{}-{}.{}.txt'.format(conjecture_number, minimum_prime, maximum_prime, minimum_omega, maximum_omega, int(time()))
    output_filepath = input('Output Filepath (<filepath>, "None") [{}]: '.format(default_output_filepath))
    output_filepath = output_filepath if output_filepath.strip().lower() not in ['none', ''] else default_output_filepath
    globals()['range_exception_{}'.format(conjecture_number)](minimum_prime, maximum_prime, minimum_omega, maximum_omega, output_filepath)

if __name__ == '__main__':
    check_req_ab_primroot_condition()
    #main(exceptions=False)
