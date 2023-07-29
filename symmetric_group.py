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

n = 6
G = SymmetricGroup(n)
C = G.gen(0)

def g(P, k):
    num = 0
    for i in range(1, n+1):
        if P(i) < P( (C**k) (i) ):
            num += 1
    return num

def N(e_0, e_1, mu, p, m):
    n = 0
    f = Mod(p, m).multiplicative_order()
    for i in range(0, f):
        print(mu * e_0 * p**i % m)
        if mu * (e_0 + e_1) * p**i % m < mu * e_0 * p**i % m:
            n += 1
    return n


p = 17
m = 5*7
f = Mod(p, m).multiplicative_order()
n = 0
s = 0
count = 0

'''
print "p = {}, m = {}, f = {}".format(p, m, f)
if f % 2 == 0:
    print "p^(f/2) ~ {} mod {}".format(Mod(p**(f/2), m), m)
'''

f = 10
M = 60
'''
print "M = {} f = {}"
for j in range(100):
    l = []
    for i in range(f):
        r = random.randint(1, M - 1)
        while r in l:
            r = random.randint(1, M - 1)
        l.append(r)

    print sorted(l), Integer(sum(l))/Integer(M)
    if Integer(sum(l))/Integer(M) == Integer(f)/Integer(2):
        print "FOUND"
'''

for j in range(1, M):
    l = set()
    for i in range(euler_phi(M)):
        l.add(Integer(Mod(j**i, M)))
    f = len(l)
    if Integer(sum(l))/Integer(M) == Integer(f)/Integer(2):
        print("FOUND")
        print(sorted(list(l)), Integer(sum(l))/Integer(M))

'''
print N(1, 19, 1, p, m)
'''
'''
for i in range(1, m):
    for j in range(1, m):
        if (i + j) % m != 0:
            num = N(i, j, 1, p, m)
            print i,j, (i + j) % m, num
            n += num
            s += num**2
            count += 1
    print "i = {}, Average = {}, Variance = {}".format(i, n * 1.0 / count, s * 1.0 / count - (n * 1.0 / count)**2)
'''
'''
for P in G:
    for k in range(1, n):
        print k, P, g(P, k),  g(P, n-k), g(P, k) +  g(P, n-k)
        if g(P, k) +  g(P, n-k) != n:
            print "FUCK FUCK FUCK"
            input()
'''
