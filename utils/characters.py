
# This file was *autogenerated* from the file /home/ben/Documents/REU2018/utils/characters.sage

# Imports #

from sage.all import *

# Parameters #

_sage_const_0 = Integer(0)
_sage_const_1 = Integer(1)
_sage_const_2 = Integer(2)

# Functions #

dic_of_sums = {}
UCF = UniversalCyclotomicField()

def add_char(u, p):
    z = UCF.gen(p)
    return z**(u.trace())

# a is in prime subfield

def char(alpha, a, p, k):
    if Mod(a, p) == 1:
        return 1
    if Mod(a, p) == 0:
        return 0
    else:
        g = primitive_root(p)
        a, g, p = Integer(a), Integer(g), Integer(p)
        l = Mod(a, p).log(g)
        z = UCF.gen(denominator(alpha))
        return z**int(numerator(alpha) * l * (p**k - 1) / (p - 1))

def gauss_sum(alpha, p, k):
    if alpha in ZZ:
        return 0
    d = denominator(alpha)
    a = numerator(alpha)
    v = Mod(p, d).multiplicative_order()
    if (alpha, p, v) in dic_of_sums:
        if k % v != 0:
            print("ORDER ERROR")
        return -(-dic_of_sums[(alpha, p, v)])**int(k/v)
    else:
        if v == 1 and p == 2:
                return 0
        else:
            result = 0
            z = UCF.gen(denominator(alpha))
            w = z**numerator(alpha)
            if v == 1:
                u = GF(p)(primitive_root(p))
            else:
                u = GF(p**v).gen()
            char = w
            element = u
            for i in range(p**v-1):
                result = result + char * add_char(element, p)
                char = char * w
                element = element * u

            dic_of_sums[(alpha, p, v)] = result
            if k % v != 0:
                print("ORDER ERROR")
            return -(-result)**int(k/v)

def get_chars(p, f, k, powers, coeffs, r):
    chars_of_coeffs = list()
    for i in range(_sage_const_0, r + _sage_const_1):
        temp_char = list()
        d_i = gcd(powers[i], (p**(f*k) - _sage_const_1))
        for j in range(_sage_const_0, d_i):
            alpha_i = j/d_i
            temp_char.append(conjugate(char(alpha_i, coeffs[i], p, f*k)))
        chars_of_coeffs.append(temp_char)
    return chars_of_coeffs

def get_sums(p, f, k, powers, r):
    gauss_sums = list()
    for i in range(_sage_const_0, r + _sage_const_1 ):
        temp_sum = list()
        d_i = gcd(powers[i], (p**(f*k) - _sage_const_1))
        for j in range(_sage_const_0, d_i):
            alpha_i = j/d_i
            gs = gauss_sum(alpha_i, p, f*k)
            temp_sum.append(gs)
        gauss_sums.append(temp_sum)
    return gauss_sums

def get_sums_size_alpha(p, f, k, powers, r):
    gauss_sums = list()
    for i in range(_sage_const_0, r + _sage_const_1 ):
        temp_sum = list()
        d_i = gcd(powers[i], (p**(f*k) - _sage_const_1))
        for j in range(_sage_const_0, d_i):
            alpha_i = j/d_i
            gs = gauss_sum(alpha_i, p, f*k)
            temp_sum.append((gs, alpha_i, p**k))
        gauss_sums.append(temp_sum)
    return gauss_sums
