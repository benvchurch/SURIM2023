from sage.all import *

p = 7
v = 1
k = GF(p**v, name='a', modulus="primitive");

a = k.gen()

def multi_order(b):
    if b == a:
        return 1
    elif b == 0:
        return 0
    else:
        return 1 + multi_order(b/a)


def multi_char(alpha, u):
    z = CyclotomicField(denominator(alpha)*p).gen()**p
    if u == 0:
        if alpha == 0:
            return 1
        else:
            return 0
    else:
        return z**(numerator(alpha)*multi_order(u))


def add_char(u,p,v,alpha):
    z = CyclotomicField(p*denominator(alpha)).gen()**denominator(alpha)
    return z**(list(u.polynomial())[0])


def gauss_sum_brute_force(p,v, alpha):
    result = 0
    z = CyclotomicField(denominator(alpha)).gen()
    w = z**numerator(alpha)
    r = w
    for i in k:
        if i !=0:
            result = result + multi_char(alpha, i) * add_char(i,p,v,alpha)

    return result


print(N(gauss_sum_brute_force(p,v,Integer(1)/Integer(3))))
