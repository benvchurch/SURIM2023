#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

'''
Special Case for Dwork Polynomial
Litt/Perry Group
Columbia University Mathematics REU 2018
Chunying Huangdai <ch3236@barnard.edu>
'''

# Imports #

from sage.all import *
from sage.matrix.operation_table import OperationTable
from itertools import product


#Math#

class CayleyTable(OperationTable):
    def __init__(self, order):
        self.order = order
        self.GF = GF(self.order)

    def multable(self):
        return OperationTable(self.GF, operator.mul)

    def addtable(self):
        return OperationTable(self.GF, operator.add)

    def mult(self, a, b):
        return self.multable.matrix_of_variables()[a][b]

    def add(self, a, b):
        return self.addtable.matrix_of_variables()[a][b]

    def list_cubes(self, elements):
        S=[self.mult(self.mult(x,x),x) for x in elements]
        return S

    def X_cubic_polynomials_in_three_variables(self,S,a):
        n=0
        for (x,y,z) in [x for x in itertools.product(S,S,S)]:
            if add(add(x,y),z)== self.mult(a,a):
                n+=1
        return n

''' for the positive integer n, f(x) = x_1^n + ..... + x_n^n
    g(x) = lambda * n * x_1 * ........... * x_n'''

def variable(order,n):
    '''list all the pt x = (x_1, ..., x_n) in Fq * ...... * Fq'''
    return [list(tuple) for tuple in product(list(range(order)),repeat=n)]

def mult(mt,a,b):
    return mt[a][b]
    

def add(at,a,b):
    return at[a][b]

def n_power(mt,n,num):
    '''the n-th power of the element'''
    power = num
    for i in range(n-1):
        power = mult(mt,power,num)
    return power


def f(mt,at,n,x):
    plst = [n_power(mt,n,x_i) for x_i in x]
    '''x is a list: x = (x_1, ..., x_n)
        compute the value of f(x) at the point x'''
    f = 0
    for i in plst:
        f = add(at,f,i)
    return f


def g(mt,at,n,x):
    base = x[0]
    '''x is a list: x = (x_1, ..., x_n)
        compute the value of g(x) at the point x'''
    for i in x[1:]:
        base = mult(mt,base,i)
    g = base
    for j in range(n-1):
        g = add(at,g,base)
    return g
    

def lst_lambda(mt,at,n,x,order):
    '''x is a list: x = (x_1, ..., x_n)'''
    fx = f(mt,at,n,x)
    gx = g(mt,at,n,x)
    lst_lambda = [lbd for lbd in range(order) if fx==mult(mt,gx,lbd)]
    print((x,fx,gx))
    print(lst_lambda)
    return x,lst_lambda

def runthrulambda(order,n):
    '''[x, its lambda list] for x in the field'''
    mt = OperationTable(GF(order), operator.mul).table()
    at = OperationTable(GF(order), operator.add).table()
    return [lst_lambda(mt,at,n,x,order) for x in variable(order,n)]


def count_lambda(lbd,order,n):
    '''the length of the list of the points for a given lambda'''
    lbdlst = runthrulambda(order,n)
    ptlst_lbd = [lst[0] for lst in lbdlst if lbd in lst[1]] 
    return len(ptlst_lbd)

def all_pt(prime, exp, n):
    '''the number of points for a given lambda s.t fx=gx'''
    order = prime**exp
    mt = OperationTable(GF(order), operator.mul).table()
    at = OperationTable(GF(order), operator.add).table()
    return [(lbd,count_lambda(lbd,order,n)) for lbd in range(n)]

def fermat(prime,exp):
    '''the number of points for fermat cases'''
    order = prime**exp
    return count_lambda(0,order,4)

if __name__ == '__main__':
    print((all_pt(2,1,2)))
    ct = CayleyTable(9)
    S = ct.list_cubes()
    print(ct.X_cubic_polynomials_in_three_variables())
