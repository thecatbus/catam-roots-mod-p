import numpy as np 
import jacobi as j
import random

x = np.poly1d([1,0])

def reducemod(f, p):
    coeffs = f.c 
    for i in range(len(coeffs)):
        coeffs[i] = coeffs[i] % p
    return np.poly1d(coeffs)

def normalise(f, p):
    a = f.c[0]
    ainv = j.modexp(a, p-2, p)
    return reducemod(int(ainv)*f, p)

def dividemod(f, g, p):
    b = g.c[0]
    binv = j.modexp(b, p-2, p)
    quo = np.poly1d([0])
    rem = reducemod(f,p)
    while rem.o >= g.o and rem != np.poly1d([0]):
        a = rem.c[0]
        m = rem.o - g.o
        tmp = int(a * binv) * (x ** m)
        rem = reducemod(rem - (tmp * g), p)
        quo = reducemod(quo + tmp, p)
    return((quo, rem))

def gcdmod(f, g, p):
    if g == np.poly1d([0]):
        return normalise(f,p)
    else:
        rem = dividemod(f, g, p)[1]
        return gcdmod(g, rem, p)

def polymodexp(f, n, g, p): # Computes f^n (mod g) in Z/pZ
    if n == 0: 
        return 1
    elif n % 2 == 0: 
        f1 = polymodexp(f, n//2, g, p)
        return dividemod(f1 ** 2, g, p)[1]
    else: 
        f1 = polymodexp(f, (n-1)//2, g, p)
        return dividemod(f * (f1 ** 2), g, p)[1]

def speedygcd(a, m, b , g, p): # gcd(a ** m + b, g)
    return gcdmod(polymodexp(a, m, g, p) + b, g, p)

def rootsmod(f, p): 
    roots = []
    f = speedygcd(x, p, -x, f, p)
    if f.o == 0:
        return("irreducible")
    elif f.o == 1: 
        soln = int((-1 * f.c[1] * j.modexp(f.c[0], p-2, p)) % p)
        roots.append(soln)
    else: 
        values = list(range(p))
        fails = []
        while True:
            v = random.choice(values)
            if v not in fails: 
                xplusv = np.poly1d([1, v])
                factor1 = speedygcd(xplusv, (p-1)//2, np.poly1d([-1]), f, p)
                if factor1.o == 0 or factor1 == f:
                    fails.append(v)
                    continue
                else: 
                    break
        print(len(fails))
        factor2 = dividemod(f, factor1, p)[0]
        roots = roots + rootsmod(factor1, p) + rootsmod(factor2, p)  
    return roots
