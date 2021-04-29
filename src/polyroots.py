import jacobi as j 
from polynomial import Polynomial 
import random

x = Polynomial(1,0)

def stdize(f,p):
    if f == Polynomial(0): 
        return f
    else: 
        g = f
        for i in range(f.degree + 1):
            g[i] = f[i] % p
        return g

def monic(f, p): 
    if f == Polynomial(0): 
        return f
    else: 
        a = f[f.degree]
        g = f * j.invmod(a, p)
        return stdize(g, p)

def dividemod(f, g, p): 
    f = stdize(f, p)
    g = stdize(g, p)
    if g.degree < 0:
        raise Exception("division by 0")
    else: 
        quo = Polynomial(0)
        rem = f
        while rem.degree >= g.degree:
            deg = rem.degree
            tmp = g * (x ** (deg - g.degree))
            coeff = (j.invmod(tmp[deg],p) * rem[deg]) % p 
            quo += coeff * (x ** (deg - g.degree))
            tmp = tmp * coeff
            rem = stdize(rem - tmp, p)
        return (quo, rem)
    
def gcdmod(f, g, p):
    if g == Polynomial(0): 
        return monic(f, p)
    else: 
        h = dividemod(f, g, p)[1]
        return gcdmod(g, h, p)

def polymodexp(f, n, g, p): # Computes f^n (mod g) in Z/pZ[X]
    result = Polynomial(1)
    tmp = dividemod(f, g, p)[1] 
    while n:
        if n % 2 != 0:
            result = dividemod(result * tmp, g, p)[1]
        tmp = dividemod(tmp * tmp, g, p)[1]
        n = n // 2
    return result

def speedygcd(f, n, h, g, p): # gcd(a^n+b, g) when n is large
    f1 = polymodexp(f, n, g, p) + h
    return gcdmod(f1, g, p)

def polyrootsmod(f, p):
    f = speedygcd(x, p, -x, f, p)
    if f.degree <= 0:
        return ([], -1)
    elif f.degree == 1:
        solution = (-1 * f[0] * j.invmod(f[1], p)) % p
        return ([solution], -1)
    else: 
        fails = []
        while True:
            v = random.randint(0,p-1)
            if v not in fails:
                v = Polynomial(v)
                factor1 = speedygcd(x + v, p//2, Polynomial(-1), f, p)
                if factor1.degree <= 0 or factor1.degree == f.degree:
                    fails.append(v[0])
                    continue
                else:
                    break
        factor2 = dividemod(f, factor1, p)[0]
        roots = polyrootsmod(factor1, p)[0] + polyrootsmod(factor2, p)[0]
        return (roots, len(fails)+1)
