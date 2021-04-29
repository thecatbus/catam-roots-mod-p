def modexp(a,n,p): 
    result = 1
    x = a % p 
    while n:
        if n % 2 != 0:
            result = (result * x) % p
        x = (x * x) % p
        n = n // 2
    return result

def invmod(a, p): 
    return modexp(a, p-2, p) 

def legendre(a,p):
    result = modexp(a, (p-1) // 2, p)
    if result <= p//2:
        return result
    else:
        return (result-p)

def jacobi(m,n):
    m = m % n
    if n == 1:
        return 1
    elif m == 0:
        return 0
    else: 
        sgn = 1
        tmp = n % 8
        while m % 2 == 0:
            m = m // 2
            if tmp == 3 or tmp == 5:
                sgn = - sgn
        if m % 4 == 3 and n % 4 == 3: 
            sgn = - sgn
        return sgn * jacobi(n,m)
