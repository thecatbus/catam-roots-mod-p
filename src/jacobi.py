def modexp(a,n,p): 
    if n == 0:
        return 1
    elif n%2==0:
        tmp = modexp(a,n // 2,p) 
        val = (tmp * tmp)%p
        if val <= (p-1) // 2:
            return val 
        else: 
            return val - p
    else:
        tmp = modexp(a, (n-1) // 2, p)
        val = (a * tmp * tmp)%p
        if val <= (p-1) // 2:
            return val 
        else: 
            return val - p

def legendre(a,p):
    return modexp(a, (p-1) // 2, p)

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
