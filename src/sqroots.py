import jacobi as j

def rootmod(a,p):
    if p % 4 == 3: 
        root = j.modexp(a, (p+1)//4, p) % p
    elif p % 8 == 5:
        if j.modexp(a, (p-1)//4, p) == 1:
            root = j.modexp(a, (p+3)//8, p) % p
        else:
            u = j.modexp(2, (p-1)//4, p)
            v = j.modexp(a, (p+3)//8, p)
            root = (u * v) % p 
    else: 
        # decompose p-1 = (2 ** alpha) * s
        s = p - 1 
        alpha = 0
        while s % 2 == 0:
            s = s // 2 
            alpha = alpha + 1
        order = 2 ** alpha
        # find non-residue n and generator b =n ** s of roots of unity
        n = 1
        while j.jacobi(n,p) == 1:
            n = n + 1
        b = j.modexp(n, s, p)
        # find binary expansion of r
        r = []
        num = j.modexp(a, s, p)
        i = 0
        while num != 1:
            if j.modexp(num, 2 ** (alpha - i - 2), p) == 1: 
                r.append(0)
            else: 
                r.append(1) 
                u = j.modexp(b, order - 2**(i+1), p)
                num = (u * num) % p
            i += 1
        exp = order 
        for i in range(len(r)):
            exp = exp - (2 ** i) * r[i]
        yinv = j.modexp(b, exp, p)
        z = j.modexp(a, (s+1)//2, p) 
        root = (z * yinv) % p
    return root
