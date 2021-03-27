import jacobi as j
import sqroots as s
import random

# Pretty printing
def monostr(number, digitcount):
    n = number
    if n == 0:
        count = 1
    else:
        count = 0
        while(n > 0):
            n = n // 10
            count = count + 1 
    return (" "*(digitcount - count) + str(number))

def legstr(n):
    if n == -1:
        return ("-1")
    else:
        return (" "+str(n))

# Question 1
p1 = 30275233
f1res = 0
f2res = 0

# clear contents before proceeding
f1 = open("../output/Legendre-100.txt", "w"); f1.close()
f2 = open("../output/Legendre-random.txt", "w"); f2.close()

f1 = open("../output/Legendre-100.txt", "a")
f2 = open("../output/Legendre-random.txt", "a")
for i in range(100):
    if i % 5 == 0 and i > 0:
        f1.write("\n")
    if i % 4 == 0 and i > 0:
        f2.write("\n")
    leg1 = j.legendre(i+1,p1)
    if leg1 == 1:
        f1res += 1
    f1.write(monostr(i+1, 3) + ": " + legstr(leg1) + " "*5)
    rnum = random.choice(range(p1)) + 1
    leg2 = j.legendre(rnum,p1)
    if leg2 == 1:
        f2res += 1
    f2.write(monostr(rnum, 8) + ": " + legstr(leg2) + " "*5)
f1.write("\n\nNumber of quadratic residues encountered: "+ str(f1res))
f2.write("\n\nNumber of quadratic residues encountered: "+ str(f2res))
f1.close(); f2.close()

# Question 2
# clear contents before proceeding
f = open("../output/Test-roots.txt", "w"); f.close()

def gentests(prime, count): 
    i = 0
    while i < count: 
        a = random.choice(range(prime))
        if j.jacobi(a, prime) == 1:
            f.write("("+str(a) 
                    + ", " 
                    + str(s.rootmod(a,prime)) + ")  ") 
            i += 1
        else:
            continue 
    f.write("\n")

f = open("../output/Test-roots.txt", "a")

f.write("mod 65537:   ")
f.write("(18612, " + str(s.rootmod(18612,65537)) + ")  ")
gentests(65537, 2)
f.write("mod 10501:   ")
gentests(10501, 3)
f.write("mod 10601:   ")
gentests(10601, 3)
f.write("mod 11311:   ")
gentests(11311, 3)
f.write("mod 11411:   ")
gentests(11411, 3)

f.close()

f = open("../output/Twenty-roots.txt", "w"); f.close()
f = open("../output/Twenty-roots.txt", "a")
p = 30275233

count = 0
for a in range(20):
    if j.jacobi(a, p) != 1:
        continue
    else: 
        count += 1
        f.write("("+str(a)+", "+str(s.rootmod(a,p))+")   ")
        if count % 5 == 0:
            count = 0
            f.write("\n")
