"""Exact symbolic audit. Dependency: python-flint (pip install python-flint)."""
from flint import fmpz_mod_mpoly_ctx

ctx = fmpz_mod_mpoly_ctx.get(['c', 'a', 'b', 'x', 'y'], 2)
c, a, b, x, y = ctx.gens()
p = a*x**2 + (a*a+b*b+b)*x + (a+c*b)*y**2 + (c*b*b+a)*y
q = ((a+c*b)*x**4 + c*b*b*x*x + (a**4+b+(1+c)**2*b**4)*x
     + ((1+c)**2*a+c*b)*y**4 + c*a*a*y*y + (c*a**4+c*b**4+a)*y)
t = a*y+b*x
s = c*t
R = a*a+c*a*b+b*b
N = a**3+a*b*b+c*b**3
S = c*N
D = 1+R+S
U = s*s+(1+R)*s+R
A0 = a*a+b*b+b
C0 = a**4+b+(1+c)**2*b**4
assert N*R == a**5+c*a**4*b+(1+c)**2*a*b**4+c*b**5
assert N*R*a*a+(a+c*b)*N*N == c*c*b**4*N
assert A0*A0+a*N == b*b*(R+1)
assert N*R*A0+C0*N == b*N*(1+R+c*N)
assert (c*c*t**4+(R*R+R+c*N)*t*t
        + N*(1+R+c*N)*t+R*p*p+N*R*p+N*q) == 0
assert U*U+D*U+D*D*s == R*((c*p)**2+1)+S*R*(c*p+1)+c*S*q
matrix = [[c*a+1, c, c*c*b], [c*b, c*a, c*(1+b)], [1, c*b, c*a]]
det = (matrix[0][0]*(matrix[1][1]*matrix[2][2]+matrix[1][2]*matrix[2][1])
       +matrix[0][1]*(matrix[1][0]*matrix[2][2]+matrix[1][2]*matrix[2][0])
       +matrix[0][2]*(matrix[1][0]*matrix[2][1]+matrix[1][1]*matrix[2][0]))
assert det == c*c*D
print('All coefficient, polar-certificate, quartic and norm identities verified')
print('exactly in F_2[alpha,a,b,x,y].')
