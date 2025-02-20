from sympy import *
a0, a1, a2, a3, a4 = symbols('a0, a1, a2, a3, a4')
EI, h, x, q = symbols('EI, h, x, q')

w = a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4

wp = diff(w, x)
wpp = diff(wp, x)
wppp = diff(wpp, x)
wpppp = diff(wppp, x)

eq = Eq(EI * wpppp - q, 0)
sol = solve(eq, a4)
w = simplify(w.subs(a4, q/(24*EI)))
# w.subs(a4, sol[0])
print(simplify(w))
wp = diff(w, x)
wpp = diff(wp, x)
sol = solve([Eq(w.subs(x, -h/2), 0), Eq(wp.subs(x, -h/2), 0), Eq(w.subs(x, h/2), 0), Eq(wp.subs(x, h/2), 0)], [a0, a1, a2, a3])
print(sol)

print(w.subs(sol))

M = diff(w.subs(sol), x, 2)
print(M)