import math as m
L = 0.12
a = 0.004
k = a * m.pi / L
H = 0.0025
W = 0.02
E = 4200000
V = 0.3
Bcr = 0.006326
mu0 = 4*m.pi*1e-7

B = (mu0 * E * m.pi * k**3) / (24 * 3**0.5 * Bcr)
B = B - (mu0 * E * m.pi**3 * k**3 * H**2) / (72 * 3**0.5 * L**2 * (1 - V**2) * Bcr)
L0 = 3 * L**3 * (1 + 0.25 *k**2) * (1-0.3**2) / (3*(1-0.3**2)*L**2-H**2*m.pi**3)


a = 7
L = 70
d = 0.1

Bs = 37.5
B = m.pi * Bs / (m.atan(a**2/(2*d*(2*a**2+4*d**2)**0.5)) - m.atan(a**2/(2*(d+L)*(2*a**2+4*(d+L)**2)**0.5)))
mr = 6.42
Rho = mr * 1e6 / (a * a * L)

'''
B = 78
Bs = B * (m.atan(a**2/(2*d*(2*a**2+4*d**2)**0.5)) - m.atan(a**2/(2*(d+L)*(2*a**2+4*(d+L)**2)**0.5))) / m.pi
'''

