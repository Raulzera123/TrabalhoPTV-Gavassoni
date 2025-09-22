
import scipy.integrate as integrate
import numpy as np
from numpy import sqrt, sin, cos, pi
from sympy import *

# Criando os simbolos 
init_printing(use_unicode=True)

alpha = symbols('alpha')

# Parâmetros usados para os cálculos
Emean = 13000000
Gmean = 650000
c = 1.2
b = 0.22
h = 1.6236
s = 6
rho_k = 390
g = 9.81
q_c = 1.5
Area = b * h
I = (b * h**3)/12
R = 43.33
theta = np.deg2rad(75)

f = R * (1 - cos(theta))
print(f"Flecha do arco: {f}")

L = 2 * R * sin(theta)
print(f"Comprimento do arco: {L}")

# Cargas atuantes no arco
Sarco = 2 * R * 75 * np.pi / 180

rho_p = (rho_k * b * h * g * Sarco) / (L * 1000)
print(f"Carga de Peso Próprio: {rho_p}")

q = rho_p + (q_c * s)
print(f"Cargas totais externas atuantes: {q}")

# Criando as variáveis reais
Va = q * L / 2
print(f"Reação vertical em A: {Va}")

H = (q * L**2) / (8 * f)
print(f"Reação de empuxo: {H}")

x = 41.8536 - 43.33 * sin(alpha)
y = 43.33 * cos(alpha) - 11.2146
expanded_x = expand(x**2)
k = q * expanded_x / 2

equa_N1 = (Va - q * x) * sin(alpha)
expanded_equa_N1 = expand(equa_N1)

equa_Q1 = (Va - q * x) * cos(alpha)
expanded_equa_Q1 = expand(equa_Q1)

# Esforços internos do arco Real
Mreal = (Va * x - k - H * y)
print(f"Equação para momento fletor do arco real: {expand(Mreal)}")

Nreal = (-H * cos(alpha) - expanded_equa_N1)
print(f"Equação para esforço normal do arco real: {expand(Nreal)}")

Qreal = (-H * sin(alpha) + expanded_equa_Q1)
print(f"Equação do esforço cortante do arco real: {expand(Qreal)}")


# Criando as variáveis virtuais
Va1 = 0.5
H1 = L / (4 * f)

# Esforços internos no arco Virtual
Mvirtual = (Va1 * x) - (H1 * y)
print(f"Equação para momento fletor do arco com carga virtual: {expand(Mvirtual)}")

Qvirtual = (-H1 * sin(alpha)) + (Va1 * cos(alpha))
print(f"Equação para esforço cortante virtual do arco: {expand(Qvirtual)}")

Nvirtual = (-H1 * cos(alpha)) - (Va1 * sin(alpha))
print(f"Equação para esforço normal virtual do arco: {expand(Nvirtual)}")

# Multiplicando Real x Virtual dos esforços
Mrv = expand(Mreal) * expand(Mvirtual)
print(f"Equação de Mreal x Mvirtual: {Mrv}")

Qrv = expand(Qreal * Qvirtual)
print(f"Equação de Qreal x Qvirtual: {Qrv}")

Nrv = expand(Nreal * Nvirtual)
print(f"Equação de Nreal x Nvirtual: {Nrv}")

# Realizando a Integral definida
MM = integrate(Mrv, (alpha, 0, theta))
v_MM = (2 * MM * R) / (Emean * I)
print(f"Flecha máxima pelo momento fletor: {v_MM.evalf()}")

NN = integrate(Nrv, (alpha, 0, theta))
v_NN = (2 * NN * R) / (Emean * Area)
print(f"Flecha máxima pela normal: {v_NN.evalf()}")

QQ = integrate(Qrv, (alpha, 0, theta))
v_QQ = (2 * QQ * R * c) / (Gmean * Area)
print(f"Flecha máxima pelo esforço cortante: {v_QQ.evalf()}")

# Verificação das flechas
vlim = L/300
v_mnq = v_MM.evalf() + v_NN.evalf() + v_QQ.evalf()
print(f"Flecha limite e flecha provocada por todos os esforços: {vlim}, {v_mnq}")

if vlim > v_mnq:
   print("Satisfaz a verificação")
else:
    print("Não satisfaz a verificação")

if f/L < 0.5:
    print("O arco é abatido")
else:
    print("O arco não é abatido")
