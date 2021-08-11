# -*- coding: utf-8 -*-
"""
Programa que permite calcular el factor de compresibilidad z basado en la
gráfica de Standing and katz. 
Basado en la presion y temperatura pseudo reducida


Created on Tue Aug  3 13:48:55 2021

@author: Carelia
"""
import math
import numpy as np

# Reduced temperature
# T and Tc have the units of °R.
def reduced_temperature(temperature, tc):
    tr = temperature / tc
    return tr

# Reduced pressure
# P and Pc have the units of psia
def reduced_pressure(pressure, pc):
    pr = pressure / pc
    return pr

# Redlich-Kwong constant A
def redlich_kwong_constant_a(pr, tr):
    a = 0.42748 * (pr / tr ** 2.5)
    return a

# Redlich-Kwong constant B
def redlich_kwong_constant_b(pr, tr):
    b = 0.08664 * (pr / tr)
    return b

# Cubic constant α
def cubic_constant_alpha(a, b):
    alpha = (1 / 3) * (3 * (a - b - b ** 2) - 1)
    return alpha

# Cubic constant β
def cubic_constant_beta(a, b):
    beta = (1 / 27) * (-2 + (9 * (a - b - b ** 2)) - (27 * a * b))
    return beta

# Discriminant
def discriminant(alpha, beta):
    d = (beta ** 2 / 4) + (alpha ** 3 / 27)
    return d

# Solution Constant A*
def solution_constant_a_star(beta, d):
    a_star = np.cbrt((-beta / 2) + np.sqrt(d))
    return a_star


# Solution Constant B*
def solution_constant_b_star(beta, d):
    b_star = np.cbrt((-beta / 2) - np.sqrt(d))
    return b_star


# Solution Constant Theta
def solution_constant_theta(beta, alpha):
    theta = math.acos(-sign(beta) * (math.sqrt((beta ** 2 / 4) / (-alpha ** 3 / 27))))
    return theta

# SIGN(β)
# If β > 0, SIGN(β) = 1
# If β = 0, SIGN(β) = 0
# If β < 0, SIGN(β) = -1
def sign(beta):
    if beta > 0:
        sign_beta = 1
    elif beta < 0:
        sign_beta = -1
    else:
        sign_beta = 0
    return sign_beta


# Trial Root Z1
def trial_root_z1(a_star, b_star):
    z1 = a_star + b_star + 1 / 3
    return z1


# Trial Root Zi
# for i = 2, 3
def trial_root_zi(a_star, b_star, i):
    zi = (-(1 / 2) * (a_star + b_star)) + ((1 / 3) * i)
    return zi


# Trial Root Zi+1
# for i = 1, 2, 3
def trial_root_zi1(alpha, theta, i):
    zi1 = (2 * math.sqrt(- alpha / 3) * math.cos((theta / 3) + (i * ((math.pi * 2) / 3)))) + (1 / 3)
    return zi1


def compressibility_factor_calc(temperature, tc, pressure, pc):
    """
    Calculates compressibility factor
    
    """
  
    tr = reduced_temperature(temperature, tc)
    pr = reduced_pressure(pressure, pc)
    a = redlich_kwong_constant_a(pr, tr)
    b = redlich_kwong_constant_b(pr, tr)
    alpha = cubic_constant_alpha(a, b)
    beta = cubic_constant_beta(a, b)
    d = discriminant(alpha, beta)
    if d < 0:
        theta = solution_constant_theta(beta, alpha)
        z1 = trial_root_zi1(alpha, theta, 1)
        z2 = trial_root_zi1(alpha, theta, 2)
        z3 = trial_root_zi1(alpha, theta, 3)
    else:
        a_star = solution_constant_a_star(beta, d)
        b_star = solution_constant_b_star(beta, d)
        if d == 0:
            z1 = trial_root_z1(a_star, b_star)
            z2 = trial_root_zi(a_star, b_star, 2)
            z3 = trial_root_zi(a_star, b_star, 3)
        else:
            z = trial_root_z1(a_star, b_star)
            return z
    z = max(z1, z2, z3)
    return z




print('Este programa calcula el factor z basado en Standing & Katz')

temperature=float(input('Ingrese la temperatura del yacimiento en R  '))
tc=float(input('Ingrese la temperatura pseudocritica en R  '))
pressure=float(input('Ingrese la presion del yacimiento en psia  '))
pc=float(input('Ingrese la presion pseudocritica en psia  '))

z=compressibility_factor_calc(temperature, tc, pressure, pc)
print('\n El factor de compresibilidad es ',round(z,3))

