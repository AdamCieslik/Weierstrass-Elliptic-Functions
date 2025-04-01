"""
weierstrass.py

Author: Adam Cieślik
Date: February 2025
Version: 1.0

Description:
This library provides Weierstrass elliptic functions and related utilities, namely:
- Omega1(g2, g3)
- Omega2(g2, g3)
- Omega3(g2, g3)
- WeierstrassInvariants(omega_i, omega_j)
- WeierstrassP(z, g2, g3)
- reduce_argument(z, g2,g3 )
- InverseWeierstrassP(z, g2, g3) (to be updated)
- WeierstrassPPrime(z, g2, g3)
- WeierstrassSigma(z, g2, g3)
- WeierstrassZeta(z, g2, g3)

A detailed description can be found at:  
[GitHub Repository](https://github.com/AdamCieslik/Weierstrass-Elliptic-Functions)

Dependencies:
- **mpmath** (for arbitrary precision arithmetic)
- **numpy** (for polynomial root finding and complex number operations)

### Usage example:
```python
import weierstrass
omega1 = weierstrass.Omega1(1.234, -4.213)
print("Omega1 =", omega1)

License:
This code is provided under the MIT License. Feel free to modify and distribute.
"""

from mpmath import *
import numpy as np

# Precision settings for mpmath
mp.dps = 25
mp.pretty = True


# wp function

## Half-periods

def Omega1(g2, g3):
    # Coefficients of the polynomial 4*x**3 - g2*x - g3 = 0
    coeffs = [4, 0, -g2, -g3]

    # Finding the roots
    roots = np.roots(coeffs)
    
    # Sorting the roots
    roots = np.sort_complex(roots)
    e1, e2, e3 = roots

    # Calculation of parameters k^2 and k'^2
    ksqr = (e2 - e3) / (e1 - e3)
    ksqrp = (e1 - e2) / (e1 - e3)


    # Calculating omega1 and omega3
    omega1 = (1 / sqrt(e1 - e3)) * ellipk(ksqr)
    omega3 = (1 / sqrt(e3 - e1)) * ellipk( ksqrp)

    # Calculation of the Im(tau)
    Imtau = (omega3 / omega1).imag

    # Condition for selecting omega1 sign
    res = omega1 if Imtau > 0 else -omega1

    return res

def Omega3(g2,g3):
    # Coefficients of the polynomial 4*x**3 - g2*x - g3 = 0
    coeffs = [4, 0, -g2, -g3]

    # Finding the roots
    roots = np.roots(coeffs)
    
    # Sorting the roots
    roots = np.sort_complex(roots)
    e1, e2, e3 = roots

    # Calculation of parameter k'^2
    ksqrp = (e1-e2)/(e1-e3)

    # Calculating omega3
    omega3= (1/sqrt(e3-e1))*ellipk(ksqrp)

    return omega3

def Omega2(g2,g3):
    res = -Omega1(g2,g3) - Omega3(g2,g3)
    return res

## Lattice invariants

def WeierstrassInvariants(omega1, omega3):
    ttau = omega3 / omega1
    if ttau.imag<0:
        print("Im(tau)<0, change Omega1 or Omega3")
        
    q = qfrom(tau=ttau)

    # Calculation of g2 and g3
    g2 = (pi**4 / (12 * omega1**4)) * (jtheta(2, 0, q)**8 - jtheta(3, 0, q)**4 * jtheta(2, 0, q)**4 + jtheta(3, 0, q)**8)
    
    g3 = (pi / (2 * omega1))**6 * ((8/27)*(jtheta(2, 0, q)**12 + jtheta(3, 0, q)**12)-(4/9)*(jtheta(2, 0, q)**4 +jtheta(3, 0, q)**4) * jtheta(2, 0, q)**4 * jtheta(3, 0, q)**4)

    return (g2, g3)

## Reduce argument

def reduce_argument(tilde_z, g2,g3 ):
    
    oomega1 = Omega1(g2, g3) 
    oomega3 = Omega3(g2,g3)
    ttau = oomega3/oomega1 
    # Step 1: Convert tilde_z to c coordinates: c = omega^{-1} * tilde_z
    c = (1/oomega1) * tilde_z 
    
    # Step 2: Separating c into its imaginary part to find beta
    Im_tau = ttau.imag
    Im_c = c.imag
    beta = (1/Im_tau) * Im_c  # beta ∈ R

    # Step 3: Calculation of alpha
    Re_tau = ttau.real
    Re_c =  c.real
    alpha = Re_c - Re_tau * beta  # alpha ∈ R
    
    # Step 4: Reducing the coordinates of mu and nu to the interval [0,1)
    mu = alpha/2- floor(alpha/2)
    nu = beta/2 - floor(beta/2)

    # Step 5: Reconstruction of the reduced z
    z_red =2* oomega1 * (mu + ttau * nu)
    
    return z_red

## Wp

def WeierstrassP(z,g2,g3):
    q = qfrom(tau=Omega3(g2,g3)/Omega1(g2,g3))
    FirstCase = (pi/(2*Omega1(g2,g3))) *jtheta(3, 0,q) *jtheta(4, 0, q) * jtheta(2, pi*z/(2*Omega1(g2,g3)),q)/jtheta(1, pi*z/(2*Omega1(g2,g3)),q)
    SecondCase = (1/12)*(pi/Omega1(g2,g3))**2 * (jtheta(2, 0,q)**4 + 2*jtheta(4, 0,q)**4)
    res = FirstCase**2 + SecondCase
    return res

## Inverse Wp (to update!)

def InverseWeierstrassP(z,g2,g3):
    # Coefficients of the polynomial 4*x**3 - g2*x - g3 = 0
    coeffs = [4, 0, -g2, -g3]

    # Finding the roots
    roots = np.roots(coeffs)
    
    # Sorting the roots
    roots = np.sort_complex(roots)
    e1, e2, e3 = roots
    
    return elliprf(z-e1, z-e2, z-e3)

## First derivative of the Wp

def WeierstrassPPrime(z,g2,g3):
    q = qfrom(tau=Omega3(g2,g3)/Omega1(g2,g3))
    numerator = jtheta(2, pi*z/(2*Omega1(g2,g3)),q)*jtheta(3, pi*z/(2*Omega1(g2,g3)),q)*jtheta(4, pi*z/(2*Omega1(g2,g3)),q)*(jtheta(1, 0, q, derivative=1)**3)
    denominator = jtheta(2, 0,q)*jtheta(3, 0,q)*jtheta(4, 0,q)*(jtheta(1, pi*z/(2*Omega1(g2,g3)),q)**3)
    res = -(1/4)*(pi/Omega1(g2,g3))**3 * (numerator/denominator)
    return res

# Sigma function

def eta1(g2,g3):
    q = qfrom(tau=Omega3(g2,g3)/Omega1(g2,g3))
    res = - (pi**2/(12*Omega1(g2,g3))) * (jtheta(1, 0, q, derivative=3)/jtheta(1, 0, q, derivative=1))
    return res

def WeierstrassSigma(z,g2,g3):
    q = qfrom(tau=Omega3(g2,g3)/Omega1(g2,g3))
    Exp = exp((eta1(g2, g3) * z**2) / (2 * Omega1(g2, g3)))
    res = 2*Omega1(g2,g3)*Exp * (jtheta(1, pi * z/(2*Omega1(g2,g3)),q)/(pi*jtheta(1, 0, q, derivative=1)))
    return res

# Zeta function

def WeierstrassZeta(z,g2,g3):
    q = qfrom(tau=Omega3(g2,g3)/Omega1(g2,g3))
    firstP = eta1(g2,g3) * z/Omega1(g2,g3)
    secondP =  pi/(2*Omega1(g2,g3)) *jtheta(1,  pi * z/(2*Omega1(g2,g3)), q, derivative=1)/jtheta(1,  pi * z/(2*Omega1(g2,g3)), q)
    return firstP + secondP