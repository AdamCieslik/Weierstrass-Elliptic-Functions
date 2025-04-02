# **Weierstrass Elliptic Functions in Python, SageMath, and Julia**  

This repository contains: 

- WeierstrassP.ipynb - A Jupyter Notebook demonstrating how to implement Weierstrass elliptic functions in Python.
- weierstrass.py - A Python library containing coded implementations of Weierstrass elliptic functions.


## **Implementation Details**  

The Jupyter Notebook illustrates how to compute the following Weierstrass elliptic functions:

- $\wp(z)$ – the Weierstrass $\wp$ function  
- $\wp'(z)$ – its first derivative  
- $\sigma(z)$ – the Weierstrass sigma function  
- $\zeta(z)$ – the Weierstrass zeta function  

Additionally, the inverse function $\wp^{-1}$ is introduced, along with essential quantities such as:  
- The half-periods $\omega_i$ derived from lattice invariants $g_2$, $g_3$.  
- Expressions for $g_2$, $g_3$ in terms of the half-periods.


The functions implemented in this notebook are equivalent to those available in *Wolfram Mathematica 14*. The code relies on standard mathematical references, including:  
- [DLMF – Digital Library of Mathematical Functions](https://dlmf.nist.gov/23)  
- [Wolfram Function Repository](https://functions.wolfram.com/EllipticFunctions)  
- [Carlson Symmetric Forms (Wikipedia)](https://en.wikipedia.org/wiki/Carlson_symmetric_form)  

## **Requirements (Python)**  
To run the notebook, the following libraries are required:  
- `mpmath`  
- `numpy`  

# To-Do List
- Code a continous expresion for $\ln{\frac{\sigma(x-y)}{\sigma(x+y)}}$
- Develop a SageMath implementation.
- Develop a Julia implementation.

## Feel Free to Contribute!  

---

**Thank you for your interest!**  
