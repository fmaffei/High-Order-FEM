"""
Programa para gerar zeros da função polinomial de Legendre. Nesse algorítimo sera usado a aprozimação de Francesco
Tricomi com erro de truncamento O(n^-4).
"""

import numpy as np
import math as mt

ordem = 4                      # Ordem do polinomio

A = np.zeros((1, ordem))       # Matriz para armazenas as raizes

for i in range(0, ordem, 1):   # Preencher a matriz em forma decrescente
    A[0][i] = (1 - (1/(8*pow(ordem, 2))) + (1/(8*pow(ordem, 3)))) * mt.cos(mt.pi*(4*(i+1) - 1)/(4*ordem + 2))
    if mt.fabs(A[0][i]) < mt.pow(10, -5):    # Filtro para eliminar valores muito pequenos
        A[0][i] = 0

print(A)
