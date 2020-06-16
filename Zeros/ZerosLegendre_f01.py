"""
Forma - 01

Programa para gerar zeros da função polinomial de Legendre. Nesse algorítimo sera usado a aprozimação de Francesco
Tricomi com erro de truncamento O(n^-4).

Parâmetros de entrada:
- Ordem do polinômio

Parâmetro se saída:
- Matriz de raizes
"""


def zeros_legendre_f01(ordem_p):

    import numpy as np
    import math as mt

    # ordem = 4                     # Ordem do polinomio - ativar para verificação
    ordem = ordem_p
    a = np.zeros((1, ordem))       # Matriz para armazenas as raizes

    for i in range(0, ordem, 1):   # Preencher a matriz em forma decrescente
        a[0][i] = (1 - (1/(8*pow(ordem, 2))) + (1/(8*pow(ordem, 3)))) * mt.cos(mt.pi*(4*(i+1) - 1)/(4*ordem + 2))
        if mt.fabs(a[0][i]) < mt.pow(10, -5):    # Filtro para eliminar valores muito pequenos
            a[0][i] = 0

    # print(a)  # Ativar para verificação
    return a
