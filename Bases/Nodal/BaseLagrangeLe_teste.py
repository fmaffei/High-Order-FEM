"""
Programa para gerar bases para elementos finitos de alta ordem. As bases são geradas no intervalo de  -1 ate 1

Base de Lagrange com zeros posicionados nos valores de zeros deu um polinômio de Legendre
"""

from Zeros.ZerosLegendre_f01 import zeros_legendre_f01
import numpy as np
import matplotlib.pyplot as plt

ordem = 5                                            # Ordem do polinómio
li = -1                                              # Limite inferior do intervalo
ls = 1                                               # Limite superior do intervalo
reso = 0.01                                          # Resolução da curva

lagrange = np.ones((int((ls-li)/reso)+1, ordem+1))   # Iniciando a produtória
x = np.zeros((1, int((ls-li)/reso)+1))               # Vetor para plotar a base

for i in range(0, int((ls-li)/reso)+1, 1):           # Preencher o vetor de plot
    x[0][i] = li + reso*i

if ordem > 1:
    C = -1*zeros_legendre_f01(ordem)                 # Matriz para gerar os pontos de zeros das bases
    print(C)
    B = np.zeros((1, ordem+1))
    B[0][0] = li
    B[0][ordem] = ls
    for z in range(1, ordem, 1):
        B[0][z] = C[0][z-1]
    print(B)

elif ordem == 1:
    B = np.zeros((1, 2))
    B[0][0] = li
    B[0][1] = ls

elif ordem == 0:
    B = np.zeros((1, 1))
    B[0][0] = 1

for k in range(0, int((ls-li)/reso)+1, 1):   # Valor de "x"
    for j in range(0, ordem+1, 1):             # Posição do zero
        for i in range(0, ordem+1, 1):         # Loop da produtória
            if i != j:
                lagrange[k][j] = lagrange[k][j]*(x[0][k]-B[0][i])/(B[0][j]-B[0][i])

fig = plt.figure()
plt.grid()
plt.plot(x[0, :], lagrange[:, :])
plt.show()
