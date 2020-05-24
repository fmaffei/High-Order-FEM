"""Programa para gerar bases para elementos finitos de alta ordem. As bases são geradas no intervalo de  -1 ate 1"""

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

# Base equiespaçada

B = np.zeros((1, ordem+1))                           # Matriz para gerar os pontos de zeros das bases

if ordem > 1:                                        # Preenchendo o vetor de zeros
    for i in range(0, ordem, 1):
        passo = 2 / ordem
        B[0][0] = li
        B[0][ordem] = ls
        B[0][0+i] = li + passo*i

elif ordem == 0:
    B[0][0] = 1

for k in range(0, int((ls-li)/reso)+1, 1):   # Valor de "x"
    for j in range(0, ordem+1, 1):           # Posição do zero
        for i in range(0, ordem+1, 1):       # Loop da produtória
            if i != j:
                lagrange[k][j] = lagrange[k][j]*(x[0][k]-B[0][i])/(B[0][j]-B[0][i])

print(B)

fig = plt.figure()
plt.grid()
plt.plot(x[0, :], lagrange[:, :])
plt.show()
