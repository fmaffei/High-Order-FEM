"""
Código para testar fundamentos teoricos da análise pro elementos finitos com malha equiespaçada

Probema simulado unidimensional:

u^2-u=0

du/dx (0) = 1
du/dx (10) = 1

onde "u" é a variável dependendete e "x" a independente
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt

ci = 1                                  # Condição de contorno inicial
cf = 1                                  # Condição de contorno final
ls = 10                                 # Limite superior do domínio
li = 0                                  # Limite inferior do domínio
numero_elementos = 15                   # Número de elementos no domínio

ta_ele = (ls-li)/numero_elementos       # Tamanho do elemento padrão

# Criando matriz de rigidez "K" para o elemento

K_e = np.zeros((2, 2))

K_e[0][0] = (ta_ele/3) + (1/ta_ele)
K_e[0][1] = (-1/ta_ele)+(ta_ele/6)
K_e[1][0] = K_e[0][1]
K_e[1][1] = K_e[0][0]

# Criando matriz de rigidez global

K = np.zeros((numero_elementos+1, numero_elementos+1))

for i in range(0, numero_elementos, 1):
    for j in range(0, 2, 1):
        for n in range(0, 2, 1):
            K[j+i][n+i] = K_e[j][n]+K[j+i][n+i]

# Criando a matrix forçante

F = np.zeros((numero_elementos+1, 1))

# Para as condições de contorno impostas

F[0][0] = -1
F[numero_elementos][0] = 1

# Solucionando o sistema trilinear pelo método de Thomas
# 1 - Reduzir a linha diagonal mais a esquerda

for i in range(1, numero_elementos+1, 1):
    K[i][i] = K[i][i] - ((K[i][i-1]*K[i-1][i])/K[i-1][i-1])

# 1 - Recalcular a coluna do termo forçante

for i in range(1, numero_elementos+1, 1):
    F[i][0] = F[i][0] - ((F[i-1][0]*K[i][i-1])/K[i-1][i-1])

# 1 - Zerar diagonal mais a esquerda

for i in range(0, numero_elementos, 1):
    K[i+1][i] = 0

# Criando matriz solução
U = np.zeros((numero_elementos+1, 1))

U[numero_elementos][0] = F[numero_elementos][0]/K[numero_elementos][numero_elementos]

for i in range(numero_elementos-1, -1, -1):
    U[i][0] = (F[i][0] - (K[i][i+1]*U[i+1][0]))/K[i][i]

# Criar vetor para plotar

xn = np.zeros((1, numero_elementos+1))
passo = (ls-li)/numero_elementos
xn[0][0] = li
xn[0][numero_elementos] = ls

x = np.linspace(0, 10, 100)
y = -(np.exp(10-x)-np.exp(x))/(1+np.exp(10))

for i in range(1, numero_elementos, 1):
    xn[0][i] = passo*i

plt.grid()
plt.plot(xn[0, :], U[:, 0], label='Numeric')
plt.plot(x, y, label='Analytic')

plt.xlabel('X')
plt.ylabel('Y')

plt.legend()
plt.show()
