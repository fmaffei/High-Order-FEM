"""
Programa para gerar bases para elementos finitos de alta ordem. As bases são geradas no intervalo de  -1 ate 1

Esse código é referente ao polinómio de Jacobi com alpha e beta igual a menos meio, ou seja, Chebyshev de primeiro tipo
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt

ordem = 5                                            # Ordem do polinómio
li = -1                                              # Limite inferior do intervalo
ls = 1                                               # Limite superior do intervalo
reso = 0.01                                          # Resolução da curva
alpha = -0.5                                            # Valores para polinomio de Legendre
beta = -0.5                                             # Valores para polinomio de Legendre

lagrange = np.ones((int((ls-li)/reso)+1, ordem+1))   # Iniciando a produtória
x = np.zeros((1, int((ls-li)/reso)+1))               # Vetor para plotar a base

for i in range(0, int((ls-li)/reso)+1, 1):           # Preencher o vetor de plot
    x[0][i] = li + reso*i

A = np.zeros((int((ls-li)/reso)+1, 1))               # Cria uma matriz para n-1
B = np.zeros((int((ls-li)/reso)+1, 1))               # Cria uma matriz para n
C = np.zeros((int((ls-li)/reso)+1, 1))               # Cria uma matriz para n+1

if ordem == 0:
    B = np.ones((int((ls-li)/reso)+1, 1))

elif ordem == 1:
    for i in range(0, int((ls-li)/reso)+1, 1):
        j = 1
        w = mt.pow(2, 2*j)*mt.pow(mt.factorial(j), 2)/mt.factorial(2*j)
        B[i][0] = 0.5*(alpha-beta+((alpha+beta+2)*x[0][i]))*w

elif ordem > 1:
    for i in range(0, int((ls-li)/reso)+1, 1):       # Preenchendo matriz para n-1
        A[i][0] = 1
    for i in range(0, int((ls-li)/reso)+1, 1):       # Preenchendo matriz para n
        B[i][0] = 0.5*(alpha-beta+((alpha+beta+2)*x[0][i]))

# Rotina para calcular n+1

    for i in range(1, ordem, 1):
        for j in range(0, int((ls-li)/reso)+1, 1):
            a1 = 2*(i+1)*(i+alpha+beta+1)*((2*i)+alpha+beta)
            a2 = ((2*i)+alpha+beta+1)*(mt.pow(alpha, 2)-mt.pow(beta, 2))
            a3 = ((2*i)+alpha+beta)*((2*i)+alpha+beta+1)*((2*i)+alpha+beta+2)
            a4 = 2*(i+alpha)*(i+beta)*((2*i)+alpha+beta+2)

            ab1 = (a2+(a3*x[0][j]))*B[j][0]
            ab2 = a4*A[j][0]

            C[j][0] = (ab1-ab2)/a1

# Essa rotina desloca os elementos para calcular a próxima ordem de polinómio

        for k in range(0, int((ls-li)/reso)+1, 1):
            A[k][0] = B[k][0]
        for k in range(0, int((ls - li) / reso) + 1, 1):
            B[k][0] = C[k][0]
        C = np.zeros((int((ls-li)/reso)+1, 1))

# Aplicação do peso correspondente para criar o polinómio de Chebtshev de primeiro tipo

w = mt.pow(2, 2*ordem)*mt.pow(mt.factorial(ordem), 2) / mt.factorial(2*ordem)
B = B*w

# Testando função polinomio de Jacobi. Essa forma é equivalente a escrita nas linhas anteriores.
Q = np.zeros((1, int((ls-li)/reso)+1))
for i in range(0, int((ls-li)/reso)+1, 1):
    from Polynomials.poly_jacobi import jacobi
    Q[0][i] = jacobi(alpha, beta, ordem, x[0][i])

# Rotina para plottar
fig1 = plt.figure()
plt.grid()
plt.plot(x[0, :], B[:, 0])

fig2 = plt.figure()
plt.grid()
plt.plot(x[0, :], Q[0, :])

plt.show()
