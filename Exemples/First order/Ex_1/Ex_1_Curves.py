"""
Código para testar fundamentos teoricos da análise pro elementos finitos com malha equiespaçada

Probema simulado unidimensional:

u^2+u=0

du/dx (0) = 1
du/dx (10) = 1

onde "u" é a variável dependendete e "x" a independente
"""

import numpy as np
import matplotlib.pyplot as plt
import math as mt

plt.rc('text', usetex=True)
plt.rcParams['font.family'] = 'sans-serif'

ci = 1                                  # Condição de contorno inicial
cf = 1                                  # Condição de contorno final
ls = 10                                 # Limite superior do domínio
li = 0                                  # Limite inferior do domínio
pontos_erro = 10                        # Números de pontos usados para calcular o erro L2
passo_e = 4                             # Passo entre uma malha e outra
n_max = 12 + passo_e                    # Número máximo de elementos na malha
n_min = 4                               # Número mínimo de elementos na malha
n_erros = int((n_max-n_min)/passo_e)    # Número de elementos da matriz de erro

Z = np.zeros((1, n_erros))              # Matriz de erros

figure1, ax1 = plt.subplots()
figure2, ax2 = plt.subplots()

for v in range(n_min, n_max, passo_e):

    erro = 0                                 # Inicia a somatória erro
    numero_elementos = v                     # Número de elementos no domínio
    ta_ele = (ls - li) / numero_elementos    # Tamanho do elemento padrão

    # Criando matriz de rigidez "K" para o elemento

    K_e = np.zeros((2, 2))

    K_e[0][0] = (ta_ele / 3) + (1 / ta_ele)
    K_e[0][1] = (-1 / ta_ele) + (ta_ele / 6)
    K_e[1][0] = K_e[0][1]
    K_e[1][1] = K_e[0][0]

    # Criando matriz de rigidez global

    K = np.zeros((numero_elementos + 1, numero_elementos + 1))

    for i in range(0, numero_elementos, 1):
        for j in range(0, 2, 1):
            for n in range(0, 2, 1):
                K[j + i][n + i] = K_e[j][n] + K[j + i][n + i]

    # Criando a matrix forçante

    F = np.zeros((numero_elementos + 1, 1))

    # Para as condições de contorno impostas

    F[0][0] = -1
    F[numero_elementos][0] = 1

    # Solucionando o sistema trilinear pelo método de Thomas
    # 1 - Reduzir a linha diagonal mais a esquerda

    for i in range(1, numero_elementos + 1, 1):
        K[i][i] = K[i][i] - ((K[i][i - 1] * K[i - 1][i]) / K[i - 1][i - 1])

    # 1 - Recalcular a coluna do termo forçante

    for i in range(1, numero_elementos + 1, 1):
        F[i][0] = F[i][0] - ((F[i - 1][0] * K[i][i - 1]) / K[i - 1][i - 1])

    # 1 - Zerar diagonal mais a esquerda

    for i in range(0, numero_elementos, 1):
        K[i + 1][i] = 0

    # Criando matriz solução
    U = np.zeros((numero_elementos + 1, 1))

    U[numero_elementos][0] = F[numero_elementos][0] / K[numero_elementos][numero_elementos]

    for i in range(numero_elementos - 1, -1, -1):
        U[i][0] = (F[i][0] - (K[i][i + 1] * U[i + 1][0])) / K[i][i]

    # Criar vetor para plotar

    xn = np.zeros((1, numero_elementos + 1))
    passo = (ls - li) / numero_elementos
    xn[0][0] = li
    xn[0][numero_elementos] = ls

    for i in range(0, numero_elementos, 1):
        xn[0][i] = passo * i

        # Criando rotina para calcular os valores de y dentro de cada elemento e subtrair do valor analitico

        xe = np.linspace((li + passo * i), (li + passo * (i + 1)), int(1000/numero_elementos))
        yn = (U[i][0]*((xe-(li + passo * i))/((li + passo * i)-(li + passo * (i + 1))))) +\
             (U[i+1][0]*((xe-(li + passo * (i + 1)))/(-(li + passo * i)+(li + passo * (i + 1)))))
        ya = -(np.exp(10 - xe) - np.exp(xe)) / (1 + np.exp(10))
        erro_1 = np.fabs(ya)-np.fabs(yn)                # Primeiro a subtração da resposta exata menos o numérico
        erro_2 = np.dot(erro_1, np.transpose(erro_1))   # Quadrado de cada diferença e sua soma
        erro = erro + mt.sqrt(erro_2)                   # Soma essa diferença para os outros elementos

    plt.figure(1)
    figure1 = plt.plot(xn[0, :], U[:, 0])

    Z[0][int((v/passo_e)-1)] = erro

# Configração de plotagem (Está escrito para três valores de elementos, caso seja necessário mias elementos, mudar os
# as linhas pertinentes)

plt.figure(1)
plt.grid()
plt.xlabel('X', fontsize=14)
plt.ylabel('Y', fontsize=14)

x = np.linspace(0, 10, 100)
y = -(np.exp(10 - x) - np.exp(x)) / (1 + np.exp(10))
plt.plot(x, y)

plt.legend(('Number of elements - 4', 'Number of elements - 8', 'Number of elements - 12', 'Analytic'),
           shadow=True, fontsize=14)
props = dict(facecolor='white', edgecolor='black', boxstyle='round')
text = '\n'.join((
    r'$\mathrm{Equation:}$',
    r'$\frac{\mathrm{d}^2\Phi }{\mathrm{d} x^2}\frac{}{}-\Phi =0$',
    r'$\mathrm{Boundary\:Conditions:}$',
    r"${\Phi }'\left ( 0 \right )=1\:and\:{\Phi }'\left ( 10 \right )=1$"
    ))
plt.text(4.2, -0.9, text, linespacing=1.5, fontsize=15, bbox=props)
ax1.tick_params(axis='both', which='major', labelsize=15)

plt.figure(2)
plt.xlabel('Number of elements', fontsize=14)
plt.ylabel('L2', fontsize=14)
ax2.tick_params(axis='both', which='major', labelsize=15)
plt.plot(range(n_min, n_max, passo_e), Z[0, :])

plt.grid()
plt.show()
plt.legend()
