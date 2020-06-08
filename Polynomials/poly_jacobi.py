"""
Função que cria polinômios de Jacobi no intervalo de -1 até 1

Parâmetros de entrada:
- Ordem do polinômio
- Alpha
- Beta
- Posição de "x"

Parâmetro se saída:
- Valo "y" do polinômio

"""

import math as mt


def jacobi(alpha, beta, ordem, x):

    # Pesos
    if alpha == 0 and beta == 0:  # Peso do polinomio de Legendre
        w = 1

    elif alpha == -0.5 and beta == -0.5:  # Peso polinômio de Chebshew de primeiro tipo
        w = mt.pow(2, 2 * ordem) * mt.pow(mt.factorial(ordem), 2) / mt.factorial(2 * ordem)

    elif alpha == 0.5 and beta == 0.5:  # Peso polinômio de Chebshew de segundo tipo
        w = mt.pow(2, 2*ordem)*mt.factorial(ordem)*mt.factorial(ordem+1) / mt.factorial(2*ordem+1)

    # Polinômios
    if ordem == 0:
        b = 1
        return b

    elif ordem == 1:
        b = 0.5*(alpha-beta+((alpha+beta+2)*x*w))
        return b

    elif ordem > 1:
        a = 1
        b = 0.5*(alpha-beta+((alpha+beta+2)*x))

    # Calculando n+1
        for i in range(1, ordem, 1):
            a1 = 2 * (i + 1) * (i + alpha + beta + 1) * ((2 * i) + alpha + beta)
            a2 = ((2 * i) + alpha + beta + 1) * (mt.pow(alpha, 2) - mt.pow(beta, 2))
            a3 = ((2 * i) + alpha + beta) * ((2 * i) + alpha + beta + 1) * ((2 * i) + alpha + beta + 2)
            a4 = 2 * (i + alpha) * (i + beta) * ((2 * i) + alpha + beta + 2)

            ab1 = (a2 + (a3 * x)) * b
            ab2 = a4 * a

            c = (ab1 - ab2) / a1

            # Essa rotina desloca os elementos para calcular a próxima ordem de polinómio

            a = b
            b = c
            c = 0
    return b*w

