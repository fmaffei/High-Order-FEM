"""
Equação para achar zero de funções definidas no intervalo de -1 a 1

Zeros encontrados usando o teorma de Bolzano

Funções implementadas, em () como chamar ela na função:

- Legendre (legendre)
- Chebyshev de primeiro tipo (chebyshev01)

Parâmetros de entrada:
- Ordem do polinômio
- Escolha do polinômio

Parâmetro se saída:
- Matriz [1] x [ordem] com valores das raizes

"""


def find_zeros(escolha, ordem):
    from Polynomials.poly_jacobi import jacobi
    import numpy as np
    import math as mt

    li = -1                                              # Limite inferior do intervalo
    ls = 1                                               # Limite superior do intervalo
    reso = 0.001                                         # Defini em quantos pontos são avaliados os polinômios
    precis = mt.pow(10, -15)                             # Precisão do zero

# Aqui definimos as função que podem ser usadas. Quando adicionar uma função, adicionar na definição a descrição

    def funcao(ind):
        legendre = jacobi(0, 0, ordem, ind)
        chebyshev01 = jacobi(-0.5, -0.5, ordem, ind)
        chebyshev02 = jacobi(0.5, 0.5, ordem, ind)
        sele = {'legendre': legendre, 'chebyshev01': chebyshev01, 'chebyshev02': chebyshev02}
        return sele[escolha]

    x = np.zeros((1, int((ls-li)/reso)+1))               # Vetor variavel independente
    y = np.zeros((1, int((ls-li)/reso)+1))               # Vetor resultado
    z = []                                               # Inicia vetor de zeros

    for i in range(0, int((ls-li)/reso)+1, 1):           # Preencher o vetor de plot
        x[0][i] = li + reso*i

    for i in range(0, int((ls-li)/reso)+1, 1):
        y[0][i] = funcao(x[0][i])

    for i in range(0, int((ls - li) / reso) + 1, 1):
        if mt.fabs(y[0][i]) < precis:
            z = np.append(z, x[0][i])

        elif 1 <= i <= int((ls - li) / reso)+1:
            if y[0][i - 1] < 0 and y[0][i] > 0:
                veri = [mt.fabs(y[0][i - 1]), mt.fabs(y[0][i])]
                error = max(veri) - min(veri)
                sup = x[0][i]
                inf = x[0][i - 1]

                while error > precis:
                    w = (sup+inf)/2
                    t = funcao(w)
                    veri = [mt.fabs(funcao(sup)), mt.fabs(funcao(inf))]
                    error = max(veri) - min(veri)

                    if t < 0:
                        inf = w

                    elif t > 0:
                        sup = w

                z = np.append(z, inf)

            elif y[0][i - 1] > 0 and y[0][i] < 0:
                veri = [mt.fabs(y[0][i - 1]), mt.fabs(y[0][i])]
                error = max(veri) - min(veri)
                sup = x[0][i]
                inf = x[0][i - 1]

                while error > precis:
                    w = (sup + inf) / 2
                    t = funcao(w)
                    veri = [mt.fabs(funcao(sup)), mt.fabs(funcao(inf))]
                    error = max(veri) - min(veri)

                    if t > 0:
                        inf = w

                    elif t < 0:
                        sup = w

                z = np.append(z, inf)
    return z
