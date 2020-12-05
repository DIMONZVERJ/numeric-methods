# encoding: utf-8
import math
import matplotlib.pyplot as plt

file_log = open('explicit_sceme.log', 'w')


def fi(x):
    return (x + 1) ** 2 * math.sin(math.pi * x)


def psi_1(t):
    return 0


def psi_2(t):
    return 1 / 12


def explicit_scheme(a, h, teta, eps):
    n = int(1 / h)
    i = 0
    u = []
    while True:
        if i == 0:
            u.append([fi(x) for x in [h * j for j in range(n + 1)]])
            file_log.write(f'Значения сеточной функции на нулевом слое: {u[i]}\n')
        else:
            lambda_var = a * teta / h ** 2
            u.append([(1 - 2 * lambda_var) * u[i - 1][j] + lambda_var * (u[i - 1][j + 1] + u[i - 1][j - 1])
                      for j in range(1, n)])
            u[i].insert(0, psi_1(i * teta))
            u[i].append(psi_2(i * teta))
            max_diff = max(map(lambda x: abs(x[0] - x[1]), zip(u[i], u[i-1])))
            file_log.write(f'Значения сеточной функции на {i}-ом слое: {u[i]}\n')
            file_log.write(f'Максимальная разность между u[i] и u[i-1]: {max_diff}\n\n')
            if max_diff < eps:
                break
        i += 1
    file_log.write(f'Ответ: {u[i-1]}; Количество итераций: {i-1}\n')
    file_log.write('Расчёты закончены')
    return u


if __name__ == "__main__":
    eps = 10 ** -4
    a = 1
    h = 0.1
    teta = h ** 2 / (2 * a)
    file_log.write(f'epsilon = {eps}; a = {a}; h = {h}; teta = {h ** 2 / (2 * a)}\n')
    file_log.write('Начинаем расчёт')
    answer = explicit_scheme(a, h, teta, eps)
    x = [i * h for i in range(int(1/h)+1)]
    plt.plot(x, answer[-1])
    plt.show()
    file_log.close()
