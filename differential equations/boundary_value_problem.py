# encoding: utf-8

from common_functions import thomas_algorithm

alpha0 = 1
alpha1 = 0
A = 1.2
a = 0.4
beta0 = 1
beta1 = 2
B = 1.4
b = 0.7
n = 6


def p(x):
    return -0.5 * x


def q(x):
    return 1


def f(x):
    return 2


if __name__ == "__main__":
    h = (b - a) / n
    a_val = [None]
    b_val = [alpha0 - alpha1 / h]
    c_val = [alpha1 / h]
    d_val = [A]
    for i in range(2, n):
        x_i = a + h * i
        a_val.append((1 - 0.5 * p(x_i) * h) / (1 + 0.5 * p(x_i) * h))
        b_val.append(-(2 - q(x_i) * h ** 2) / (1 + 0.5 * p(x_i) * h))
        c_val.append(1)
        d_val.append((h ** 2 * f(x_i)) / (1 + 0.5 * p(x_i) * h))
    a_val.append(- beta1 / h)
    b_val.append(beta0 + beta1 / h)
    c_val.append(None)
    d_val.append(B)

    y = thomas_algorithm(a_val, b_val, c_val, d_val)
    with open('boundary_value_problem.log', 'w') as log_file:
        log_file.write('Решение дифференциального уравнения, удовлетворющее краевым условиям\n')
        log_file.write(f'{y}')

