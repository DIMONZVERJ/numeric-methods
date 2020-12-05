# encoding: utf-8
import math
import numpy as np
from common_functions import solve_eq

log_file = open('fredholm_volterra_eq.log', 'w')
save_txt_attr = dict(delimiter=' | ', fmt='%12.5e')


def k(x_1, x_2):
    return x_2 * math.sin(2 * math.pi * x_1)


def f(x):
    return x


def fredholm_eq(lambda_var, a, b, n):
    h = (b - a) / n
    coefficients = []
    for i in range(n + 1):
        coefficients.append(list(map(lambda x: - x * lambda_var,
                                     [h / 2 * int((j != 0 and j != n) + 1) * k(a + i * h, a + j * h)
                                      for j in range(n + 1)])))
        coefficients[i][i] += 1
        coefficients[i].append(f(i))
    ans = np.array(solve_eq(np.array(coefficients), log_file, **save_txt_attr))
    log_file.write('Ответ: ')
    np.savetxt(log_file, [ans], **save_txt_attr)
    log_file.write('\n\n\n')
    return ans


def volterra_eq(lambda_var, a, b, n):
    y = []
    h = (b - a) / n
    for i in range(n + 1):
        x_i = a + i * h
        x_0 = a
        if i == 0:
            y.append(f(x_i))
        else:
            y.append((f(x_i) + lambda_var * h * (k(x_i, x_0) + 2 * sum([k(x_i, a + j * h) * y[j] for j in range(i)]) / 2)) /
                     (1 - lambda_var * h * k(x_i, x_i) / 2))
    log_file.write('\n')
    log_file.write('Ответ: ')
    np.savetxt(log_file, [y], **save_txt_attr)
    return y


if __name__ == "__main__":
    lambda_var = 1
    a = 0
    b = 1
    n = 10
    log_file.write('Уравнение Фредгольма 2-ого рода\n')
    fredholm_eq(lambda_var, a, b, n)
    log_file.write('Уравнение Вольтерра 2-ого рода\n')
    volterra_eq(lambda_var, a, b, n)
    log_file.close()