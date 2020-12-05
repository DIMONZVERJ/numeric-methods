# encoding: utf-8

from matplotlib import pyplot as plt
import numpy as np
import math


log_file = open('differential_equations.log', 'w')


def euler_method(function, x0, xn, h, y0):
    res = dict(x=[x0], y=[y0])
    for i in np.arange(x0+h, xn+h, h):
        x = i
        y = res['y'][-1] + h * function(res['x'][-1], res['y'][-1])
        res['x'].append(x)
        res['y'].append(y)
    return res


def euler_cauchy_method(function, x0, xn, h, y0):
    log_file.write('Info. Метод Эйлера-Коши\n')
    log_file.write(f"y' = x^2 - y, f({x0}) = {y0}, xn = {xn}, h = {h}\n")
    res = dict(x=[x0], y=[y0])
    j = 0
    for i in np.arange(x0+h, xn+h, h):
        x = i
        y_temp = res['y'][-1] + h * function(res['x'][-1], res['y'][-1])
        y = res['y'][-1] + h / 2 * (func(res['x'][-1], res['y'][-1]) + func(res['x'][-1], y_temp))
        log_file.write(f'x_{j+1} = {i}, y_temp = {y_temp}, y_{j+1} = {y} \n')
        j += 1
        res['x'].append(x)
        res['y'].append(y)
    log_file.write(f'Вычисление закончено, общее количество шагов: {j}\n\n')
    return res


def rynge_rule(y1, y2, k):
    return y2 + (y2 - y1) / (2 ** k - 1)


def func(x, y):
    return x ** 2 - y


def common_solution(x, c):
    return c * math.exp(-x) + x ** 2 - 2 * x + 2


def array_values(function, x0, xn, h, c):
    res = dict(x=[], y=[])
    for i in np.arange(x0, xn+h, h):
        res['x'].append(i)
        res['y'].append(function(i, c))
    return res


if __name__ == "__main__":
    kwargs = dict(x0=0, xn=1, h=0.1)
    res_euler_cauchy_1 = euler_cauchy_method(func, y0=2, **kwargs)
    res_solve = array_values(common_solution, c=0, **kwargs)
    plt.title(f'Решение задачи Коши: y\'=x^2 - y, x0 = {kwargs["x0"]}, xn = {kwargs["xn"]}, y0 = 2, h = {kwargs["h"]}')
    plt.plot(res_euler_cauchy_1['x'], res_euler_cauchy_1['y'], color='blue')
    plt.plot(res_solve['x'], res_solve['y'], color='green')
    plt.legend(["Метод Эйлера-Коши", "Аналитическое решение"])
    plt.show()

    kwargs = dict(x0=0, xn=1, h=0.1 / 2)
    res_euler_cauchy_2 = euler_cauchy_method(func, y0=2, **kwargs)
    res_solve = array_values(common_solution, c=0, **kwargs)
    plt.title(f'Решение задачи Коши: y\'=x^2 - y, x0 = {kwargs["x0"]}, xn = {kwargs["xn"]}, y0 = 2, h = {kwargs["h"]}')
    plt.plot(res_euler_cauchy_2['x'], res_euler_cauchy_2['y'], color='blue')
    plt.plot(res_solve['x'], res_solve['y'], color='green')
    plt.legend(["Метод Эйлера-Коши", "Аналитическое решение"])
    plt.show()

    log_file.write('Info. Оценка погрешности при помощи правила Рунге\n')
    y1, y2 = res_euler_cauchy_1['y'][-1], res_euler_cauchy_2['y'][-1]
    val = rynge_rule(y1, y2, 2)
    log_file.write(f'y_h = {y1}, y_h/2 = {y2}, res = {val}')
    log_file.close()