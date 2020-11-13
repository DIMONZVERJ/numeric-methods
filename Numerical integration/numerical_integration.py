# coding=utf8
from typing import Optional

import pandas as pd
from common_functions import *

quadrature_coefficients = pd.DataFrame(index=[1, 2, 3, 4, 5],
                                       columns=['t', 'C'],
                                       data=[[[0], [2]],
                                             [[-1 / 3 ** 0.5, 1 / 3 ** 0.5], [1, 1]],
                                             [[-(3 / 5) ** 0.5, 0, (3 / 5) ** 0.5], [5 / 9, 8 / 9, 5 / 9]],
                                             [[-(3 / 7 - 2 / 7 * (6 / 5) ** 0.5) ** 0.5,
                                               -(3 / 7 + 2 / 7 * (6 / 5) ** 0.5) ** 0.5,
                                               +(3 / 7 - 2 / 7 * (6 / 5) ** 0.5) ** 0.5,
                                               +(3 / 7 + 2 / 7 * (6 / 5) ** 0.5) ** 0.5],
                                              [(18 + 30 ** 0.5) / 36,
                                               (18 - 30 ** 0.5) / 36,
                                               (18 + 30 ** 0.5) / 36,
                                               (18 - 30 ** 0.5) / 36]],
                                             [[0,
                                               -1 / 3 * (5 - 2 * (10 / 7) ** 0.5) ** 0.5,
                                               -1 / 3 * (5 + 2 * (10 / 7) ** 0.5) ** 0.5,
                                               +1 / 3 * (5 - 2 * (10 / 7) ** 0.5) ** 0.5,
                                               +1 / 3 * (5 + 2 * (10 / 7) ** 0.5) ** 0.5],
                                              [128 / 225,
                                               (322 + 13 * 70 ** 0.5) / 900,
                                               (322 - 13 * 70 ** 0.5) / 900,
                                               (322 + 13 * 70 ** 0.5) / 900,
                                               (322 - 13 * 70 ** 0.5) / 900]]
                                             ])


def write_log(function: callable(float), a, b, n, h: Optional, res):
    with open('numerical_integration.log', 'a') as file_log:
        file_log.write(f'{function.__doc__} a = {a}, b = {b}, n = {n}\n')
        if h is not None:
            file_log.write(f'h = {h}\n')
        file_log.write(f'Результат: {res}\n\n')


def log_error(error_p):
    with open('numerical_integration.log', 'a') as file_log:
        file_log.write(f'Погрешность: {error_p}\n')


def left_triangle(function: callable(float), a: float, b: float, n: int):
    """Метод левых прямоугольников"""
    assert n > 0 and a < b
    h = (b - a) / n
    res = sum([h * function(a + i * h)
               for i in range(n)])
    write_log(left_triangle, a, b, n, h, res)
    return res


def right_triangle(function: callable(float), a: float, b: float, n: int):
    """Метод правых прямоугольников"""
    assert n > 0 and a < b
    h = (b - a) / n
    res = sum([h * function(a + i * h)
               for i in range(1, n+1)])
    write_log(right_triangle, a, b, n, h, res)
    return res


def middle_triangle(function: callable(float), a: float, b: float, n: int):
    """Метод средних прямоугольников"""
    assert n > 0 and a < b
    h = (b - a) / n
    res = sum([h * function((a + i * h) + h / 2)
               for i in range(n)])
    write_log(middle_triangle, a, b, n, h, res)
    return res


def trapezoidal_rule(function: callable(float), a: float, b: float, n: int):
    """Формула трапеций"""
    assert n > 0 and a < b
    h = (b - a) / n
    res = h * ((function(a) + function(b))/2 + sum([function(a + i * h) for i in range(1, n)]))
    write_log(trapezoidal_rule, a, b, n, h, res)
    return res


def simpsons_rule(function: callable(float), a: float, b: float, n: int):
    """Формула Симпосона"""
    assert n > 0 and a < b and n % 2 == 0
    h = (b - a) / n
    res = h / 3 * (
            function(a) + function(b) +
            2 * sum([function(a + i * h) for i in range(2, n-1, 2)]) +
            4 * sum([function(a + i * h) for i in range(1, n, 2)])
    )
    write_log(simpsons_rule, a, b, n, h, res)
    return res


def gaussian_quadrature(function: callable(float), a: float, b: float, n: int):
    """Формула Гаусса"""
    assert n > 0 and a < b
    try:
        res = (b - a) / 2 * sum([C * function((b + a) / 2 + (b - a) / 2 * t)
                                 for C, t
                                 in zip(quadrature_coefficients.iloc[n-1]['C'], quadrature_coefficients.iloc[n-1]['t'])])
        write_log(gaussian_quadrature, a, b, n, None, res)
        return res
    except IndexError:
        with open('numerical_integration.log', 'a') as file_log:
            file_log.write(f'Выход за пределы индекса. Убедитесь, что для n = {n} присутствуют коэффиценты\n\n')
        return None


def calc_method(function: callable(float), method: callable(object), a: float, b: float):
    accuracy = 10 ** (-5)
    method_10 = method(function, a, b, 10)
    method_20 = method(function, a, b, 20)
    error = abs(method_20 - method_10)
    log_error(error)
    if error < accuracy:
        with open('numerical_integration.log', 'a') as file_log:
            file_log.write('Необходимая точность достигнута при количестве частей 10 и 20\n\n')
    else:
        method_40 = method(function, a, b, 40)
        error = abs(method_40 - method_20)
        log_error(error)
        if error < accuracy:
            with open('numerical_integration.log', 'a') as file_log:
                file_log.write('Необходимая точность достигнута при количестве частей 20 и 40\n\n')
        else:
            # Remainder order value
            m = None
            if method == trapezoidal_rule:
                m = 2
            elif method == simpsons_rule:
                m = 4
            # Runge rule
            value_runge = method_40 + method_20 ** 2 / (method_40 ** 2 - method_20 ** 2) * (method_40 - method_20)
            with open('numerical_integration.log', 'a') as file_log:
                file_log.write('Уточнённое значение по правилу Рунге '
                               f'(использовались значения полученные при количестве шагов 20 и 40): {value_runge}')


def gaussian_quadrature_run(function: callable(object), a: float, b: float):
    res_4 = gaussian_quadrature(function, a, b, 4)
    res_5 = gaussian_quadrature(function, a, b, 5)
    log_error(abs(res_5 - res_4))


if __name__ == "__main__":

    func = parse_function('function.txt')
    a_in = 2.6
    b_in = 3.4
    with open('numerical_integration.log', 'w'):
        pass
    calc_method(func, trapezoidal_rule, a_in, b_in)
    calc_method(func, simpsons_rule, a_in, b_in)
    gaussian_quadrature_run(func, a_in, b_in)
