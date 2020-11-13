from scipy.misc import derivative
import time
import pandas as pd

from common_functions import *

pd.set_option('display.precision', 16)
pd.set_option('expand_frame_repr', False)
EPS = 10 ** (-5)
A = 0
B = 1


def write_log(function_name: str, res_df: pd.DataFrame):
    with open(f'{function_name}.log', 'w') as f:
        f.write(res_df.__repr__())
        f.write('\n\n\n\n')
        f.write(f'count iterations: {res_df.index.stop}')


def bisection_method(function: callable(object), a, b, eps):
    assert a < b
    x = (a + b) / 2
    res_df = pd.DataFrame(columns=['a', 'c', 'b', 'f(a)', 'f(c)', 'f(b)', 'time executing'])
    while abs(function(x)) >= eps:
        start = time.time()
        res_df = res_df.append({'a': a, 'c': x, 'b': b, 'f(a)': function(a), 'f(c)': function(x), 'f(b)': function(b)},
                               ignore_index=True)
        if function(a) == 0:
            x = a
        if function(b) == 0:
            x = b

        if function(a) * function(x) < 0:
            b = x
        elif function(x) * function(b) < 0:
            a = x
        x = (a + b) / 2
        res_df.loc[:, 'time executing'].iloc[-1] = time.time() - start

    write_log(bisection_method.__name__, res_df)
    return res_df


def secant_method(function: callable(object), a, b, eps):
    assert a < b
    x_n = b
    res_df = pd.DataFrame(columns=['x_i', 'f(x_i)'])
    if function(a) == 0:
        x_n = a
    while abs(function(x_n)) >= eps:
        x_n = a - function(a) / (function(x_n) - function(a)) * (x_n - a)
        res_df = res_df.append({'x_i': x_n, 'f(x_i)': function(x_n)}, ignore_index=True)

    write_log(secant_method.__name__, res_df)
    return res_df


def newtons_method(function: callable(object), a, b, eps):
    assert a < b
    x_n = a
    if function(b) == 0:
        x_n = b

    res_df = pd.DataFrame(columns=['x_i', 'f(x_i)'])
    while abs(function(x_n)) >= eps:
        x_n = x_n - function(x_n) / derivative(function, x_n, dx=1)
        res_df = res_df.append({'x_i': x_n, 'f(x_i)': function(x_n)}, ignore_index=True)
    write_log(newtons_method.__name__, res_df)
    return res_df


def combine_secant_tangent_method(function: callable(object), a, b, eps):
    assert a < b
    x_n = a
    x_n_line = b

    if function(a) == 0:
        x_n_line = a
    if function(b) == 0:
        x_n = b

    res_df = pd.DataFrame(columns=['x_n', 'x_n_line', 'f(x_n)', 'f(x_n_line)', 'x_n_line - x_n'])
    while abs(x_n - x_n_line) >= eps:
        x_n, x_n_line = x_n - function(x_n) / derivative(function, x_n, dx=1), \
                        x_n - function(x_n) / (function(x_n_line) - function(x_n)) * (x_n_line - x_n)
        res_df = res_df.append({'x_n': x_n,
                                'x_n_line': x_n_line,
                                'f(x_n)': function(x_n),
                                'f(x_n_line)': function(x_n_line),
                                'x_n_line - x_n': x_n_line - x_n
                                },
                               ignore_index=True)
    write_log(combine_secant_tangent_method.__name__, res_df)
    return res_df


def iteration_method(function: callable(object), a, b, eps):
    k = 6
    x = 0
    if function(a) == 0:
        return a
    if function(b) == 0:
        return b
    res_df = pd.DataFrame(columns=['x', 'f(x)'])
    while abs(function(x)) >= eps:
        x = x - function(x) / k
        res_df = res_df.append({'x': x, 'f(x)': function(x)}, ignore_index=True)
    write_log(iteration_method.__name__, res_df)
    return res_df


if __name__ == "__main__":
    f = parse_function('function.txt')
    for method in [bisection_method, secant_method, newtons_method, combine_secant_tangent_method, iteration_method]:
        method(f, A, B, EPS)
        print(f'{method.__name__} is executed!')
