# encoding: utf-8
import numpy as np
from numpy import unravel_index


def parse_function(file_name):
    with open(file_name, 'r') as file:
        return eval(f'lambda x: {file.readline()}')


def thomas_algorithm(a_val: list, b_val: list, c_val: list, d_val: list):
    """Метод прогонки"""
    p_val = []
    q_val = []

    # init p_1 and q_1
    p_val.append(-c_val[0] / b_val[0])
    q_val.append(d_val[0] / b_val[0])

    x_val = []
    for i in range(1, len(a_val)):
        # find p_i and q_i for i from 2 to n - 1.
        # NOTE: for q_i we find also q_n to further find x_n
        if i < len(a_val) - 1:
            p_val.append(-c_val[i] / (a_val[i] * p_val[i - 1] + b_val[i]))
        q_val.append((d_val[i] - a_val[i] * q_val[i - 1]) / (a_val[i] * p_val[i - 1] + b_val[i]))

    # find x_n
    x_val.append(q_val.pop())
    # execute a reverse step to find x_(n-1)..x_1
    for i in range(len(a_val) - 1, 1, -1):
        x_val.append(p_val.pop() * x_val[-1] + q_val.pop())
    # find x_0
    x_val.append(1 / b_val[0] * (-c_val[0] * x_val[1] + d_val[0]))
    # reverse list to normal order: x_0, x_1..x_n
    x_val.reverse()
    return x_val


def solve_eq(arr, file_obj, delimiter, fmt):
    # положение максимального по модулю элемента в массиве (строка, столбец)
    row, col = unravel_index(abs(arr[:, :-1]).argmax(), shape=arr[:, :-1].shape)
    file_obj.write('Исходная матрица:\n')
    np.savetxt(file_obj, arr, delimiter=delimiter, fmt=fmt)
    file_obj.write(f'Максимальный элемент находится в {row + 1} строке и {col + 1} столбце\n')
    # максимальный элемент
    max_el = (arr[:, :-1])[row, col]
    file_obj.write(f'Максимальный по модулю элемент: {max_el}\n')
    # если система уравнений сведена к виду m*x = n, то находим x = n/m
    if arr.shape[0] == 1:
        return [arr[0, 1] / arr[0, 0]]
    # для каждого элемента в столбце (за исключением главного элемента)
    for i in range(arr.shape[0]):
        if i != row:
            # находим m_i
            m = -arr[i, col] / max_el
            file_obj.write(f'Множитель для {i + 1}-ой строки равен: {m}\n')
            # умножаем m_i на главную строку и суммируем i-ой строкой
            arr[i] = arr[row] * m + arr[i]
    file_obj.write('Матрица после умножения:\n')
    np.savetxt(file_obj, arr, delimiter=delimiter, fmt=fmt)
    # удаляем главную строку и главный столбец
    new_arr = np.delete(np.delete(arr, row, 0), col, 1)
    file_obj.write('Матрица без главной строки и главного столбца:\n')
    np.savetxt(file_obj, new_arr, delimiter=delimiter, fmt=fmt)
    file_obj.write('\n\n\n')
    # находим решение для новой матрицы
    res = solve_eq(new_arr, file_obj, delimiter, fmt)
    file_obj.write(f'Решение системы из {new_arr.shape[0]} уравнений: {res}\n\n')
    # находим значение x в главной строке
    res.insert(col, (arr[row, -1] - (np.delete(arr[row, :-1], col, 0) * res).sum()) / max_el)
    return res
