import numpy as np
from numpy import unravel_index


def solve_eq(arr, file_obj):
    # положение максимального по модулю элемента в массиве (строка, столбец)
    row, col = unravel_index(abs(arr[:, :-1]).argmax(), shape=arr[:, :-1].shape)
    file_obj.write('Исходная матрица:\n')
    np.savetxt(file_obj, arr, delimiter=' | ', fmt='%10.5f')
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
    np.savetxt(file_obj, arr, delimiter=' | ', fmt='%10.5f')
    # удаляем главную строку и главный столбец
    new_arr = np.delete(np.delete(arr, row, 0), col, 1)
    file_obj.write('Матрица без главной строки и главного столбца:\n')
    np.savetxt(file_obj, new_arr, delimiter=' | ', fmt='%10.5f')
    file_obj.write('\n\n\n')
    # находим решение для новой матрицы
    res = solve_eq(new_arr, file_obj)
    file_obj.write(f'Решение системы из {new_arr.shape[0]} уравнений: {res}\n\n')
    # находим значение x в главной строке
    res.insert(col, (arr[row, -1] - (np.delete(arr[row, :-1], col, 0) * res).sum()) / max_el)
    return res


if __name__ == "__main__":
    arr = np.loadtxt('file1.txt')
    with open('gaussian_method_log.txt', 'w') as log_file:
        ans = solve_eq(arr, log_file)
        log_file.write('Ответ: ')
        np.savetxt(log_file, ans, fmt='%10.5f', newline=' | ')
