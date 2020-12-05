import matplotlib.pyplot as plt
from common_functions import thomas_algorithm

file_log = open('result.txt', 'w')


def h_value(x_next, x_now):
    return x_next - x_now


def l_value(x_prev, x_now, x_next):
    return h_value(x_next, x_now) / (h_value(x_next, x_now) + h_value(x_now, x_prev))


def get_m_values(x_values):
    file_log.write('Находим значения производной m_i в точках x_i\n\n\n')
    # a_now = 1 - l_now; a_0 = None
    a = [None] + [1 - l_value(x[i - 1], x[i], x[i + 1]) for i in range(len(x_values) - 1)]
    b = [2] * len(x_values)

    # c_now = l_now; c_n = None
    c = [l_value(x[i - 1], x[i], x[i + 1]) for i in range(len(x_values) - 1)] + [None]

    # d_now = 6 * y_prev / (h_prev * (h_now + h_prev)) - 6 * y_now / (h_now * h_prev) + 6 * y_next / (h_now * (h_now + h_prev))
    d = [0] + \
        [6 * (y[i - 1] / (h_value(x[i], x[i - 1]) * (h_value(x[i + 1], x[i]) + h_value(x[i], x[i - 1]))) -
              y[i] / (h_value(x[i + 1], x[i]) * h_value(x[i], x[i - 1])) +
              y[i + 1] / (h_value(x[i + 1], x[i]) * (h_value(x[i + 1], x[i]) + h_value(x[i], x[i - 1]))))
         for i in range(1, len(x_values) - 1)] + [0]
    m_values = thomas_algorithm(a, b, c, d)
    file_log.write(f'Коэффициенты a, b, c, d для метода прогонки \na = {a}\nb = {b}\nc = {c}\nd = {d}\n\n\n')
    file_log.write(f'Значения коэффициентов m_i = {m_values}\n\n\n')
    return m_values


def spline_interpolation(x_values, y_values):
    assert len(x_values) == len(y_values)
    m_values = get_m_values(x_values)
    a_values = [0] * (len(x_values) - 1)
    b_values = [0] * (len(x_values) - 1)
    c_values = [0] * (len(x_values) - 1)
    d_values = [0] * (len(x_values) - 1)
    for i in range(len(x_values) - 1):
        a_values[i] = y[i]
        c_values[i] = m_values[i] / 2
        d_values[i] = (m_values[i + 1] - m_values[i]) / (6 * h_value(x[i + 1], x[i]))
        b_values[i] = (y[i + 1] - y[i]) / h_value(x[i + 1], x[i]) - 1 / 3 * h_value(x[i + 1], x[i]) * m_values[
            i] - 1 / 6 * h_value(x[i + 1], x[i]) * m_values[i + 1]

    functions_str = [f'S_i(x) = {a_values[i]} + {b_values[i]} * (x - {x_values[i]}) + {c_values[i]} * (x - {x_values[i]}) ** 2 + {d_values[i]} * (x - {x_values[i]}) ** 3'
                     for i in range(len(a_values))]
    for func in functions_str:
        file_log.write(f'{func}\n\n\n')

    return [eval(f'lambda x: {a_values[i]} + {b_values[i]} * (x - {x_values[i]}) + {c_values[i]} * (x - {x_values[i]}) ** 2 + {d_values[i]} * (x - {x_values[i]}) ** 3')
            for i in range(len(a_values))]


def plot_spline(x, y, spline_functions,count_step):
    new_x = []
    new_y = []
    for i in range(len(x) - 1):
        step_size = (x[i + 1] - x[i]) / count_step
        for j in range(count_step):
            new_x.append(x[i] + j * step_size)
            new_y.append(spline_functions[i](x[i] + j * step_size))
    new_x.append(x[-1])
    new_y.append(y[-1])

    plt.plot(new_x, new_y)
    plt.plot(x, y, marker='o', color='red')
    plt.show()


if __name__ == "__main__":
    x = [1.415,
         1.420,
         1.425,
         1.430,
         1.435,
         1.440]

    y = [0.888551,
         0.889599,
         0.890637,
         0.891667,
         0.892687,
         0.893698]
    x_value = 1.4258
    x = [1, 3, 4, 5, 6, 8]
    y = [1, 3, 5, 2, 0, 4]
    spline_functions = spline_interpolation(x, y)
    file_log.write(f'X = {x_value}, Y = {spline_functions[list(map(lambda x_: 1 if x_ > x_value else 0, x)).index(1) - 1](x_value)}')
    file_log.close()
    plot_spline(x, y, spline_functions, 5)
