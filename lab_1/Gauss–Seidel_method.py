import numpy as np


def seidel(c, b, eps):
    if max(np.sum(c, axis=1)) < 1:
        print('Метод Зейделя сходится')
        # начальное приближение
        X_n = np.zeros(c.shape[0])
        j = 1
        while True:
            # копия прошлых значений для проверки точности
            temp = np.copy(X_n)
            # получение нового приближения
            for i in range(c.shape[0]):
                X_n[i] = c[i] @ X_n + b[i]
            acc = abs(max(temp - X_n))
            print(f'Итерация {j}: {X_n} точность: {acc}')
            if acc < eps:
                print('Необходимая точность достигнута')
                return X_n
            j += 1
    else:
        print('Система не удовлетворяет критерию сходимости')
        return None


if __name__ == "__main__":
    arr = np.loadtxt('file2.txt')
    seidel(arr[:, :-1], arr[:, -1], 0.00001)
