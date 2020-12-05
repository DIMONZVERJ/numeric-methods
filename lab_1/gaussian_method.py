import numpy as np
from common_functions import solve_eq


if __name__ == "__main__":
    arr = np.loadtxt('file1.txt')
    with open('gaussian_method_log.txt', 'w') as log_file:
        ans = solve_eq(arr, log_file, delimiter=' | ', fmt='%12.5e')
        log_file.write('Ответ: ')
        np.savetxt(log_file, ans, fmt='%10.5f', newline=' | ')
