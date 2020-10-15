import numpy as np
import matplotlib.pyplot as plt


def beck_criteria(A_rad, A_mid):
    A_mid = np.abs(np.linalg.inv(A_mid))
    A = np.matmul(A_mid, A_rad)
    eig = np.abs(np.linalg.eigvals(A))
    spec_rad = max(eig)
    if spec_rad < 1:
        return True
    else:
        return False


def create_two_matrices(eps, n):
    N_m = np.ones(shape=(n, n)) * eps / 2
    N_r = np.ones(shape=(n, n)) * eps / 2
    for i in range(n):
        N_m[i][i] = 1
        N_r[i][i] = 0
    return N_r, N_m


if __name__ == '__main__':
    ######## First Task ########
    A = [[1, 1], [1.1, 1]]
    E = np.ones(shape=(2, 2))
    U, s, V = np.linalg.svd(A, full_matrices=True)
    U1, s_e, V = np.linalg.svd(E, full_matrices=True)
    a = min(s) / max(s_e)
    print('Eps <', a)

    ######## Second Task #######
    plt.figure()
    eps = 0
    step = 0.0001
    y = []
    for n in range(2, 100):
        eps = 0
        M_rad, M_mid = create_two_matrices(eps, n)
        while beck_criteria(M_rad, M_mid):
            eps += step
            M_rad, M_mid = create_two_matrices(eps, n)
        # print('n = ', n, 'eps = ', eps)
        y += [eps]
    plt.plot(range(2, 100), y)
    plt.title('Зависимость параметра eps от размерности матрицы')
    plt.xlabel('n')
    plt.ylabel('eps')
    plt.show()
