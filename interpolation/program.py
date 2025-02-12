import numpy as np


def lagrange_interpolation(X, Y, N, XX, m):
    """
    Интерполирует функцию с помощью многочлена Лагранжа степени m на неравномерной сетке.

    Args:
        X (np.ndarray): Вектор значений аргументов (узлов интерполяции).
        Y (np.ndarray): Вектор значений функции в узлах интерполяции.
        N (int): Количество узлов интерполяции.
        XX (float): Значение аргумента, при котором вычисляется интерполяционное значение.
        m (int): Степень многочлена Лагранжа.

    Returns:
        float: Вычисленное интерполяционное значение функции в точке XX.
    """
    if N < m + 1:
        print("N < m + 1")
        exit(1)

    prev_YY = 0.0
    eps_i_prev = float("inf")

    for i in range(1, m + 1):
        distances = np.abs(X - XX)
        indices = np.argsort(distances)[:i]
        indices.sort()

        l_values = np.zeros(i)
        for j in range(i):
            l_j = 1.0
            for k in range(i):
                if k != j:
                    l_j *= (XX - X[indices[k]]) / (X[indices[j]] - X[indices[k]])
            l_values[j] = l_j

        YY = np.sum(Y[indices] * l_values)

        eps_i = np.abs(YY - prev_YY)
        print(eps_i, eps_i_prev, YY, prev_YY)

        if i >= 2 and eps_i > eps_i_prev:
            print(i)
            return prev_YY

        prev_YY = YY
        eps_i_prev = eps_i

    return prev_YY


if __name__ == "__main__":
    # x = np.linspace(np.e, np.e**4, 100)
    # f_x = np.log(x)

    f_x = np.logspace(np.log10(0.000001), np.log10(3.14), 50)
    print(f_x)

    X = f_x
    Y = np.sin(X)

    N = len(X)
    m = 10

    XX = 0.1
    YY = lagrange_interpolation(X, Y, N, XX, m)

    print(f"Интерполяционное значение в точке {XX}: {YY}")
