import numpy as np

np.set_printoptions(suppress=True)


def lagrange_interpolation(X, Y, N, XX, m, eps):
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

    print(Y)
    prev_YY = 0.0
    eps_i_prev = float("inf")

    for i in range(1, m + 1):
        if i == 2:
            prev_YY = 0.0
            eps_i_prev = float("inf")

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
        print(f"m={i}")
        print(f"Epsilon: {eps_i}")
        print(f"Epsilon prev: {eps_i_prev}")
        print(f"Diff: {abs(eps_i - eps_i_prev)}")
        print(f"YY: {YY}")
        print(f"YY_prev: {prev_YY}")

        if eps_i < eps:
            return YY

        # TODO: dunno how to fix atm
        if i >= 3 and eps_i > eps_i_prev:
            print(i)
            return prev_YY

        prev_YY = YY
        eps_i_prev = eps_i

    return prev_YY


def ff(x):
    return (1 / 10) * x**3 + x**2 + (1 / 2) * x


def my_func(arr):
    return np.asarray(list(map(ff, arr)), dtype=np.float64)


if __name__ == "__main__":
    # x = np.linspace(np.e, np.e**4, 100)
    # f_x = np.log(x)

    f_x = np.sort(
        np.concatenate(
            (
                -np.logspace(np.log10(0.000001), np.log10(10.14), 70),
                np.logspace(np.log10(0.000001), np.log10(10.14), 70),
            )
        )
    )
    print(f_x)

    # print(my_func(f_x))

    X = f_x
    Y = my_func(f_x)

    N = len(X)
    m = 7

    XX = -0.26
    eps = 1e-10

    YY = lagrange_interpolation(X, Y, N, XX, m, eps)

    print(f"Интерполяционное значение в точке {XX}: {YY}")
