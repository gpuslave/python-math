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
            print(f"return m={i}")
            return YY

        # TODO: dunno how to fix atm
        if i >= 3 and eps_i > eps_i_prev:
            print(f"return m={i-1}")
            return prev_YY

        prev_YY = YY
        eps_i_prev = eps_i

    return prev_YY


def func_a(x):
    """Function y = (1/10)x^3 + x^2 + (1/2)x"""
    return (1 / 10) * x**3 + x**2 + (1 / 2) * x


def func_b(x):
    """Function y = (1/2)x^4 + 2x^3 + (1/2)x^2 + (1/5)x"""
    return (1 / 2) * x**4 + 2 * x**3 + (1 / 2) * x**2 + (1 / 5) * x


def func_to_np(arr, func):
    return np.asarray(list(map(func, arr)), dtype=np.float64)


if __name__ == "__main__":
    # x_logsp = np.sort(
    #     np.concatenate(
    #         (
    #             -np.logspace(np.log10(0.000001), np.log10(10.14), 70),
    #             np.logspace(np.log10(0.000001), np.log10(10.14), 70),
    #         )
    #     )
    # )
    x = np.unique(func_to_np(np.linspace(-4, 1.5, 10), func_b))

    # X = x_logsp
    X = np.sort(x)

    Y = func_to_np(X, func_a)

    N = len(X)
    m = 7

    XX = -12.5
    eps = 1e-12

    YY = lagrange_interpolation(X, Y, N, XX, m, eps)

    print(f"Интерполяционное значение в точке {XX}: {YY}")
    print(f"Значение: {func_a(XX)}, разность: {abs(YY-func_a(XX))}")
