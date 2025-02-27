import numpy as np


def runge_kutta_2nd_order(f, y_current, x_current, h):
    """
    Решает один шаг задачи Коши методом Рунге-Кутта второго порядка.

    Аргументы:
    f : функция, представляющая dy/dx = f(x, y)
    y_current : текущее значение y
    x_current : текущее значение x
    h : размер шага

    Возвращает:
    y_next : следующее значение y
    """
    k1 = h * f(x_current, y_current)
    k2 = h * f(x_current + h / 2, y_current + k1 / 2)
    y_next = y_current + k2
    return y_next


def runge_kutta_3rd_order(f, y_current, x_current, h):
    """
    Решает один шаг задачи Коши методом Рунге-Кутта третьего порядка.

    Аргументы:
    f : функция, представляющая dy/dx = f(x, y)
    y_current : текущее значение y
    x_current : текущее значение x
    h : размер шага

    Возвращает:
    y_next : следующее значение y
    """

    k1 = h * f(x_current, y_current)
    k2 = h * f(x_current + h / 2, y_current + k1 / 2)
    k3 = h * f(x_current + h, y_current - k1 + 2 * k2)
    y_next = y_current + (k1 + 4 * k2 + k3) / 6

    return y_next


def runge_kutta_4th_order(f, y_current, x_current, h):
    """
    Решает один шаг задачи Коши методом Рунге-Кутта четвертого порядка.

    Аргументы:
    f : функция, представляющая dy/dx = f(x, y)
    y_current : текущее значение y
    x_current : текущее значение x
    h : размер шага

    Возвращает:
    y_next : следующее значение y
    """
    k1 = h * f(x_current, y_current)
    k2 = h * f(x_current + h / 2, y_current + k1 / 2)
    k3 = h * f(x_current + h / 2, y_current + k2 / 2)
    k4 = h * f(x_current + h, y_current + k3)
    y_next = y_current + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return y_next


def solve_cauchy_adaptive_step(f, y0, x_start, x_end, h_initial, h_min, tolerance):
    """
    Решает задачу Коши с автоматическим выбором шага, используя методы Рунге-Кутта 2-го и 4-го порядков.

    Аргументы:
    f : функция, представляющая dy/dx = f(x, y)
    y0 : начальное значение y в x_start
    x_start : начальное значение x
    x_end : конечное значение x
    h_initial : начальный размер шага
    tolerance : желаемая точность

    Возвращает:
    x_values : список значений x
    y_values : список соответствующих значений y
    """
    x_current = x_start
    y_current = y0
    h = h_initial
    h_min_flag = False

    x_values = [x_start]
    y_values = [y0]

    with open("output.txt", "a") as output_file:
        output_file.truncate(0)

        print(
            f"x:{x_current:.5f}\t\ty:{y_current:.5f}\t\th:{h:.5f}\t\te:{tolerance:.5f}\n"
        )
        output_file.write(
            f"x:{x_current:.5f}\t\ty:{y_current:.5f}\t\th:{h:.5f}\t\te:{tolerance:.5f}\n"
        )

        while x_current < x_end:
            if x_current + h > x_end:
                h = x_end - x_current  # Чтобы точно достичь x_end

            error_estimation = tolerance + 1  # Инициализация для входа в цикл while

            while error_estimation > tolerance and not h_min_flag:
                if h < h_min:
                    print(
                        f"Предупреждение: Размер шага стал очень маленьким ({h}). Точность может быть не достигнута."
                    )
                    h_min_flag = True
                    break  # Выход из внутреннего цикла, если шаг слишком мал

                # Вычисление решения с шагом h, используя 2-й порядок
                y_rk2_step = runge_kutta_2nd_order(f, y_current, x_current, h)

                # Вычисление решения с шагом h, используя 3-й порядок (для оценки ошибки)
                y_rk3_step = runge_kutta_3rd_order(f, y_current, x_current, h)

                # Оценка ошибки как разница между 2-м и 4-м порядками
                error_estimation = np.abs(y_rk2_step - y_rk3_step)

                if error_estimation > tolerance:
                    h *= 0.5  # Уменьшение шага вдвое
                else:
                    break  # Шаг приемлем

            if error_estimation <= tolerance:
                y_current = y_rk2_step  # Принимаем решение 2-го порядка
                x_current += h

                print(
                    f"x:{x_current:.5f}\t\ty:{y_current:.5f}\t\th:{h:.5f}\t\te:{tolerance:.5f}\n"
                )
                output_file.write(
                    f"x:{x_current:.5f}\t\ty:{y_current:.5f}\t\th:{h:.5f}\t\te:{tolerance:.5f}\n"
                )
                # x_values.append(x_current)
                # y_values.append(y_current)

                if h * 2 <= (x_end - x_current) and h * 2 <= h_initial * 2:
                    h *= 2.0  # Увеличение шага вдвое для следующего шага, если это возможно и не превышает начальный шаг удвоенный.

            else:
                if h <= h_min and x_current < x_end:
                    y_current = runge_kutta_2nd_order(f, y_current, x_current, h)
                    x_current += h
                    print(
                        f"x:{x_current:.5f}\t\ty:{y_current:.5f}\t\th:{h:.5f}\t\te:{tolerance:.5f}\n"
                    )
                    output_file.write(
                        f"x:{x_current:.5f}\t\ty:{y_current:.5f}\t\th:{h:.5f}\t\te:{tolerance:.5f}\n"
                    )
                    # x_values.append(x_current)
                    # y_values.append(y_current)

    return x_values, y_values


def read_vars_from_file(filename):
    file = open(filename, mode="r")
    readline = file.readlines()

    x_start = float(readline[0].strip())
    x_end = float(readline[1].strip())
    C = float(readline[2].strip())
    y0 = float(readline[3].strip())
    h_min = tolerance = float(readline[4].strip())
    tolerance = float(readline[5].strip())

    return x_start, x_end, C, y0, h_min, tolerance


if __name__ == "__main__":
    input_filename = "input.txt"

    # Пример использования: dy/dx = y - x^2 + 1, y(0) = 0.5
    def f_example(x, y):
        return y - x**2 + 1

    # y0_example = 0.5
    # x_start_example = 0
    # x_end_example = 6
    # tolerance_example = 1e-5

    x_start_example, x_end_example, C, y0_example, h_min, tolerance_example = (
        read_vars_from_file(input_filename)
    )

    h_initial_example = (x_end_example - x_start_example) / 10.0
    print(
        x_start_example,
        x_end_example,
        y0_example,
        h_min,
        tolerance_example,
        h_initial_example,
    )

    x_vals, y_vals = solve_cauchy_adaptive_step(
        f_example,
        y0_example,
        x_start_example,
        x_end_example,
        h_initial_example,
        h_min,
        tolerance_example,
    )

    # with open("output.txt", "w") as output_file:
    #     print("x\t\ty")
    #     for x, y in zip(x_vals, y_vals):
    #         output_file.write(f"{x:.3f}\t\t{y:.6f}\n")

    import matplotlib.pyplot as plt

    plt.plot(
        x_vals, y_vals, marker="o", linestyle="-", label="Решение с адаптивным шагом"
    )

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Решение задачи Коши методом Рунге-Кутта с адаптивным шагом")
    plt.legend()
    plt.grid(True)
    plt.show()
    print(
        "[Image of Plot of solution of Cauchy problem using Runge-Kutta method with adaptive step]"
    )
