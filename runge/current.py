import numpy as np

def runge_kutta(f, x_current, y_current, h, direction=1):
  """
  returns:
  (rk2, rk3)
  """
  h = direction * h

  k1 = h * f(x_current, y_current)
  k2 = h * f(x_current + h / 2, y_current + k1 / 2)
  k3 = h * f(x_current + h, y_current - k1 + k2 * 2)

  rk2 = y_current + k2
  rk3 = y_current + (k1 + k2 * 4 + k3) / 6
  return rk2, rk3
   

# def runge_kutta_3rd_order(f, y_current, x_current, h):
#     """
#     Решает один шаг задачи Коши методом Рунге-Кутта третьего порядка.

#     Аргументы:
#     f : функция, представляющая dy/dx = f(x, y)
#     y_current : текущее значение y
#     x_current : текущее значение x
#     h : размер шага

#     Возвращает:
#     y_next : следующее значение y
#     """

#     k1 = h * f(x_current, y_current)
#     k2 = h * f(x_current + h / 2, y_current + k1 / 2)
#     k3 = h * f(x_current + h, y_current - k1 + 2 * k2)
#     y_next = y_current + (k1 + 4 * k2 + k3) / 6

#     return y_next

# def runge_kutta_2nd_order(f, y_current, x_current, h):
#     """
#     Решает один шаг задачи Коши методом Рунге-Кутта второго порядка.

#     Аргументы:
#     f : функция, представляющая dy/dx = f(x, y)
#     y_current : текущее значение y
#     x_current : текущее значение x
#     h : размер шага

#     Возвращает:
#     y_next : следующее значение y
#     """
#     k1 = h * f(x_current, y_current)
#     k2 = h * f(x_current + h / 2, y_current + k1 / 2)
#     y_next = y_current + k2
#     return y_next



def solve_cauchy_adaptive_step(f, y0, C, x_start, x_end, h_initial, h_min, tolerance):
    """
    Решает задачу Коши с автоматическим выбором шага, используя методы Рунге-Кутта 2-го и 3-го порядков.

    Аргументы:
    f : функция, представляющая dy/dx = f(x, y)
    y0 : начальное значение y в x_start
    x_start : начальное значение x
    x_end : конечное значение x
    h_initial : начальный размер шага
    tolerance : желаемая точность
    """

    x_current = C
    y_current = y0
    h = h_initial
    h_min_flag = False

    direction = 1 if C == x_start else -1
    if direction == -1:
      x_start, x_end = x_end, x_start 

    print(x_current, direction, x_start, x_end)

    with open("output.txt", "a", encoding="utf-8") as output_file:
        output_file.truncate(0)

        output_file.write(
            f"x:{x_current:.5f}\t\ty:{y_current:.5f}\t\th:{h:.5f}\t\te:{tolerance:.5f}\n"
        )

        while x_current * direction < x_end * direction:

            if direction * (x_current + h * direction) > x_end * direction:
                h = np.abs(x_end - x_current)  # Чтобы точно достичь x_end

            error_estimation = tolerance + 1  # Инициализация для входа в цикл while

            while error_estimation > tolerance and not h_min_flag:
                if h < h_min:
                    output_file.write(
                        f"Предупреждение: Размер шага стал очень маленьким ({h}). Точность может быть не достигнута."
                    )
                    h_min_flag = True
                    break  # Выход из внутреннего цикла, если шаг слишком мал

                y_rk2_step, y_rk3_step = runge_kutta(f, x_current, y_current, h, direction)

                # Оценка ошибки как разница между 2-м и 3-м порядками
                error_estimation = np.abs(y_rk2_step - y_rk3_step)

                if error_estimation > tolerance:
                    h *= 0.5  # Уменьшение шага вдвое
                else:
                    break  # Шаг приемлем

            if error_estimation <= tolerance:
                y_current = y_rk2_step  # Принимаем решение 2-го порядка
                x_current += direction * h

                output_file.write(
                    f"x:{x_current:.5f}\t\ty:{y_current:.5f}\t\th:{h:.5f}\t\te:{tolerance:.8f}\t\test:{error_estimation:.8f}\t\tdiff:{(error_estimation - tolerance):.8f}\n"
                )

                if h * 2 <= np.abs(x_end - x_current) and h * 2 <= h_initial * 2:
                    h *= 2.0  # Увеличение шага вдвое для следующего шага, если это возможно и не превышает начальный шаг удвоенный.

            else:
                if h <= h_min and x_current * direction < x_end * direction:
                    y_current, _ = runge_kutta(f, x_current, y_current, h)
                    x_current += direction * h

                    output_file.write(
                        f"x:{x_current:.5f}\t\ty:{y_current:.5f}\t\th:{h:.5f}\t\te:{tolerance:.8f}\t\test:{error_estimation:.8f}\t\tdiff:{(error_estimation - tolerance):.8f}\n"
                    )

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


def parse_output_file(filename="output.txt"):
    """Parse the output file from Runge-Kutta solver"""
    x_vals, y_vals, h_vals, err_vals = [], [], [], []

    with open(filename, "r", encoding="utf-8") as file:
        for line in file:
            if line.startswith("x:"):
                parts = line.strip().split("\t\t")
                x_vals.append(float(parts[0].split(":")[1]))
                y_vals.append(float(parts[1].split(":")[1]))
                h_vals.append(float(parts[2].split(":")[1]))

                # Extract error estimation if available
                if len(parts) > 4 and "est:" in parts[4]:
                    err_vals.append(float(parts[4].split(":")[1]))
                else:
                    err_vals.append(0)

    return x_vals, y_vals, h_vals, err_vals


def visualize_results():
    """Plot the numerical solution and step size from output.txt"""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    # Parse data
    x_vals, y_vals, h_vals, err_vals = parse_output_file()

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Solution plot
    ax1.plot(x_vals, y_vals, "b-", linewidth=2)
    ax1.set_title("Numerical Solution y(x)")
    ax1.set_ylabel("y")
    ax1.grid(True)

    # Step size plot (log scale for better visualization)
    ax2.plot(x_vals[1:], h_vals[1:], "r-", linewidth=1.5)
    # ax2.set_ylim((min(h_vals) * 0.001, max(h_vals) * 0.3))
    # y_ticks = np.linspace(min(h_vals) * 0.9, max(h_vals) * 1.1, 10)
    # y_ticks = np.logspace(np.log10(min(h_vals) * 0.9), np.log10(max(h_vals) * 1.1), 10)
    # ax2.set_yticks(y_ticks)

    ax2.set_title("Adaptive Step Size")
    ax2.set_xlabel("x")
    ax2.set_ylabel("Step Size (h)")
    ax2.grid(True)

    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax2.yaxis.set_major_formatter(formatter)

    plt.tight_layout()
    plt.savefig("runge_kutta_solution.png", dpi=300)
    plt.close()


if __name__ == "__main__":

    # Пример использования: dy/dx = y - x^2 + 1, y(0) = 0.5
    def f_example(x, y):
        return y - x**2 + 1
    input_filename = "input.txt"

    # def f_example(x, y):
    #     return 2 * x
    # input_filename = "input2.txt"

    x_start_example, x_end_example, C, y0_example, h_min, tolerance_example = (
        read_vars_from_file(input_filename)
    )

    h_initial_example = (x_end_example - x_start_example) / 10.0

    print(
        x_start_example,
        x_end_example,
        C,
        y0_example,
        h_min,
        tolerance_example,
        h_initial_example,
    )

    solve_cauchy_adaptive_step(
        f_example,
        y0_example,
        C,
        x_start_example,
        x_end_example,
        h_initial_example,
        h_min,
        tolerance_example,
    )

    visualize_results()
