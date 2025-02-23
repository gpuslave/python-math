import numpy as np


def runge_kutta_2nd_order(f, y_current, x_current, h):
    """
    Решает один шаг задачи Коши методом Рунге-Кутта второго порядка.
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
    """

    k1 = h * f(x_current, y_current)
    k2 = h * f(x_current + h / 2, y_current + k1 / 2)
    k3 = h * f(x_current + h / 2, y_current + k2 / 2)
    k4 = h * f(x_current + h, y_current + k3)
    y_next = y_current + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return y_next


def solve_cauchy_adaptive_step(
    f, y0, x_start, x_end, h_initial, tolerance, output_filename="output.txt"
):
    """
    Решает задачу Коши с автоматическим выбором шага, используя методы Рунге-Кутта 2-го и 4-го порядков.
    Записывает каждое вычисленное значение решения сразу в файл, не сохраняя все значения в памяти.
    Возвращает только статистику.
    """

    x_current = x_start
    y_current = y0
    h = h_initial
    h_min = h_initial / 1000.0  # Минимальный размер шага

    integration_points_count = 0
    accuracy_failed_points_count = 0
    min_step_count = 0

    try:
        with open(
            output_filename, "w"
        ) as output_file:  # Открываем файл для записи в начале функции
            output_file.write(
                f"X-координата Приближенное_значение_Y Локальная_погрешность\n"
            )  # Заголовок в файл

            output_file.write(
                f"{x_start:.6f} {y0:.6e} 0.000000e+00\n"
            )  # Записываем начальную точку

            while x_current < x_end:
                if x_current + h > x_end:
                    h = x_end - x_current  # Чтобы точно достичь x_end

                error_estimation = tolerance + 1  # Инициализация для входа в цикл while

                while error_estimation > tolerance:
                    if h < h_min:
                        print(
                            f"Предупреждение: Размер шага стал очень маленьким ({h}). Точность может быть не достигнута в точке x={x_current:.3f}."
                        )
                        accuracy_failed_points_count += 1
                        min_step_count += 1  # Счет минимальных шагов даже если точность не достигнута.
                        break  # Выход из внутреннего цикла, если шаг слишком мал

                    # Вычисление решения с шагом h, используя 2-й порядок
                    y_rk2_step = runge_kutta_2nd_order(f, y_current, x_current, h)
                    # Вычисление решения с шагом h, используя 4-й порядок (для оценки ошибки)
                    y_rk4_step = runge_kutta_4th_order(f, y_current, x_current, h)
                    # Оценка ошибки как разница между 2-м и 4-м порядками
                    error_estimation = np.abs(y_rk2_step - y_rk4_step)

                    if error_estimation > tolerance:
                        h *= 0.5  # Уменьшение шага вдвое
                    else:
                        break  # Шаг приемлем

                if error_estimation <= tolerance:
                    y_current = y_rk2_step  # Принимаем решение 2-го порядка
                    x_current += h

                    output_file.write(
                        f"{x_current:.6f} {y_current:.6e} {error_estimation:.6e}\n"
                    )  # Записываем в файл сразу после вычисления
                    integration_points_count += 1

                    if h * 2 <= (x_end - x_current) and h * 2 <= h_initial * 2:
                        h *= 2.0  # Увеличение шага вдвое для следующего шага, если это возможно и не превышает начальный шаг удвоенный.
                else:
                    if h <= h_min:
                        min_step_count += 1  # Если вышли из-за минимального шага, тоже считаем минимальный шаг
                        y_current = y_rk2_step  # Даже если точность не достигнута, все равно двигаемся вперед с минимальным шагом, используя 2й порядок
                        x_current += h
                        output_file.write(
                            f"{x_current:.6f} {y_current:.6e} {error_estimation:.6e}\n"
                        )  # Записываем в файл даже если ошибка больше допуска
                        integration_points_count += 1
                        break  # Выходим из внешнего цикла, если достигли мин шага и точность не достигнута

    except Exception as e:
        print(f"Ошибка при записи в файл '{output_filename}': {e}")

    return integration_points_count, accuracy_failed_points_count, min_step_count


def read_input_parameters_from_file(filename="input.txt"):
    """
    Читает параметры задачи Коши из файла.
    Файл должен содержать каждый параметр на новой строке в следующем порядке:
    x_start, x_end, y0, h_initial.
    """
    try:
        with open(filename, "r") as f:
            lines = f.readlines()
            x_start = float(lines[0].strip())
            x_end = float(lines[1].strip())
            y0 = float(lines[2].strip())
            h_initial = float(lines[3].strip())
        return x_start, x_end, y0, h_initial
    except FileNotFoundError:
        print(f"Ошибка: Файл '{filename}' не найден.")
        return None, None, None, None
    except ValueError:
        print(
            f"Ошибка: Неверный формат данных в файле '{filename}'. Убедитесь, что x_start, x_end, y0, h_initial - числа и каждый на новой строке."
        )
        return None, None, None, None
    except IndexError:
        print(
            f"Ошибка: Недостаточно параметров в файле '{filename}'. Ожидаются x_start, x_end, y0, h_initial."
        )
        return None, None, None, None


def write_statistics_to_file(
    filename, integration_points_count, accuracy_failed_points_count, min_step_count
):
    """
    Дописывает статистику решения задачи Коши в конец файла.
    Предполагается, что файл уже открыт и в него записаны результаты интегрирования.
    """
    try:
        with open(filename, "a") as f:  # Открываем файл в режиме добавления ('a')
            f.write(
                f"{integration_points_count} {accuracy_failed_points_count} {min_step_count}\n"
            )
        print(f"Статистика дописана в файл '{filename}'")
    except Exception as e:
        print(f"Ошибка при записи статистики в файл '{filename}': {e}")


if __name__ == "__main__":
    # Имя файла для входных параметров
    input_filename = "input.txt"
    # Имя файла для выходных результатов
    output_filename = "output.txt"

    # Чтение параметров из файла
    x_start_example, x_end_example, y0_example, h_initial_example = (
        read_input_parameters_from_file(input_filename)
    )

    if x_start_example is not None:  # Проверка, что параметры были успешно прочитаны
        # Пример использования: dy/dx = y - x^2 + 1, y(0) = 0.5
        def f_example(x, y):
            return y - x**2 + 1

        tolerance_example = 1e-5  # Заданная точность

        # Решение задачи Коши с адаптивным шагом и записью в файл в реальном времени
        integration_points, accuracy_failed_points, min_steps = (
            solve_cauchy_adaptive_step(
                f_example,
                y0_example,
                x_start_example,
                x_end_example,
                h_initial_example,
                tolerance_example,
                output_filename,
            )
        )

        # Запись статистики в конец файла
        write_statistics_to_file(
            output_filename, integration_points, accuracy_failed_points, min_steps
        )

        print("\nСтатистика:")
        print(f"Число точек интегрирования: {integration_points}")
        print(f"Число точек, где не достигнута точность: {accuracy_failed_points}")
        print(f"Общее количество минимальных шагов интегрирования: {min_steps}")

        print(
            f"Результаты интегрирования и статистика записаны в файл '{output_filename}'"
        )
