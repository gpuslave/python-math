import numpy as np
import math
import matplotlib.pyplot as plt

def parse_possible_ln_or_num(value):
  value = value.strip()  
  if value.startswith("ln"):
    return math.log(float(value[2:]))
  return float(value)

def read_vars(filename):
  with open(filename, mode='r') as file:
    readline = file.readlines()    

    ## constraint section
    constraints = list(map(parse_possible_ln_or_num, readline[0].split()))
    # print(constraints)

    ## limit section
    limits = readline[1].split()
    limit1 = parse_possible_ln_or_num(limits[0])
    limit2 = parse_possible_ln_or_num(limits[1])
    amount = int(limits[2])
    # print(limit1, limit2, amount)

    ## grid section
    grid = []
    for i in range(amount):
      temp = readline[i+2].split()[:2]
      temp = map(lambda x: x.replace(",", "."), temp)
      temp = list(map(lambda x: float(x), temp))
      grid.append(temp)
    # print(grid)

    return constraints, limit1, limit2, amount, grid

def runge_kutta_2(f_list, init_cond, x_start, h, x_end, amount, grid):
    # x = np.arange(x_start, x_end+h, h)
    # x = np.linspace(x_start, x_end, amount)
    x = grid
    y = np.zeros((len(x), len(init_cond)))

    y[0] = init_cond

    for i in range(len(x) - 1):
      h = x[i+1] - x[i]
      k1 = np.array([ h * func(x[i], y[i]) for func in f_list])
      k2 = np.array([ h * func(x[i] + 0.5*h, y[i] + 0.5*k1) for func in f_list])
      y[i+1] = y[i] + k2

    print()
    print(x, "\n", y)
    return x, y

def runge_kutta_2_reversed(f_list, init_cond, x_start, h, x_end, amount, grid):
    # x = np.arange(x_end, x_start-h, -h)
    # x = np.linspace(x_start, x_end, amount)[::-1]
    x = grid[::-1]
    y = np.zeros((len(x), len(init_cond)))

    y[0] = init_cond

    for i in range(len(x) - 1):
      h = x[i+1] - x[i]
      k1 = np.array([ -h * func(x[i], y[i]) for func in f_list])
      k2 = np.array([ -h * func(x[i] - 0.5*h, y[i] + 0.5*k1) for func in f_list])
      y[i+1] = y[i] + k2

    print()
    print(x, "\n", y)
    return x, y

def solve_linear_algebraic_system(u, v, w, alpha_coeff, beta_coeff, gamma_coeff):
    """Решает систему линейных алгебраических уравнений (3.10)."""
    det = u * (-beta_coeff) - (-v) * alpha_coeff
    if det == 0:
        return None, None  # Определитель равен нулю, нет единственного решения
    else:
        y_prime = (w * (-beta_coeff) - (-v) * gamma_coeff) / det
        y = (u * gamma_coeff - w * alpha_coeff) / det
        return y, y_prime

if __name__ == "__main__":
  fname = "test-1.txt"

  constraints, limit1, limit2, amount, grid = read_vars(fname)
  alpha1, beta1, gamma1, alpha2, beta2, gamma2 = constraints

  ## Шаг
  h = (limit2 - limit1) / (amount - 1)


  grid_x = [ temp[1] for temp in grid ]
  amount = len(grid_x)
  ## Вспомогательные функции
  def p_func(x): return -1
  def q_func(x): return 0
  def f_func(x): return 0

  ## Прогонка в первую сторону
  def part_u(x, y):
    return p_func(x) * y[0] + y[1]
  
  def part_v(x, y):
    return q_func(x) * y[0]
  
  def part_w(x, y):
    return f_func(x) * y[0]

  ## Разложение по компонентам
  f_list_1 = [part_u, part_v, part_w]
  initial_coniditions_1 = np.array([alpha1, -beta1, gamma1])

    

  x_vals_1, results_1 = runge_kutta_2(f_list_1, initial_coniditions_1, limit1, h, limit2, amount-1, grid_x)
  u_vals = results_1[:, 0]
  v_vals = results_1[:, 1]
  w_vals = results_1[:, 2]

  ## Прогонка во вторую сторону
  def part_alpha(x, y):
    return -p_func(x) * y[0] - y[1]
  
  def part_beta(x, y):
    return -q_func(x) * y[0]
  
  def part_gamma(x, y):
    return -f_func(x) * y[0]

  f_list_2 = [part_alpha, part_beta, part_gamma]
  initial_coniditions_2 = np.array([alpha2, -beta2, gamma2])

  x_vals_2, results_2 = runge_kutta_2_reversed(f_list_2, initial_coniditions_2, limit1, h, limit2, amount-1, grid_x)
  print(len(x_vals_2), len(x_vals_1))
  ## NOTE: надо ли развернуть results_2?
  x_vals_2 = x_vals_2[::-1]
  results_2 = results_2[::-1]
  # print(results_2)

  alpha_vals = results_2[:, 0]
  beta_vals = results_2[:, 1]
  gamma_vals = results_2[:, 2]

  y_solution = []
  y_prime_solution = []

  for i in range(len(x_vals_1)):
    u = u_vals[i]  
    v = v_vals[i]
    w = w_vals[i]
    alpha_coeff = alpha_vals[i]
    beta_coeff = beta_vals[i]
    gamma_coeff = gamma_vals[i]

    y_val, y_prime_val = solve_linear_algebraic_system(u, v, w, alpha_coeff, beta_coeff, gamma_coeff)
    if y_val is not None:
        y_solution.append(y_val)
        y_prime_solution.append(y_prime_val)
    else:
        print(f"Нет решения в точке x = {x_vals_1[i]:.4f}")
        y_solution.append(0)
        y_prime_solution.append(0)
  print(y_solution, "\n", y_prime_solution)
    

  # # Построение графиков (опционально)
  plt.figure(figsize=(10, 6))
  plt.plot(grid_x, y_solution, label='Численное решение')
  plt.xlabel('x')
  plt.ylabel('y')
  plt.title('Решение краевой задачи методом универсальной прогонки')
  plt.grid(True)
  plt.legend()

  plt.tight_layout()
  plt.savefig("solution.png", dpi=600)
  plt.close()