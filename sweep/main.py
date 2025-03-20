import numpy as np
import math

def parse_possible_ln_or_num(input):
  if not input.startswith("ln"):
    return float(input)
  else:
    return math.log(float(input[2:]))
  
def read_vars(filename):
  with open(filename, mode='r') as file:
    readline = file.readlines()    

    ## constraint section
    constraints = readline[0].split()
    # print(alfa1, alfa2, gamma2)

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

def solve_cauchy(f, x_current, ):
  pass

def sweep():
  pass

if __name__ == "__main__":
  fname = "test-1.txt"

  constraints, limit1, limit2, amount, grid = read_vars(fname)

  alfa1, beta1, gamma1, alfa2, beta2, gamma2 = constraints