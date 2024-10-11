import time
import numpy as np

class MatrixSolver:
    def __init__(self, dtype = np.float64):
        self.dtype = dtype
        self.shape = (10, 10)
        self.size = self.shape[0]
        self.k = 0

        self.qq = np.zeros(self.shape, dtype=self.dtype)
        self.f = np.zeros(self.shape[0], dtype=self.dtype)

        self.a = np.zeros(self.size-1, dtype=self.dtype)
        self.b = np.zeros(self.size,   dtype=self.dtype)
        self.c = np.zeros(self.size-1, dtype=self.dtype)
        self.p = np.zeros(self.size,   dtype=self.dtype)

    def scan_matrix(self, filename):
        with open(filename, encoding="utf-8", mode="r") as file:
            i = 0
            for line in file:
                line = line.split()
                if line:
                    self.f[i] = line[-1]
                    # test.append(line[:-1])

                    self.qq[i,:] = line[:-1] 
                    i += 1
            print(self.qq)

    def input_matrix(self, arr, f):
        self.set_matrix(arr)
        self.set_f(f)

    def set_matrix(self, arr):
        self.qq = np.copy(arr)

        self.shape = self.qq.shape
        self.size = self.shape[0]

        self.f = np.zeros(self.size, dtype=self.dtype)

        self.a = np.zeros(self.size-1, dtype=self.dtype)
        self.b = np.zeros(self.size,   dtype=self.dtype)
        self.c = np.zeros(self.size-1, dtype=self.dtype)
        self.p = np.zeros(self.size,   dtype=self.dtype)

    def set_f(self, f):
        self.f = np.copy(f) 

    def build_vectors(self):
        self.k = 1

        for i in range(1, self.size-1):
            col = np.delete(self.qq[:, i], [i-1, i, i+1])
            if np.any(col!=0):
                self.k = i

        # print("k = ", self.k)

        self.a = self.qq[np.arange(1, self.size), np.arange(self.size-1)]

        self.b = self.qq[np.arange(self.size), np.arange(self.size)]

        self.c = self.qq[np.arange(self.size - 1), np.arange(1, self.size)]

        self.p = self.qq[:, self.k].copy()


        # print("a = ", self.a)
        # print("b = ", self.b)
        # print("c = ", self.c)
        # print("p = ", self.p)
        # print("f = ", self.f)
        # print("---")

    def view_vectors(self):
        print("a = ", self.a)
        print("b = ", self.b)
        print("c = ", self.c)
        print("p = ", self.p)
        print("f = ", self.f)
        print("---")
    
    def view_matrix(self):
        print(self.qq)

    def subtract_from(self, from_idx, source_idx):
        if from_idx < source_idx:
            # print(from_idx, source_idx)

            # value = (-1.) * self.c[from_idx]
            value = self.c[from_idx]

            if source_idx == self.k + 2:
                self.a[from_idx-1] -= self.p[source_idx] * value

            self.c[from_idx] -= value
            self.b[from_idx] -= self.a[source_idx-1] * value
            self.p[from_idx] -= self.p[source_idx] * value
            self.f[from_idx] -= self.f[source_idx] * value
        else:
            # value = (-1.) * self.a[from_idx-1]
            value = self.a[from_idx-1]

            self.a[from_idx-1] -= value

            if source_idx == self.k - 2:
                self.c[from_idx] -= self.p[source_idx] * value

            self.b[from_idx] -= self.c[source_idx] * value 
            self.p[from_idx] -= self.p[source_idx] * value
            self.f[from_idx] -= self.f[source_idx] * value

    def gauss_down(self):
        for i in range(self.size-1):
            self.div_line_by(i, self.b[i])            
            self.subtract_from(i+1, i)

    def gauss_up(self):
        for i in range(self.size-1, 0, -1):
            self.subtract_from(i-1, i)

    def gauss_until_k(self):
        # print("gauss_until_k")
        for i in range(self.k):
            self.div_line_by(i, self.b[i])            
            self.subtract_from(i+1, i)

    def reverse_gauss_until_k(self):
        # print("reverse_gauss_until_k")
        for i in range(self.k-1, 0, -1):
            # print(i-1, i)
            self.subtract_from(i-1, i)

    def gauss_past_k(self):
        # print("gauss_past_k")
        for i in range(self.size-1, self.k, -1):
            self.div_line_by(i, self.b[i])            
            self.subtract_from(i-1, i)

    def reverse_gauss_past_k(self):
        # print("reverse_gauss_past_k")
        for i in range(self.k, self.size-1):
            # print(i+1, i)
            self.subtract_from(i+1, i)

    def gauss_line_k(self):
        # print("gauss_line_k")
        self.div_line_by(self.k, self.b[self.k])
        f_k = self.f[self.k]
        # print(f_k)

        self.a[self.k] = 0
        self.c[self.k-1] = 0
        for i in range(self.k):
            # value = (-1.) * self.p[i]
            value = self.p[i]
            self.p[i] -= value
            self.f[i] -= value * f_k

        for i in range(self.size-1, self.k, -1):
            # value = (-1.) * self.p[i]
            value = self.p[i]
            self.p[i] -= value
            self.f[i] -= value * f_k

    def div_line_by(self, line_idx, divisor):
        if divisor == 0:
            return None

        if line_idx > 0:
            self.a[line_idx - 1] /= divisor

        if line_idx < self.size - 1:
            self.c[line_idx] /= divisor

        self.b[line_idx] /= divisor
        self.p[line_idx] /= divisor
        self.f[line_idx] /= divisor

    def solve(self):
        solver.build_vectors()

        solver.gauss_until_k()
        # solver.view_vectors()

        solver.gauss_past_k()
        # solver.view_vectors()

        solver.gauss_line_k()
        # solver.view_vectors()

        solver.reverse_gauss_until_k()
        # solver.view_vectors()

        solver.reverse_gauss_past_k()
        # solver.view_vectors()
        return self.f

    def multiply_test(self):
        a = np.ones(self.size, dtype=self.dtype)
        # print(a)
        return np.matmul(self.qq , a)

    def multiply_test_aboba(self, x):
        return np.matmul(self.qq , x)
    
    def max_sub_ones(self):
        return np.max(np.abs(self.f - np.ones(self.size, dtype=self.dtype)))

class SpecialMatrixGen:

    def __init__(self, dtype=np.double):
        self.dtype = dtype 
   
    def gen_matrix(self, dim, low, high):
        """dimensionality = 10^dim"""
        rng = np.random.default_rng()
        size = rng.integers(10**dim, 10**(dim+1))
        k = rng.integers(1, size-1)
        shape = (size, size)
        arr = np.zeros(shape, dtype=self.dtype)
        f = rng.uniform(low, high, size)

        arr[np.arange(1, size), np.arange(size-1)] = rng.uniform(low, high, size-1)

        b_uniform = rng.uniform(low, high, size)
        # print(np.any(b_uniform==0.0), b_uniform)
        while np.any(b_uniform==0.0):
            b_uniform = rng.uniform(low, high, size)
        arr[np.arange(size), np.arange(size)] = b_uniform

        arr[np.arange(size - 1), np.arange(1, size)] = rng.uniform(low, high, size-1)
        arr[np.arange(size), k] = rng.uniform(low, high, size)
        
        # for i in range(1, size-1):
        #     col = np.delete(arr[:, i], [i-1, i, i+1])
        #     if np.any(col, where=lambda x: x!=0)):

        return arr, f

if __name__ == "__main__":
    start_time = time.time()
    np.set_printoptions(precision=5, linewidth=150, floatmode="unique", suppress=True)
    num_gen = np.random.default_rng()

    solver = MatrixSolver(np.double)
    gen = SpecialMatrixGen()

    solver.scan_matrix("matrix.txt")
    by_one = solver.multiply_test()

    f = solver.solve()
    print(f)

    solver.set_f(by_one)

    f = solver.solve()
    # q = maxx
    maxx = solver.max_sub_ones()
    print(f)
    print(maxx)

    # mat = gen.gen_matrix(1, -10, 10)

    result_2d= []

    for i in range(1, 4):
        for j in [10, 100, 1000]:
            result = []
            result.append(10**i)

            mat = gen.gen_matrix(i, -j, j)

            result.append(mat[0].shape)
            result.append((-j, j))

            solver.input_matrix(*mat)
            by_one = solver.multiply_test()
            solver.set_f(by_one)

            solver.solve()
            maxx = solver.max_sub_ones()
            # print("maxx = ", maxx)
            result.append(maxx)

            ## x = mat[1]
            random_f = solver.multiply_test_aboba(mat[1])
            solver.set_f(random_f)
            x_star = solver.solve()
            # print([mat[1][i] - x_star[i] if x_star[i] <= maxx else (mat[1][i] - x_star[i])/maxx for i in range(mat[1].shape[0])])
            last_maxx = np.max([np.abs(mat[1][i] - x_star[i]) if np.abs(x_star[i]) <= maxx else np.abs(mat[1][i] - x_star[i])/maxx for i in range(mat[1].shape[0])])
            # print("last maxx: ", last_maxx) 
            result.append(last_maxx)
            result_2d.append(result)

        #     break
        # break




    # solver.input_matrix(*mat) 
    # f = solver.solve()

    # print(f)
    for i in result_2d:
        print(i)

    print(f"--- {(time.time() - start_time)} seconds ---")

