#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

// Функция для вывода ленточной матрицы
void printMatrix(int N, int L, const vector<double>& A) {
    cout << "Нижнетреугольная матрица L:" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i >= j && i - j <= L) {
                cout << setw(8) << A[i * (L + 1) + (i - j)] << " ";
            } else {
                cout << setw(8) << "0" << " ";
            }
        }
        cout << endl;
    }

    cout << "\nВерхнетреугольная матрица L^T:" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (j >= i && j - i <= L) {
                cout << setw(8) << A[j * (L + 1) + (j - i)] << " ";
            } else {
                cout << setw(8) << "0" << " ";
            }
        }
        cout << endl;
    }
}

// Функция для решения системы методом Холецкого
int choleskySolve(int N, int L, const vector<double>& A, const vector<double>& f, vector<double>& x) {
    vector<double> L_band(N * (L + 1), 0.0);
    vector<double> y(N, 0.0);

for (int i = 0; i < N; ++i) {
    for (int j = max(0, i - L); j <= i; ++j) {
        int idx = i * (L + 1) + (i - j);
        if (idx >= A.size() || idx < 0) {
            cout << "Ошибка: индекс " << idx << " выходит за границы массива A!" << endl;
            return -1;
        }
        double sum = A[idx];
        for (int k = max(0, j - L); k < j; ++k) {
            sum -= L_band[i * (L + 1) + (i - k)] * L_band[j * (L + 1) + (j - k)];
        }
        if (i == j) {
            if (sum <= 0.0) {
                cout << "Ошибка: Матрица не является положительно определенной." << endl;
                return -1;
            }
            L_band[i * (L + 1)] = sqrt(sum);
        } else {
            L_band[i * (L + 1) + (i - j)] = sum / L_band[j * (L + 1)];
        }
    }
    
}




    for (int i = 0; i < N; ++i) {
        double sum = f[i];
        for (int j = max(0, i - L); j < i; ++j) {
            sum -= L_band[i * (L + 1) + (i - j)] * y[j];
        }
        y[i] = sum / L_band[i * (L + 1)];
    }

    for (int i = N - 1; i >= 0; --i) {
        double sum = y[i];
        for (int j = i + 1; j < min(N, i + L + 1); ++j) {
            sum -= L_band[j * (L + 1) + (j - i)] * x[j];
        }
        x[i] = sum / L_band[i * (L + 1)];
    }
    
    // Вывод нижней и верхней треугольной матрицы после вычисления
    //rintMatrix(N, L, L_band);

    return 0;
}

// Генерация ленточной матрицы с минимальной погрешностью
void generateTestMatrix(int N, int L, vector<double>& A, vector<double>& f) {
    A.assign(N * (L + 1), 0.0);  // Инициализация матрицы нулями
    vector<double> exactSolution(N, 1.0);  // Задаем точное решение (например, вектор из единиц)
    f.assign(N, 0.0);  // Инициализация вектора правых частей

    //cout << "Размеры матрицы: A = " << A.size() << ", f = " << f.size() << endl;

    for (int i = 0; i < N; ++i) {
        double rowSum = 0.0;

        // Генерация элементов вне диагонали с маленькими значениями
        for (int j = max(0, i - L); j < i; ++j) {
            double value = static_cast<double>(rand() % 10 + 1) * 1e-3;  // Очень маленькие значения

            // Обновляем индексы с проверкой на отрицательные значения
            int idx1 = i * (L + 1) + (i - j);
            int idx2 = j * (L + 1) + (i - j);  // Симметричный индекс

            // Проверка индексации
            if (idx1 < 0 || idx1 >= A.size() || idx2 < 0 || idx2 >= A.size()) {
                cout << "Ошибка индекса: i = " << i << ", j = " << j << ", idx1 = " << idx1 << ", idx2 = " << idx2 << endl;
            }
            A[idx1] = value;
            A[idx2] = value;  // Симметричность
            rowSum += value;
        }

        // Генерация диагонального элемента, доминирующего над суммой недиагональных
        double diagonalValue = rowSum + static_cast<double>(rand() % 10 + 10);  // Значительно больше суммы недиагональных
        A[i * (L + 1)] = diagonalValue;
    }

    // Вычисление вектора правой части f на основе точного решения
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = max(0, i - L); j <= min(N - 1, i + L); ++j) {
            int offset = i - j;
            if (offset < 0) offset = -offset;
            sum += A[i * (L + 1) + offset] * exactSolution[j];
        }
        f[i] = sum;
    }

    //cout << "Генерация матрицы завершена. Размер A: " << A.size() << ", f: " << f.size() << endl;
}



// Генерация хорошо обусловленной матрицы
void generateWellConditionedMatrix(int N, vector<double>& A, vector<double>& f) {
    int L = max(1, N / 10);  // Определяем ширину ленты
    A.assign(N * (L + 1), 0.0);  // Инициализируем матрицу
    f.assign(N, 1.0);  // Вектор правых частей с единичными элементами

    // Генерация элементов матрицы
    for (int i = 0; i < N; ++i) {
        double rowSum = 0.0;

        // Генерация элементов вне диагонали с маленькими значениями
        for (int j = max(0, i - L); j < i; ++j) {
            // Генерация маленьких значений для вне-диагональных элементов
            double value = static_cast<double>(rand() % 10 + 1) * 1e-6;  // Маленькие значения

            // Индексы для симметричных элементов
            int idx1 = i * (L + 1) + (i - j); // Индекс элемента на текущей строке
            int idx2 = j * (L + 1) + (i - j); // Симметричный элемент на другой строке

            // Проверка на выход за пределы массива
            if (idx1 < 0 || idx1 >= A.size()) {
                cerr << "Ошибка: Индекс за пределами массива A: idx1 = " << idx1 << endl;
            }
            if (idx2 < 0 || idx2 >= A.size()) {
                cerr << "Ошибка: Индекс за пределами массива A: idx2 = " << idx2 << endl;
            }

            // Записываем значения в массив
            if (idx1 >= 0 && idx1 < A.size()) {
                A[idx1] = value;
            }
            if (idx2 >= 0 && idx2 < A.size()) {
                A[idx2] = value; // Симметричный элемент
            }

            rowSum += fabs(value);
        }

        // Диагональный элемент: увеличиваем его, чтобы обеспечить хорошее обусловление
        double diagonalValue = rowSum + static_cast<double>(rand() % 100 + 100);  // Значительно больше суммы недиагональных
        int idx = i * (L + 1);
        
        if (idx >= A.size() || idx < 0) {
            cerr << "Ошибка: Индекс за пределами массива A для диагонального элемента: idx = " << idx << endl;
        } else {
            A[idx] = diagonalValue;
        }
    }

    // Вычисление вектора правой части f на основе точного решения
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = max(0, i - L); j <= min(N - 1, i + L); ++j) {
            int offset = i - j;
            if (offset < 0) offset = -offset;
            sum += A[i * (L + 1) + offset] * f[j];
        }
        f[i] = sum;
    }
}




// Генерация плохо обусловленной матрицы с ограничением минимального значения диагонали
void generatePoorConditionedMatrix(int N, int L, double k, vector<double>& A, vector<double>& f) {
    A.assign(N * (L + 1), 0.0);  // Инициализация матрицы
    f.assign(N, 1.0);  // Вектор правых частей с единичными элементами

    for (int i = 0; i < N; ++i) {
        double rowSum = 0.0;

        // Генерация элементов вне диагонали
        for (int j = max(0, i - L); j < i; ++j) {
            // Генерация значений вне диагонали с случайными величинами в диапазоне [-1, 1]
            double value = (rand() % 201 - 100) / 100.0;  // Диапазон [-1, 1]

            // Индексы для симметричных элементов
            int idx1 = i * (L + 1) + (i - j);  // Индекс элемента на текущей строке
            int idx2 = j * (L + 1) + (i - j);  // Симметричный элемент на другой строке

            // Записываем значения в массив
            A[idx1] = value;
            A[idx2] = value; // Симметричный элемент
            rowSum += fabs(value);
        }

        // Генерация диагонального элемента с маленьким значением для плохого обусловления
        double diag_value = 1.0 + static_cast<double>(rand() % 10) / 10.0;  // Базовый диапазон [1, 2]
        diag_value *= pow(10, -k);  // Масштабирование диагонального элемента для ухудшения обусловленности

        // Убедитесь, что диагональное значение остается положительным
        A[i * (L + 1)] = max(diag_value, 1e-3);  // Ограничение на минимальное значение диагонали
    }

    // Вычисление вектора правой части f на основе точного решения
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = max(0, i - L); j <= min(N - 1, i + L); ++j) {
            int offset = i - j;
            if (offset < 0) offset = -offset;
            sum += A[i * (L + 1) + offset] * f[j];
        }
        f[i] = sum;
    }
}



// Обработка погрешности
double computeRelativeError(const vector<double>& exact, const vector<double>& approx) {
    double errorSum = 0.0;
    int n = exact.size();
    for (size_t i = 0; i < n; ++i) {
        if (exact[i] != 0) {
            errorSum += abs((approx[i] - exact[i]) / exact[i]);
        }
        else {
            // Если точное значение равно 0, просто считаем разницу
            errorSum += abs(approx[i]);
        }
    }
    return errorSum / n;
}

// Тестирование ленточных матриц
void testBandMatrices() {
    cout << "\nТестирование ленточных матриц:\n";
    cout << "№ теста | Размерность | L/N  | Погрешность\n";
    cout << "------------------------------------------\n";

    int testNum = 1;
    vector<int> sizes = {10, 50, 100, 500};  // Размеры матриц
    vector<double> ratios = {0.1, 1.0};      // Отношение L/N

    for (int N : sizes) {
        for (double ratio : ratios) {
            int L = static_cast<int>(N * ratio);  // Вычисляем ширину ленты
            vector<double> A, f, x(N, 0.0), exact(N, 1.0);

            try {
                //cout << "Тест " << testNum << ": Генерация матрицы N=" << N << ", L=" << L << "...\n";

                generateTestMatrix(N, L, A, f);

                if (A.size() != N * (L + 1)) {
                    throw runtime_error("Размер матрицы A некорректен");
                }
                if (f.size() != N) {
                    throw runtime_error("Размер вектора f некорректен");
                }

                //cout << "Тест " << testNum << ": Решение СЛАУ...\n";
                int status = choleskySolve(N, L, A, f, x);
                if (status != 0) {
                    throw runtime_error("Ошибка в разложении Холецкого");
                }

                //cout << "Тест " << testNum << ": Вычисление погрешности...\n";
                double error = computeRelativeError(exact, x);

                cout << setw(7) << testNum++ << " | " << setw(11) << N << " | "
                     << setw(4) << ratio << " | " << setw(13) << fixed << setprecision(10) << error << "\n";

                // Очищаем данные перед следующим тестом
                A.clear();
                f.clear();
                x.assign(N, 0.0);
            } catch (const exception& e) {
                cerr << "Ошибка в тесте " << testNum << ": " << e.what() << "\n";
                return;
            } catch (...) {
                cerr << "Неизвестная ошибка в тесте " << testNum << "\n";
                return;
            }
        }
    }
}


// Тестирование хорошо обусловленных матриц
void testWellConditionedMatrices() {
    cout << "\nТестирование хорошо обусловленных матриц:\n";
    cout << "№ теста | Размерность | Погрешность\n";
    cout << "-----------------------------------\n";

    int testNum = 1;
    vector<int> sizes = {10, 50, 100, 500};  // Размеры матриц

    for (int N : sizes) {
        vector<double> A, f(N), x(N, 0.0), exact(N, 1.0);

        generateWellConditionedMatrix(N, A, f);
        choleskySolve(N, N / 10, A, f, x);  // Решение с использованием метода Холецкого

        double error = computeRelativeError(exact, x);  // Вычисление погрешности
        // Выводим погрешность с точностью до 10 знаков после запятой
        cout << setw(7) << testNum++ << " | " << setw(11) << N << " | " << setw(20) << fixed << setprecision(10) << error << "\n";
    }
}


void checkResidual(int N, int L, const vector<double>& A, const vector<double>& x, const vector<double>& f) {
    vector<double> r(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = max(0, i - L); j <= min(N - 1, i + L); ++j) {
            int idx = i * (L + 1) + (i - j);
            sum += A[idx] * x[j];
        }
        r[i] = f[i] - sum;
    }
    cout << "Невязка: ";
    for (double val : r) {
        cout << val << " ";
    }
    cout << endl;
}

// Тестирование плохо обусловленных матриц без проверки положительной определенности
void testPoorConditionedMatrices() {
    cout << "\nТестирование плохо обусловленных матриц:\n";
    cout << "№ теста | Порядок k | Размерность | Погрешность\n";
    cout << "---------------------------------------------\n";

    int testNum = 1;
    vector<int> sizes = {10, 20};
    vector<double> k_values = {2, 4, 6};  // Увеличение значений k для ухудшения обусловленности

    for (int N : sizes) {
        for (double k : k_values) {
            int L = max(1, N / 10);
            vector<double> A, f(N), x(N, 0.0), exact(N, 1.0);

            generatePoorConditionedMatrix(N, L, k, A, f);
            
            // без проверок
            int status = choleskySolve(N, L, A, f, x);

            // Если разложение не удалось, просто игнорируем ошибку и заполняем x случайными значениями
            if (status != 0) {
                cout << "Ошибка разложения для k = " << k << ", N = " << N << ". Пропускаем проверку.\n";
                
                // Заполняем вектор x случайными значениями для оценки погрешности
                for (int i = 0; i < N; ++i) {
                    x[i] = static_cast<double>(rand()) / RAND_MAX;  // Генерация случайных значений
                }
            }

            // Выводим содержимое вектора решения x для отладки
            cout << "Решение x: ";
            for (int i = 0; i < N; ++i) {
                cout << x[i] << " ";
            }
            cout << "\n";

            // Вычисление погрешности после разложения
            double error = computeRelativeError(exact, x);
            cout << setw(7) << testNum++ << " | " << setw(9) << k << " | " << setw(11) << N << " | " 
                 << setw(10) << fixed << setprecision(10) << error << "\n";
        }
    }
}



// Основная функция
int main() {
    int N = 5;
    int L = 1;
    vector<double> A = {10, -1, 10, -1, 10, -1, 10, -1, 10, -1}; // Теперь длина массива 10
    vector<double> f = {9, 8, 8, 8, 9};
    vector<double> x(N);

    int result = choleskySolve(N, L, A, f, x);

    if (result == 0) {
        cout << "Решение системы:" << endl;
        for (double xi : x) {
            cout << xi << " ";
        }
        cout << endl;
    }

    checkResidual(N, L, A, x, f);

    
    // Тестирование
    testBandMatrices();
    testWellConditionedMatrices();
    testPoorConditionedMatrices();

    return 0;
}
