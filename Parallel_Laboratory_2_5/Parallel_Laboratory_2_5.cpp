// Parallel_Laboratory_2_5.cpp
// Лабораторная работа №2, Вариант №5
// Тема: Работа с библиотекой MPI. Написать MPI-приложение для расчета числа PI по формуле ряда: 1 - 1/3 + 1/5 - 1/7 + 1/9 - 1/11 + 1/13 - ...
// Автор: Мартынов В. В., группа ТИИ-2
// Запуск программы:
//   mpiexec -n <число_процессов> <путь_к_exe> <число_итераций>
// Пример:
//   mpiexec -n 4 "c:\...\parallel_laboratory_2_5.exe" 10000000

#define _USE_MATH_DEFINES   // Позволяет использовать константы вроде M_PI
#include <mpi.h>
#include <iostream>
#include <cmath>            // математические функции и константы
#include <iomanip>          // для форматированного вывода (setprecision)

using namespace std;

// Функция вычисления части ряда Лейбница для расчета числа PI
double compute_partial_pi(long long start, long long end) {
    double local_sum = 0.0;
    for (long long i = start; i < end; ++i) {
        double term = ((i % 2 == 0) ? 1.0 : -1.0) / (2.0 * i + 1.0);
        local_sum += term;
    }
    return local_sum;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);                     // инициализация MPI среды

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // определяем номер текущего процесса
    MPI_Comm_size(MPI_COMM_WORLD, &size);       // определяем общее число процессов

    // Проверка аргументов
    if (argc < 2) {
        if (rank == 0)
            cerr << "Usage: mpirun -np <processes> ./mpi_pi_series <number_of_iterations>\n";

        // Завершение работы MPI
        MPI_Finalize();
        return 1;
    }

    // Преобразуем аргумент командной строки в число
    long long n = atoll(argv[1]);
    double start_time = MPI_Wtime();

    // Распределение итераций между процессами
    long long chunk = n / size;
    long long start = rank * chunk;
    long long end = (rank == size - 1) ? n : start + chunk;

    // Каждый процесс вычисляет свою часть
    double local_sum = compute_partial_pi(start, end);

    // Сбор всех частичных сумм
    double global_sum = 0.0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();

    // Главный процесс выводит результат
    if (rank == 0) {
        double pi = 4.0 * global_sum;
        cout << fixed << setprecision(12);
        cout << "Computed PI = " << pi << endl;
        cout << "Error      = " << fabs(M_PI - pi) << endl;
        cout << "Iterations = " << n << endl;
        cout << "Processes  = " << size << endl;
        cout << "Execution time: " << (end_time - start_time) << " seconds\n";
    }

    // Завершение работы MPI
    MPI_Finalize();
    return 0;
}