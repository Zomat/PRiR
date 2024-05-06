#include <iostream>
#include <cmath>
#include <mpi.h>

double calculate_partial_sum(int start, int end) {
    double partial_sum = 0;
    for (int i = start; i <= end; ++i) {
        partial_sum += 1.0 / i;
    }
    return partial_sum;
}

void run(int n, int p) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        // Server Process
        double total_sum = 0;
        for (int i = 1; i < size; ++i) {
            double partial_sum;
            MPI_Recv(&partial_sum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum += partial_sum;
        }

        if (size == 1) {
            total_sum = calculate_partial_sum(1, n);
        }
        double gamma = total_sum - log(n);
        std::cout << "Stala Gamma Eulera: " << gamma << std::endl;
    } else {
        // Worker Processes
        int start = (n / (size - 1)) * (rank - 1) + 1;
        int end = std::min((n / (size - 1)) * rank, n);
        double partial_sum = calculate_partial_sum(start, end);
        MPI_Send(&partial_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int n = 1000000;
    if (argc >= 2) {
        n = std::atoi(argv[1]);
    }

    int p, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (p == 0) {
        std::cout << "Uruchomiono " << size << " procesow" << std::endl;
    }

    double start_time = MPI_Wtime();
    run(n, p);
    double end_time = MPI_Wtime();

    double max_end_time;
    MPI_Reduce(&end_time, &max_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (p == 0) {
        std::cout << "Czas wykonania: " << max_end_time - start_time << " sekund" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
