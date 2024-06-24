#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <limits.h>
#include <fstream>
#include <map>
#include <vector>
#include <chrono>
#include <algorithm>
#ifdef __APPLE__
#include <__numeric/iota.h>
#else
#include <numeric>
#endif
#include <random>
#include <cstdio>
#include <mpi.h>

using namespace std;

int N, M;
vector<vector<int>> p;
vector<int> order;

void loadFromFile(const string& filename, string dataNumber = "000") {
    ifstream file(filename);
    if (!file) {
        cerr << "Nie można otworzyć pliku " << filename << endl;
        return;
    }

    string s;
    do {
      file >> s;
    } while (s.compare("data." + dataNumber + ":"));

    file >> N >> M;
    cout << "Maszyn: " << M << " Zadań: " << N << endl;
    p.resize(N, vector<int>(M));
    order.resize(N);
    int temp_p;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
        file >> temp_p;
        p[i][j] = temp_p;
      }
    }

    file.close();
}

void displayP(const vector<vector<int>>& p) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            cout << p[i][j] << " ";
        }
        cout << endl;
    }
}

void displayOrder(const vector<int>& order) {
    for (int i = 0; i < order.size(); i++) {
        cout << order[i] + 1 << " ";
    }
    cout << endl;
}

int calculateCMax(const vector<vector<int>>& p, const vector<int>& order) {
    int N = order.size();
    int M = p[0].size();
    vector<vector<int>> c(N + 1, vector<int>(M + 1, 0));

    for (int i = 1; i <= N; i++) {
        c[i][1] = c[i - 1][1] + p[order[i - 1]][0];
    }

    for (int i = 1; i <= M; i++) {
        c[1][i] = c[1][i - 1] + p[order[0]][i - 1];
    }

    for (int i = 2; i <= N; i++) {
        for (int j = 2; j <= M; j++) {
            c[i][j] = max(c[i - 1][j], c[i][j - 1]) + p[order[i-1]][j-1];
        }
    }

    return c[N][M];
}

void cooldown(int i, double *t, double factor) {
    *t = i % 100 == 0 ? *t * factor : *t;
}

double calculate_initial_temperature(const vector<double>& deltas, double target_acceptance_probability) {
    double mean_delta = accumulate(deltas.begin(), deltas.end(), 0.0) / deltas.size();
    return -mean_delta / log(target_acceptance_probability);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
        if (argc < 3) {
            cerr << "Użycie: " << argv[0] << " <data_number> <auto_tune>" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    string dataNumber;
    bool auto_tune;

    if (world_rank == 0) {
        dataNumber = argv[1];
        loadFromFile("data.txt", dataNumber);
        auto_tune = argv[2][0] == '1';
    }

    // Rozsyłanie flagi auto_tune do wszystkich procesów
    MPI_Bcast(&auto_tune, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    // Rozsyłanie danych wejściowych do wszystkich procesów
    if (world_rank == 0) {
        for (int i = 1; i < world_size; ++i) {
            MPI_Send(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&M, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            for (int j = 0; j < N; ++j) {
                MPI_Send(p[j].data(), M, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
    } else {
        MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&M, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        p.resize(N, vector<int>(M));
        for (int j = 0; j < N; ++j) {
            MPI_Recv(p[j].data(), M, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    vector<int> order(N);
    iota(order.begin(), order.end(), 0);

    random_device rd;
    mt19937 g(rd());

    // Strojenie początkowej temperatury
    double initial_temp_calibration = 1000;
    double cooling_factor_calibration = 1;
    vector<double> deltas;

    if (auto_tune) {
        for (int i = 0; i < 1000; i++) {
            vector<int> tempOrder = order;
            int index1 = uniform_int_distribution<int>(0, N - 1)(g);
            int index2 = uniform_int_distribution<int>(0, N - 1)(g);
            swap(tempOrder[index1], tempOrder[index2]);

            int accCmax = calculateCMax(p, tempOrder);
            int cmax = calculateCMax(p, order);
            double delta = accCmax - cmax;
            deltas.push_back(delta);

            if (delta < 0 || exp(-delta / initial_temp_calibration) > uniform_real_distribution<double>(0.0, 1.0)(g)) {
                order = tempOrder;
            }
        }
    }

    double target_acceptance_probability = 0.8;
    double initial_temp = auto_tune ? calculate_initial_temperature(deltas, target_acceptance_probability) : 2500;
    double cooling_factor = 0.9930;
    double min_temp = 0.0001;

    if (auto_tune && world_rank == 0) {
        cout << "Dostrojona początkowa temperatura: " << initial_temp << endl;
    }

    // Właściwa faza symulowanego wyżarzania
    double t = initial_temp;

    vector<int> best_order = order;
    int best_cmax = calculateCMax(p, order);

    for (int i = 0; i < 100000 && t > min_temp; i++) {
        cooldown(i, &t, cooling_factor);

        vector<int> tempOrder = order;
        int index1 = uniform_int_distribution<int>(0, N - 1)(g);
        int index2 = uniform_int_distribution<int>(0, N - 1)(g);
        swap(tempOrder[index1], tempOrder[index2]);

        if (i % 1000 == 0) {
            cout << "RANK: " << world_rank << " Iteracja: " << i << " Temp: " << t << " Kolejnosc: ";
            displayOrder(tempOrder);
        }

        int accCmax = calculateCMax(p, tempOrder);
        int cmax = calculateCMax(p, order);

        if (accCmax < cmax) {
            order = tempOrder;
            cmax = accCmax;
        } else {
            double acceptanceProbability = exp(-(accCmax - cmax) / t);
            uniform_real_distribution<double> distribution(0.0, 1.0);
            double randomValue = distribution(g);
            if (randomValue < acceptanceProbability) {
                order = tempOrder;
                cmax = accCmax;
            }
        }

        if (cmax < best_cmax) {
            best_cmax = cmax;
            best_order = order;

            ofstream output_file("cmax_values_rank_" + to_string(world_rank) + ".txt", ofstream::out | ofstream::app);
            if (output_file.is_open()) {
                output_file << best_cmax << endl;
                output_file.close();
            } else {
                cerr << "Nie można otworzyć pliku do zapisu dla ranka " << world_rank << endl;
            }
        }
    }

    int global_best_cmax;
    vector<int> global_best_order(N);

    MPI_Allreduce(&best_cmax, &global_best_cmax, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (best_cmax == global_best_cmax) {
        MPI_Send(best_order.data(), N, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    if (world_rank == 0) {
        for (int i = 0; i < world_size; ++i) {
            int recv_best_cmax;
            MPI_Status status;
            MPI_Recv(global_best_order.data(), N, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            recv_best_cmax = calculateCMax(p, global_best_order);
            if (recv_best_cmax == global_best_cmax) {
                break;
            }
        }

        cout << "Global CMAX: " << global_best_cmax << endl;
        cout << "Global ORDER: ";
        displayOrder(global_best_order);
    }

    MPI_Finalize();

    return 0;
}