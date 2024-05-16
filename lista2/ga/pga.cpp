#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <fstream>
#include <map>
#include <random>
#include <mpi.h>

using namespace std;

#define GEN_NUM 1000

class TSP {
private:
    static const int POP_SIZE = 1500;
    const double mutation_probability = 0.15;
    vector<string> cities;
    int CITIES;

    // Structure of a genome
    struct individual {
        vector<string> genome;
        int fitness;
    };

    vector<individual> population;

public:
    TSP(const vector<string> &c) {
        population.resize(POP_SIZE);
        this->cities = c;
        this->CITIES = c.size();
    }

    int rand_num(int start, int end) {
        int r = end - start;
        int rnum = start + rand() % r;
        return rnum;
    }

    bool repeat(string s, char ch) {
        for (int i = 0; i < s.size(); i++) {
            if (s[i] == ch)
                return true;
        }
        return false;
    }

    vector<individual> tournament_selection(int num_parents) {
        vector<individual> parents(num_parents);
        int tournament_size = 5;

        for (int i = 0; i < num_parents; i++) {
            vector<individual> tournament(tournament_size);
            for (int j = 0; j < tournament_size; j++) {
                int ind = rand() % POP_SIZE;
                tournament[j] = population[ind];
            }
            sort(tournament.begin(), tournament.end(), lessthan);
            parents[i] = tournament[0];
        }

        return parents;
    }

    void mutate(vector<string>& genome) {
        int num_swaps = genome.size() * 0.1;  // Oblicz 10% długości genomu

        for (int i = 0; i < num_swaps; i++) {
            // Losowo wybierz dwa różne indeksy w genomie
            int index1 = rand() % genome.size();
            int index2 = rand() % genome.size();
            while (index1 == index2) {  // Upewnij się, że indeksy są różne
                index2 = rand() % genome.size();
            }

            // Zamień miejscami miasta na wybranych indeksach
            swap(genome[index1], genome[index2]);
        }
    }

    vector<string> create_genome() {
        vector<string> genome = this->cities;
        random_device rd;
        mt19937 g(rd());
        shuffle(genome.begin(), genome.end(), g);
        return genome;
    }

    int cal_fitness(vector<string> genome, map<pair<string, string>, int> &distances) {
        int f = 0;
        for (int i = 0; i < genome.size() - 1; i++) {
            pair<string, string> city_pair = make_pair(genome[i], genome[i + 1]);
            if (distances[city_pair] == 0)
                return INT_MAX;
            f += distances[city_pair];
        }
        pair<string, string> city_pair = make_pair(genome[genome.size() - 1], genome[0]);
        if (distances[city_pair] == 0)
            return INT_MAX;
        f += distances[city_pair];

        return f;
    }

    double cooldown(double temp) {
        return temp * 0.9945;
    }

    static bool lessthan(struct individual t1, struct individual t2) {
        return t1.fitness < t2.fitness;
    }

    void display_population() {
      for (int i = 0; i < POP_SIZE; i++) {
          cout << "Genome: ";
          for (const auto& city : population[i].genome) {
              cout << city << ' ';
          }
          cout << ", Fitness: " << population[i].fitness << endl;
      }
    }

    vector<string> crossover(vector<string> parent1, vector<string> parent2) {
        int crossPoint = rand_num(1, this->CITIES - 1);

        // Copy first part of parent1
        vector<string> child(parent1.begin(), parent1.begin() + crossPoint);

        // Copy unique elements from parent2
        for (int i = 0; i < parent2.size(); i++) {
            if (find(child.begin(), child.end(), parent2[i]) == child.end()) {
                child.push_back(parent2[i]);
            }
        }

        return child;
    }

    vector<individual> roulette_selection(int num_parents) {
        vector<individual> parents(num_parents);
        int total_fitness = 0;
        for (int i = 0; i < POP_SIZE; i++) {
            total_fitness += population[i].fitness;
        }

        for (int i = 0; i < num_parents; i++) {
            int pick = rand() % total_fitness;
            int current = 0;
            for (int j = 0; j < POP_SIZE; j++) {
                current += population[j].fitness;
                if (current > pick) {
                    parents[i] = population[j];
                    break;
                }
            }
        }

        return parents;
    }

    void genomeToInt(int* genomeIndexes, vector<string> genome) {
        for (int i = 0; i < genome.size(); i++) {
            genomeIndexes[i] = std::find(this->cities.begin(), this->cities.end(), genome[i]) - this->cities.begin();
        }
    }

    vector<string> intToGenome(const int* genomeIndexes, int size) {
        vector<string> genome(size);
        for (int i = 0; i < size; i++) {
            genome[i] = this->cities[genomeIndexes[i]];
        }
        return genome;
    }

    void TSPUtil(map<pair<string, string>, int> &distances) {
        // Generation Number
        int gen = 1;
        // Number of Gene Iterations
        int gen_thres = GEN_NUM;

        int world_rank, world_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        struct individual temp;

        // Populating the genome pool.
        for (int i = 0; i < POP_SIZE; i++) {
            temp.genome = create_genome();
            temp.fitness = cal_fitness(temp.genome, distances);
            population[i] = temp;
        }

        bool found = false;
        double temperature = 100;

        const int exchange_interval = 5;

        // Iteration to perform population crossing and gene mutation.
        while (temperature > 0.1 <= gen_thres) {
            sort(population.begin(), population.end(), lessthan);
            cout << "Generation: " << gen << endl;
            cout << "Best solution: ";
            for (const auto& city : population[0].genome) {
                cout << city << ' ';
            }
            cout << ", Fitness: " << population[0].fitness << endl;

            int elite_size = 5;
            vector<individual> new_population(POP_SIZE);
            copy(population.begin(), population.begin() + elite_size, new_population.begin());

            for (int i = elite_size; i < POP_SIZE; i++) {
                vector<individual> parents = tournament_selection(2);
                struct individual p1 = parents[0];
                struct individual p2 = parents[1];

                vector<string> new_g = crossover(p1.genome, p2.genome);
                struct individual new_genome;
                if ((double)rand() / RAND_MAX < mutation_probability) {
                    mutate(new_g);
                }
                new_genome.genome = new_g;
                new_genome.fitness = cal_fitness(new_genome.genome, distances);

                new_population[i] = new_genome;
            }

            if (gen % exchange_interval == 0) {
                // Przygotuj dane do wysłania
                int send_count = elite_size * this->CITIES; // Liczba miast na genom * liczba genomów
                int* send_data = new int[send_count];

                // Przekształć genom na format liczbowy
                for (int i = 0; i < elite_size; i++) {
                    genomeToInt(send_data + i * this->CITIES, new_population[i].genome);
                }

                // Przygotuj bufor na odbierane dane
                int recv_count = send_count * world_size;
                int* recv_data = new int[recv_count];

                // Wykonaj Allgather
                MPI_Allgather(send_data, send_count, MPI_INT, recv_data, send_count, MPI_INT, MPI_COMM_WORLD);

                // Znajdź najgorsze rozwiązania w populacji
                sort(new_population.begin(), new_population.end(), lessthan);
                int worst_index = POP_SIZE - recv_count / this->CITIES;
                // Zastąp najgorsze rozwiązania najlepszymi z recv_data
                for (int i = 0; i < recv_count; i += this->CITIES) {
                    vector<int> genomeInt(recv_data + i, recv_data + i + this->CITIES);
                    vector<string> genome = intToGenome(genomeInt.data(), this->CITIES);
                    int fitness = cal_fitness(genome, distances);

                    // Zastąp najgorsze rozwiązanie (ostatnie w posortowanej populacji)
                    int worst_index = new_population.size() - 1 - i;
                    cout << "ZASTEPOWANIE " << new_population[worst_index].fitness << " przez " << fitness << endl;
                    new_population[worst_index].genome = genome;
                    new_population[worst_index].fitness = fitness;
                }

                // Pamiętaj o zwolnieniu pamięci
                delete[] send_data;
                delete[] recv_data;
            }

            //temperature = cooldown(temperature);
            population = new_population;
            gen++;
        }
    }
};

void display_map(const map<pair<string, string>, int> &distances) {
    for (const auto &pair : distances) {
        cout << "Odległość między " << pair.first.first << " i " << pair.first.second << " wynosi " << pair.second << endl;
    }
}

void load_map_from_file(map<pair<string, string>, int> &distances, const std::string& filename, vector<string> cities) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Nie można otworzyć pliku: " << filename << std::endl;
        return;
    }

    vector<vector<int>> distance_matrix(cities.size(), vector<int>(cities.size()));
    for (int i = 0; i < cities.size(); ++i) {
        for (int j = 0; j < cities.size(); ++j) {
            file >> distance_matrix[i][j];
        }
    }

    for (int i = 0; i < cities.size(); ++i) {
        for (int j = i + 1; j < cities.size(); ++j) {
            distances[{cities[i], cities[j]}] = distance_matrix[i][j];
            distances[{cities[j], cities[i]}] = distance_matrix[j][i];
        }
    }
}

void load_cities_from_file(vector<string> &cities, const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Nie można otworzyć pliku: " << filename << std::endl;
        return;
    }

    string city;
    while (getline(file, city)) {
        cities.push_back(city);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    vector<string> cities;
    load_cities_from_file(cities, "codes30.txt");
    std::map<std::pair<std::string, std::string>, int> distances;
    load_map_from_file(distances, "dist30.txt", cities);
    //display_map(distances);
    TSP tsp(cities);

    tsp.TSPUtil(distances);
    MPI_Finalize();
}