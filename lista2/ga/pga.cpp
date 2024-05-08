#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <mpi.h>

using namespace std;

class TSP {
private:
    // Initial population size for the algorithm
    static const int POP_SIZE = 100;

    // Structure of a genome
    struct individual {
        string genome;
        int fitness;
    };

    vector<individual> population;

public:
    TSP() {
        population.resize(POP_SIZE);
    }

    static const int CITIES = 50;

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

    string mutatedGene(string genome) {
        while (true) {
            int r = rand_num(1, CITIES);
            int r1 = rand_num(1, CITIES);
            if (r1 != r) {
                char temp = genome[r];
                genome[r] = genome[r1];
                genome[r1] = temp;
                break;
            }
        }
        return genome;
    }

    string create_genome() {
        string genome = "0";
        while (true) {
            if (genome.size() == CITIES) {
                genome += genome[0];
                break;
            }
            int temp = rand_num(1, CITIES);
            if (!repeat(genome, (char)(temp + 48)))
                genome += (char)(temp + 48);
        }
        return genome;
    }

    int cal_fitness(string genome, int map[CITIES][CITIES]) {
        int f = 0;
        for (int i = 0; i < genome.size() - 1; i++) {
            if (map[genome[i] - 48][genome[i + 1] - 48] == INT_MAX)
                return INT_MAX;
            f += map[genome[i] - 48][genome[i + 1] - 48];
        }
        return f;
    }

    int cooldown(int temp) {
        return (90 * temp) / 100;
    }

    static bool lessthan(struct individual t1, struct individual t2) {
        return t1.fitness < t2.fitness;
    }

    void display_population() {
        for (int i = 0; i < POP_SIZE; i++) {
            cout << "genome: " << population[i].genome << ", Fitness: " << population[i].fitness << endl;
        }
    }

    void TSPUtil(int map[CITIES][CITIES], int rank, int worldSize) {
        // Generation Number
        int gen = 1;
        // Number of Gene Iterations
        int gen_thres = 5;

        struct individual temp;

        // Populating the genome pool.
        for (int i = rank; i < POP_SIZE; i += worldSize) {
            temp.genome = create_genome();
            temp.fitness = cal_fitness(temp.genome, map);
            population[i] = temp;
        }

        bool found = false;
        int temperature = 10000;

        // Iteration to perform population crossing and gene mutation.
        while (temperature > 1000 && gen <= gen_thres) {
            sort(population.begin(), population.end(), lessthan);
            cout << "Generation: " << gen << endl;
            display_population();
            vector<struct individual> new_population(POP_SIZE);

            for (int i = rank; i < POP_SIZE; i += worldSize) {
                struct individual p1 = population[i];

                while (true) {
                    string new_g = mutatedGene(p1.genome);
                    struct individual new_genome;
                    new_genome.genome = new_g;
                    new_genome.fitness = cal_fitness(new_genome.genome, map);

                    if (new_genome.fitness <= population[i].fitness) {
                        new_population[i] = new_genome;
                        break;
                    }
                    else {
                        // Accepting the rejected children at a possible probability above threshold.
                        float prob = pow(2.7, -1 * ((float)(new_genome.fitness - population[i].fitness) / temperature));
                        if (prob > 0.5) {
                            new_population[i] = new_genome;
                            break;
                        }
                    }
                }
            }

            temperature = cooldown(temperature);

            MPI_Allgather(&new_population[0], POP_SIZE, MPI_INT, &population[0], POP_SIZE, MPI_INT, MPI_COMM_WORLD);

            gen++;
        }
    }
};

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int worldSize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int map[TSP::CITIES][TSP::CITIES] = { { 0, 2, INT_MAX, 12, 5 },
                                { 2, 0, 4, 8, INT_MAX },
                                { INT_MAX, 4, 0, 3, 3 },
                                { 12, 8, 3, 0, 10 },
                                { 5, INT_MAX, 3, 10, 0 } };

    TSP tsp;
    tsp.TSPUtil(map, rank, worldSize);

    MPI_Finalize();
    return 0;
}
