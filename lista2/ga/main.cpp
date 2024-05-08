#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <fstream>
#include <map>
#include <random>
using namespace std;

#define GEN_NUM 1000

class TSP {
private:
    // Initial population size for the algorithm
    static const int POP_SIZE = 1000;
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

    vector<string> mutatedGene(vector<string> genome) {
        int num_swaps = rand_num(1, this->CITIES / 2);
        for (int i = 0; i < num_swaps; i++) {
            while (true) {
                int r = rand_num(1, this->CITIES);
                int r1 = rand_num(1, this->CITIES);
                if (r1 != r) {
                    swap(genome[r], genome[r1]);
                    break;
                }
            }
        }
        return genome;
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

    void TSPUtil(map<pair<string, string>, int> &distances) {
        // Generation Number
        int gen = 1;
        // Number of Gene Iterations
        int gen_thres = GEN_NUM;

        struct individual temp;

        // Populating the genome pool.
        for (int i = 0; i < POP_SIZE; i++) {
            temp.genome = create_genome();
            temp.fitness = cal_fitness(temp.genome, distances);
            population[i] = temp;
        }

        bool found = false;
        int temperature = 10000;

        // Iteration to perform population crossing and gene mutation.
        while (temperature > 1000 && gen <= gen_thres) {
            sort(population.begin(), population.end(), lessthan);
            cout << "Generation: " << gen << endl;
            cout << "Best solution: ";
            for (const auto& city : population[0].genome) {
                cout << city << ' ';
            }
            cout << ", Fitness: " << population[0].fitness << endl;

            int elite_size = POP_SIZE / 10; // 10% populacji
            vector<individual> new_population(POP_SIZE); // inicjalizacja nowej populacji z pełnym rozmiarem
            copy(population.begin(), population.begin() + elite_size, new_population.begin()); // skopiuj elitę

            for (int i = elite_size; i < POP_SIZE; i++) {
                vector<individual> parents = tournament_selection(2); // wybierz tylko dwóch rodziców
                struct individual p1 = parents[0];
                struct individual p2 = parents[1];

                vector<string> new_g = crossover(p1.genome, p2.genome);
                struct individual new_genome;
                new_genome.genome = new_g;
                new_genome.fitness = cal_fitness(new_genome.genome, distances);

                // Simulated annealing
                if (new_genome.fitness < population[i].fitness) {
                    new_population[i] = new_genome;
                } else {
                    double acceptance_prob = exp((population[i].fitness - new_genome.fitness) / temperature);
                    if (acceptance_prob > ((double) rand() / (RAND_MAX))) {
                        new_population[i] = new_genome;
                    } else {
                        new_population[i] = population[i];
                    }
                }
            }

            temperature = cooldown(temperature);
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

int main() {
    vector<string> cities;
    load_cities_from_file(cities, "codes30.txt");
    std::map<std::pair<std::string, std::string>, int> distances;
    load_map_from_file(distances, "dist30.txt", cities);
    //display_map(distances);
    TSP tsp(cities);
    tsp.TSPUtil(distances);
}