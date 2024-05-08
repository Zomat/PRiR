/**
 * Komputerowo Zintegrowane Wytwarzanie - Laboratorium 2
 * Program rozwiązujący problem szeregowania zadań z określonymi czasami przygotowania, wykonania i dostarczenia.
 * Autor: Mateusz Zolisz
 */
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <limits.h>
#include <fstream>
#include <map>
#include <vector>
#include <chrono>
#include <algorithm>

using namespace std;

int N;
vector<int> p, d, w;
vector<int> f;

vector<bool> toBinary(int n, int length) {
    vector<bool> binary(length, 0);
    int index = length - 1;
    while (n > 0) {
        binary[index--] = (n % 2 == 1);
        n /= 2;
    }
    return binary;
}

string boolVectorToString(vector<bool> const& v) {
    string s = "";
    for (bool b : v) {
        s += (b ? "1" : "0");
    }
    return s;
}

int C(int i) {
    int c = 0;
    for (int j = 0; j < N; j++) {
        if ((i & (1 << j))) {
            c += p[j];
        }
    }
    return c;
}

map<pair<int, int>, int> memoK;

int K(int a, int b) {
    pair<int, int> args = make_pair(a, b);
    if (memoK.count(args)) {
      return memoK[args];
    } else {
        int result = w[a] * max(0, (b - d[a]));
        memoK[args] = result;
        return result;
    }
}

vector<int> getBestOrder(int i, int c, int *minF, bool debug);

void process(int i, int c, int *minF, vector<int>* bestOrder, bool debug = false) {
    if (minF == nullptr) {
        *minF = numeric_limits<int>::max();
    }

    int tempF, tempI;
    vector<int> candidates;
    for (int j = N - 1; j >= 0; j--) {
        if ((i & (1 << j))) {
            tempI = i ^ (1 << j);
            vector<bool> tempBinI = toBinary(tempI, N);
            tempF = f[tempI] + K(j, c);
            if (debug) {
                cout << "\t \t F(" << boolVectorToString(tempBinI) << ") + K_" << j + 1 << "(" << c << ") =  " << tempF << endl;
            }

            if (tempF < *minF) {
                *minF = tempF;
                candidates.clear();
                candidates.push_back(j);
            } else if (tempF == *minF) {
                candidates.push_back(j);
            }
        }
    }

    if (bestOrder != nullptr && !candidates.empty()) {
        sort(candidates.begin(), candidates.end(), greater<int>());
        *bestOrder = getBestOrder(i ^ (1 << candidates[0]), C(i ^ (1 << candidates[0])), minF, debug);
        bestOrder->push_back(candidates[0]);
    }
}

vector<int> getBestOrder(int i, int c, int *minF, bool debug = false) {
    vector<int> bestOrder;

    if (i == 0) {
        return bestOrder;
    }

    process(i, c, minF, &bestOrder, debug);

    return bestOrder;
}

void loadFromFile(const string& filename, int dataNumber = 4) {
    ifstream file(filename);
    if (!file) {
        cerr << "Nie można otworzyć pliku " << filename << endl;
        return;
    }

    string s;
    string jString = to_string(dataNumber);

    do {
      file >> s;
    } while (s.compare("data." + jString + ":"));

    file >> N;
   
    int x, y, z;
    for (int i = 0; i < N; i++) {
      file >> x >> y >> z;
      p.push_back(x);
      w.push_back(y);
      d.push_back(z);
    }

    f.resize(1 << N);
    file.close();
}

int main(int argc, char* argv[]) {
    loadFromFile("data.txt", stoi(argv[1]));

    auto start1 = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < (1 << N); i++) {
        vector<bool> binI = toBinary(i, N);
        int c = C(i);
        cout << "C(" << boolVectorToString(binI) << ") = " << c << endl;

        if (i == 0) {
            f[i] = 0;
        } else {
            int minF = numeric_limits<int>::max();
            process(i, c, &minF, nullptr, true);

            f[i] = minF;
        } 

        cout << "F(" <<  boolVectorToString(binI) << ") =  " << f[i] << endl;
        cout << "-----------------------------------" << endl;
    }

    auto end1 = chrono::high_resolution_clock::now();
    auto start2 = std::chrono::high_resolution_clock::now();
    int tempMinF = numeric_limits<int>::max();
    vector<int> bestOrder = getBestOrder(((1 << N) - 1), C(((1 << N) - 1)), &tempMinF, false);

    cout << "Optymalne F: " << f[((1 << N) - 1)] << '\n';
    cout << "Permutacja: ";
    for (int task : bestOrder) {
        cout << task + 1 << ' ';
    }
    cout << '\n';

    auto end2 = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed1 = end1 - start1;
    chrono::duration<double> elapsed2 = end2 - start2;
    cout << "Czas wykonania algorytmu: " << elapsed1.count() << " sekund\n";
    cout << "Czas odtworzenia kolejnosci: " << elapsed2.count() << " sekund\n";
    
    return 0;
}