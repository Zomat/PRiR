#include <iostream>
#include <thread>
#include <mutex>
#include <chrono>
#include <semaphore.h>

const int k = 10;  // liczba palaczy
const int l = 4;  // liczba ubijaczy
const int m = 3;  // liczba pudełek zapałek

sem_t fillers;
sem_t matchbox;
std::mutex coutMutex;

void message(const std::string& msg) {
    coutMutex.lock();
    std::cout << msg << std::endl;
    coutMutex.unlock();
}

void smoker(int id) {
    while (true) {
        message("Palacz " + std::to_string(id) + " czeka na ubijacza.");
        sem_wait(&fillers);

        message("Palacz " + std::to_string(id) + " używa ubijacza.");
        std::this_thread::sleep_for(std::chrono::seconds(1));
        message("Palacz " + std::to_string(id) + " oddaje ubijacz.");

        sem_post(&fillers);
        
        message("Palacz " + std::to_string(id) + " czeka na pudełko zapałek.");
        sem_wait(&matchbox);
        message("Palacz " + std::to_string(id) + " zapala fajkę.");
        std::this_thread::sleep_for(std::chrono::seconds(1));
        message("Palacz " + std::to_string(id) + " oddaje pudełko zapałek.");
        sem_post(&matchbox);
        
        message("Palacz " + std::to_string(id) + " pali fajkę.");
        std::this_thread::sleep_for(std::chrono::seconds(2));
    }
}

int main() {
    sem_init(&fillers, 0, l);
    sem_init(&matchbox, 0, m);

    std::thread smokers[k];
    for (int i = 0; i < k; ++i) {
        smokers[i] = std::thread(smoker, i);
    }

    for (int i = 0; i < k; ++i) {
        smokers[i].join();
    }

    sem_destroy(&fillers);
    sem_destroy(&matchbox);

    return 0;
}