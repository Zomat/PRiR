#include <iostream>
#include <thread>
#include <mutex>
#include <chrono>
#include <semaphore.h>

const int k = 10;  // liczba palaczy
const int l = 4;  // liczba ubijaczy
const int m = 4;  // liczba pudełek zapałek

sem_t ubijacze;
sem_t pudełkaZapałek;
std::mutex coutMutex;

void komunikat(const std::string& msg) {
    coutMutex.lock();
    std::cout << msg << std::endl;
    coutMutex.unlock();
}

void palacz(int palaczID) {
    while (true) {
        komunikat("Palacz " + std::to_string(palaczID) + " czeka na ubijacza.");
        sem_wait(&ubijacze);
        komunikat("Palacz " + std::to_string(palaczID) + " używa ubijacza.");
        std::this_thread::sleep_for(std::chrono::seconds(1)); // symulacja czasu potrzebnego na ubijanie
        komunikat("Palacz " + std::to_string(palaczID) + " oddaje ubijacz.");
        sem_post(&ubijacze);
        
        komunikat("Palacz " + std::to_string(palaczID) + " czeka na pudełko zapałek.");
        sem_wait(&pudełkaZapałek);
        komunikat("Palacz " + std::to_string(palaczID) + " zapala fajkę.");
        std::this_thread::sleep_for(std::chrono::seconds(1)); // symulacja czasu potrzebnego na zapalenie fajki
        komunikat("Palacz " + std::to_string(palaczID) + " oddaje pudełko zapałek.");
        sem_post(&pudełkaZapałek);
        
        komunikat("Palacz " + std::to_string(palaczID) + " pali fajkę.");
        std::this_thread::sleep_for(std::chrono::seconds(2)); // symulacja czasu trwania palenia fajki
    }
}

int main() {
    // Inicjalizacja semaforów
    sem_init(&ubijacze, 0, l);
    sem_init(&pudełkaZapałek, 0, m);

    // Tworzenie i uruchamianie wątków palaczy
    std::thread palacze[k];
    for (int i = 0; i < k; ++i) {
        palacze[i] = std::thread(palacz, i);
    }

    // Oczekiwanie na zakończenie wątków palaczy
    for (int i = 0; i < k; ++i) {
        palacze[i].join();
    }

    // Zwalnianie semaforów
    sem_destroy(&ubijacze);
    sem_destroy(&pudełkaZapałek);

    return 0;
}