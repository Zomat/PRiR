#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <ncurses.h>
#include <string.h>
#include <stdbool.h>
#include <ctime>

#define MAX_WAIT_TIME 3000000
#define EAT_TIME 3

enum states
{
    THINKING, EATING, HUNGRY
};

struct DiningPhils
{
    pthread_cond_t cond[10];
    pthread_mutex_t lock;
    int count[10];
    bool fork[10];
    long total_wait_time[10];
    enum states state[10];
};

struct DiningPhils dp;

pthread_t visualization;

long wait_temp[10];

int phils_number = 5;

const char *getStateName(enum states s)
{
    switch (s)
    {
    case THINKING:
        return "MYSLI    ";
    case EATING:
        return "JE       ";
    case HUNGRY:
        return "GLODNY   ";
    default:
        return "???";
    }
}

const char *getForkStateName(bool s)
{
    if (s)
    {
        return "WOLNY ";
    }

    return "ZAJETY";
}

// Pobranie aktualnego czasu zegara
long get_time()
{
    struct timespec ts;

    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
    {
        return (long)(ts.tv_sec * 1000000000 + ts.tv_nsec);
    }

    return 0;
}

/*
 * Pozwol innemu filozofowi zjesc, jezeli jest glodny
 */
int give_way(int i)
{
    long now = get_time();
    long n1 = 0, n2 = 0, cur = wait_temp[i] + now;
    if (dp.state[(i + 1) % phils_number] == HUNGRY)
    {
        n1 = wait_temp[(i + 1) % phils_number] + now;
        if (n1 > cur && n1 > MAX_WAIT_TIME)
            return 0;
    }

    if (dp.state[(i + phils_number - 1) % phils_number] == HUNGRY)
    {
        n2 = wait_temp[(i + phils_number - 1) % phils_number] + now;
        if (n2 > cur && n2 > MAX_WAIT_TIME)
            return 0;
    }
    return 1;
}

void check_forks(int phil_id)
{

    if (dp.state[phil_id] == HUNGRY &&
        dp.state[(phil_id + 1) % phils_number] != EATING &&
        dp.state[(phil_id + phils_number - 1) % phils_number] != EATING &&
        give_way(phil_id))
    {
        // Mozna jesc
        dp.state[phil_id] = EATING;

        // Paleczki zajete
        dp.fork[phil_id] = false;
        dp.fork[(phil_id + 1) % phils_number] = false;

        // Zapisz czas oczekiwania na jedzenie
        wait_temp[phil_id] += get_time();
        dp.total_wait_time[phil_id] += wait_temp[phil_id];

        // Wyslij sygnal ze mozna jesc
        pthread_cond_signal(&dp.cond[phil_id]);
    }
}

void return_forks(int phil_id)
{
    pthread_mutex_lock(&dp.lock);
    dp.state[phil_id] = THINKING;
    dp.count[phil_id]++;

    dp.fork[phil_id] = true;
    dp.fork[(phil_id + 1) % phils_number] = true;

    // Poinformuj filozofa po lewej i prawej, żeby sprawdził widelce
    check_forks((phil_id + phils_number - 1) % phils_number);
    check_forks((phil_id + 1) % phils_number);

    pthread_mutex_unlock(&dp.lock);
}

// Filozof jest glodny
// Sprawdza wolne widelce, jezeli nie sa wolne czeka na sygnal
void grab_forks(int phil_id)
{
    pthread_mutex_lock(&dp.lock);
    // Sekcja krytyczna
    dp.state[phil_id] = HUNGRY;

    // zapisanie startu czasu uczekiwania
    wait_temp[phil_id] = (-1) * get_time();

    check_forks(phil_id);

    if (dp.state[phil_id] != EATING)
    {
        // Odblokuj locka i czekaj na sygnal ze wzial widelce
        pthread_cond_wait(&dp.cond[phil_id], &dp.lock);
    }

    pthread_mutex_unlock(&dp.lock);
}

void *philosopher(void *param)
{
    int phil_id = *((int *)param);

    while (1)
    {
        // Losuj czas jedzenie i myslenia filozofa
        srand(time(NULL) + phil_id);
        int think_time = (rand() % EAT_TIME) + 1;
        int eat_time = (rand() % EAT_TIME) + 1;

        dp.state[phil_id] = THINKING;
        sleep(think_time);

        // Sprobuj wziazc widelce
        grab_forks(phil_id);
        // Jezeli ten kod sie wykonuje, to znaczy udalo sie wziazc 2 widelce
        // Jedz...
        sleep(eat_time);

        // Zwroc widelce i zwolnij zasoby
        return_forks(phil_id);
    }
}

void *visualize(void*)
{
    initscr();

    while (1)
    {
        clear();

        long wait_time = 0;
        char status[20];

        printw("|---------------------------------------------------|\n");
        printw("|------------- UCZTUJACY FILOZOFOWE ----------------|\n");
        printw("|---------------------------------------------------|\n\n");

        for (int i = 0; i < phils_number; i++)
        {
            wait_time += (long)dp.total_wait_time[i] / 1000000;
            printw("Stan Filozofa %d: %s | ", i, getStateName(dp.state[i]));
            printw("Stan Zasobu %d: %s  |\n", i, getForkStateName(dp.fork[i]));
            printw("Jadl %d razy\n", dp.count[i]);
            printw("----------------------------------------------------|\n");
        }

        refresh();
        usleep(25000);
    }

    endwin();
}

int main(int argc, char *argv[])
{
    int i;

    if (argc > 1 && atoi(argv[1]) > 4 && atoi(argv[1]) < 11)
    {
        phils_number = atoi(argv[1]);
    }

    // Inicjalizacja mutexa
    pthread_mutex_init(&dp.lock, nullptr);

    pthread_t tid[phils_number];
    int id[phils_number];
    for (i = 0; i < phils_number; i++)
    {
        id[i] = i;
        // Inicjalizacja zmiennych warunkowych
        pthread_cond_init(&dp.cond[i], nullptr);

        // Inicjalizacja stanow filozofa (domyslnie mysli)
        dp.state[i] = THINKING;

        // Inicjalizacja licznika zjedzonych posilkow
        dp.count[i] = 0;

        // Inicjalizacja watkow
        pthread_create(&tid[i], nullptr, philosopher, &id[i]);

        // Domyslnie widelce sa wolne
        dp.fork[i] = true;
    }

    // Inicjalizuj watek wizualizacji
    pthread_create(&visualization, nullptr, visualize, nullptr);

    // main moze sie zakonczyc tylko wtedy,
    // kiedy zakoncza sie wszystkie watki podane ponizej
    for (i = 0; i < phils_number; i++)
        pthread_join(tid[i], nullptr);

    return 0;
}
