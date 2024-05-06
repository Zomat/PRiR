#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>

using namespace std;

#define THREADS_NUM 5

class monitor {
private:
    int reader_count = 0;
    int writer_count = 0;

    int waiting_readers = 0;
    int waiting_writers = 0;

    int data = 0;
    int reads_count = 0;

    condition_variable canRead;
    condition_variable canWrite;
    
    mutex coutLock;
    mutex condLock;

public:
    void message(string const &msg) {
        unique_lock<mutex> lock(coutLock);
        cout << msg << endl;
    }

    void beginRead(int i) {
        unique_lock<mutex> lock(condLock);

        if (writer_count == 1 || waiting_writers > 0) {
            waiting_readers++;
            canRead.wait(lock, [&](){ return writer_count == 0; });
            waiting_readers--;
        }

        reader_count++;
        reads_count++;
        message("Reader " + to_string(i) + " is reading: " + to_string(data));
        canRead.notify_all();
    }

    void endRead(int i) {
        unique_lock<mutex> lock(condLock);

        if (--reader_count == 0)
            canWrite.notify_one();
    }

    void beginWrite(int i) {
        unique_lock<mutex> lock(condLock);
        canWrite.wait(lock, [&](){ return reads_count >= 3 || data == 0; });
        
        if (writer_count == 1) {
            ++waiting_writers;
            canWrite.wait(lock, [&](){ return reader_count == 0 && writer_count == 0; });
            --waiting_writers;
        }
        writer_count = 1;
        data = rand() % 100 + 1;
        reads_count = 0;
        message("----------------------------------------------------------");
        message("Writer " + to_string(i) + " is writing: " + to_string(data));
        message("----------------------------------------------------------");
    }

    void endWrite(int i) {
        unique_lock<mutex> lock(condLock);
        writer_count = 0;

        if (waiting_readers > 0)
            canRead.notify_all();
        else
            canWrite.notify_one();
    }
};

monitor M;

void reader(int id) {
    while (1) {
        this_thread::sleep_for(chrono::seconds(1));
        M.beginRead(id);
        M.endRead(id);
    }
}

void writer(int id) {
    while(1) {
        M.beginWrite(id);
        M.endWrite(id);
        this_thread::sleep_for(chrono::seconds(3));
    }
}

int main() {
    vector<thread> readers(THREADS_NUM);
    vector<thread> writers(THREADS_NUM);

    for (int i = 0; i < THREADS_NUM; i++) {
        readers[i] = thread(reader, i);
        writers[i] = thread(writer, i);
    }

    for (int i = 0; i < 5; i++) {
        readers[i].join();
        writers[i].join();
    }

    return 0;
}