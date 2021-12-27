#ifndef SAFE_QUEUE_H
#define SAFE_QUEUE_H

#include <mutex>
#include <queue>

template<typename T>
class SafeQueue {
private:
    std::queue<T> queue;
    std::mutex mutex;

public:
    SafeQueue() {}

    bool empty() {
        std::unique_lock<std::mutex> lock(mutex);
        return queue.empty();
    }

    int size() {
        std::unique_lock<std::mutex> lock(mutex);
        return queue.size();
    }

    void enqueue(T& t) {
        std::unique_lock<std::mutex> lock(mutex);
        queue.push(t);
    }

    bool dequeue(T& t) {
        std::unique_lock<std::mutex> lock(mutex);

        if (queue.empty()) {
            return false;
        }

        t = std::move(queue.front());
        queue.pop();

        return true;
    }
};

#endif