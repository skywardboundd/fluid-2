#include <iostream>
#include <vector>
#include <thread>
#include <queue>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <atomic>

class ThreadPool {
public:
    explicit ThreadPool(size_t threadCount) : stop(false) {
        for (size_t i = 0; i < threadCount; ++i) {
            workers.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(this->queueMutex);
                        this->condition.wait(lock, [this]() {
                            return this->stop || !this->tasks.empty();
                        });

                        if (this->stop && this->tasks.empty())
                            return; 

                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }

                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers)
            worker.join();
    }

    template <class F, class... Args>
    void enqueue(F&& f, Args&&... args) {
        {
            std::unique_lock<std::mutex> lock(queueMutex);

            tasks.emplace([f = std::forward<F>(f), ...args = std::forward<Args>(args)]() {
                f(args...);
            });
        }
        condition.notify_one();
    }

private:
    std::vector<std::thread> workers; 
    std::queue<std::function<void()>> tasks; 

    std::mutex queueMutex; 
    std::condition_variable condition; 
    std::atomic<bool> stop; 
};

// int main() {
//     ThreadPool pool(4);

//     for (int i = 0; i < 10; ++i) {
//         pool.enqueue([i]() {
//             std::cout << "Task " << i << " is running on thread " << std::this_thread::get_id() << std::endl;
//             std::this_thread::sleep_for(std::chrono::milliseconds(500));
//         });
//     }

//     std::this_thread::sleep_for(std::chrono::seconds(3)); 
//     return 0;
// }
