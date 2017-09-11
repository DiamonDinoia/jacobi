//
// Created by marco on 9/6/17.
//

#ifndef JACOBI_BARRIER_H
#define JACOBI_BARRIER_H

#include <atomic>

/**
 * simple barrier used to have a lock-free synchronization between threads
 */
struct spinning_barrier {
    unsigned const count;
    std::atomic<unsigned> spaces;
    std::atomic<unsigned> generation;

    explicit spinning_barrier(unsigned count_) :
            count(count_), spaces(count_), generation(0) {}

    inline  __attribute__((always_inline))
    void wait() {
        unsigned const my_generation = generation;
        if (!--spaces) {
            spaces = count;
            ++generation;
        } else {
            while (generation == my_generation);
        }
    }
};

#endif //JACOBI_BARRIER_H
