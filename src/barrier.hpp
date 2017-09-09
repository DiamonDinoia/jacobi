//
// Created by marco on 9/6/17.
//

#ifndef JACOBI_BARRIER_H
#define JACOBI_BARRIER_H

#include <atomic>


class spinning_barrier {
public:
    explicit spinning_barrier(unsigned long n) : n_(n), nwait_(0), step_(0) {}

    inline bool wait() {
        unsigned long step = step_.load();

        if (nwait_.fetch_add(1) == n_ - 1) {
            /* OK, last thread to come.  */
            nwait_.store(0); // XXX: maybe can use relaxed ordering here ??
            step_.fetch_add(1);
            return true;
        } else {
            /* Run in circles and scream like a little girl.  */
            while (step_.load() == step);
            return false;
        }
    }

protected:
    /* Number of synchronized threads. */
    const unsigned long n_;

    /* Number of threads currently spinning.  */
    std::atomic<unsigned long> nwait_;

    /* Number of barrier synchronizations completed so far,
     * it's OK to wrap.  */
    std::atomic<unsigned long> step_;
};

#endif //JACOBI_BARRIER_H
