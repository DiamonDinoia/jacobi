//
// Created by marco on 9/2/17.
//

#include "jacobi_map.hpp"
#include "ff/map.hpp"
#include <iostream>


using namespace std;

namespace {

    typedef struct task_t {
        float *solutions;
        float *old_solutions;
        float error;
    } task_t;


    float **coefficients = nullptr;

    float *terms = nullptr;

    float tolerance = 0.f;

    float errors = 0.f;

    ulong max_iterations = 0;
    ulong size = 0;


    float *solution = nullptr;
    ulong nworkers = 8;

    atomic_flag flag;

    /**
     * same as the serial implementation, the only difference is that each worker computes a subset of the soluton
     */
    struct mapWorker : ff::ff_Map<task_t, task_t, float> {


        task_t *svc(task_t *input) override {
            // computes the solutions
            auto reduce = [&](float &var, const ulong i) {
                var += abs(input->solutions[i] - input->old_solutions[i]);
            };
            auto reduce2 = [&](const ulong i, float &var) {
                var += abs(input->solutions[i] - input->old_solutions[i]);

            };
            this->parallel_for_static(0, size, 1, 0, [&input](const ulong i) {
                float tmp = terms[i];
#pragma simd
                for (ulong j = 0; j < size; ++j) {
                    tmp -= (input->old_solutions[j] * coefficients[i][j]);;
                }
                tmp += (input->old_solutions[i] * coefficients[i][i]);
                input->solutions[i] = tmp / (coefficients[i][i]);
            }, nworkers);

//
//            this->parallel_reduce_static(input->error, 0.f, 0, size, 1, 0, reduce2, reduce, nworkers);
//
//            this->parallel_for_static(0, size, 1, 0, [&input](const ulong i) {
//                input->old_solutions[i] = input->solutions[i];
//            }, nworkers);
            // calculate the error
            this->parallel_for_static(0, size, 1, 0, [&input](const ulong i, float error = 0.f) {
                error += abs(input->solutions[i] - input->old_solutions[i]);
                input->old_solutions[i] = input->solutions[i];
                while (!flag.test_and_set(std::memory_order_relaxed)) {}
                input->error += error;
                flag.clear(std::memory_order_relaxed);
            }, nworkers);
            ff_send_out(input);
            return GO_ON;

        }
    };


    /**
     * Just pass the results back to the emitter.
     * I think there is a way to not use this stage I did not find it.
     */
    struct receiver : ff::ff_monode_t<task_t> {
        task_t *svc(task_t *input) override {
            ff_send_out_to(input, 0);
            return GO_ON;
        }
    };


    ulong iteration = 0;

    auto start = Time::now();
    auto end_time = Time::now();


    /**
     * The emitter:
     *  First time initializes the data structures and send it to the workers
     *  Each iteration checks the error and the number of iterations
     *  Last time save the results, clean up the head and terminates the other stage
     */
    struct generator : ff::ff_node_t<task_t> {

        task_t *task = nullptr;

        int svc_init() override {
            task = new task_t();
            task->solutions = new float[size]{(tolerance - tolerance)};
            task->old_solutions = new float[size]{(tolerance - tolerance)};
            task->error = 0.f;
            start = Time::now();
            flag.clear();
        }

        task_t *svc(task_t *input) override {
            // first time initialization
            if (input == nullptr) return task;
            // check error and iterations
            ++iteration;
            errors = input->error / (float) size;
            if (errors > tolerance && iteration < max_iterations) {
                input->error = 0.f;
                return input;
            }
            // save the solution and terminates
            end_time = Time::now();
            solution = input->solutions;
            delete (input->old_solutions);
            delete (input);
            return EOS;
        }
    };
}


float *jacobi_map(float **_coefficients, float *_terms, const ulong _size,
                  const ulong _iterations, const float _tolerance, const ulong _nworkers) {

    //setting up global data structure
    coefficients = _coefficients;
    terms = _terms;
    tolerance = _tolerance;
    max_iterations = _iterations;
    size = _size;
    iteration = 0;
    nworkers = _nworkers;
    // creating the worker

    std::vector<std::unique_ptr<ff::ff_node>> workers;
    for (ulong i = 0; i < 1; ++i) {
        workers.push_back(ff::make_unique<mapWorker>());
    }

    // fastflow core code
    ff::ff_Farm<task_t> farm(std::move(workers));
    generator myEmitter;
    receiver myReceiver;
    ff::ff_Pipe<task_t> pipe(myEmitter, farm, myReceiver);
    pipe.wrap_around();

    if (pipe.run_and_wait_end() < 0)
        std::cerr << "map jacobi ERROR running pipe" << std::endl;
    std::cout << "map jacobi | iterations computed: " << iteration << " error: " << errors << std::endl;
    std::cout << "map jacobi | computation time: " << dsec(end_time - start).count() << std::endl;
    return solution;
}

