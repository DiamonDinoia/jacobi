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

    struct mapWorker : ff::ff_Map<task_t, task_t, float> {
        task_t *svc(task_t *input) override {
            float error = 0.f;
            this->parallel_for(0, size, [&input, &error](const ulong i) {
                float tmp = terms[i];
#pragma ivdep
                for (ulong j = 0; j < size; ++j) {
                    if (i == j) continue;
                    tmp -= (input->old_solutions[j]) * (coefficients[i][j]);;
                }
                input->solutions[i] = tmp / (coefficients[i][i]);
                error += abs(input->solutions[i] - input->old_solutions[i]);
                input->old_solutions[i] = input->solutions[i];
            }, nworkers);
            input->error += error;

            ff_send_out(input);
            return GO_ON;

        }
    };



    struct receiver : ff::ff_monode_t<task_t> {
        task_t *svc(task_t *input) override {
            ff_send_out_to(input, 0);
            return GO_ON;
        }
    };


    ulong iteration = 0;

    auto start = Time::now();
    auto end_time = Time::now();



    struct generator : ff::ff_node_t<task_t> {
        task_t *svc(task_t *input) override {
            if (input == nullptr) {
                auto task = new task_t();
                task->solutions = new float[size]{(tolerance - tolerance)};
                task->old_solutions = new float[size]{(tolerance - tolerance)};
                task->error = 0.f;
                start = Time::now();
                flag.clear();
                return task;
            }

            ++iteration;
            errors = input->error / (float) size;
            if (errors > tolerance && iteration < max_iterations) {
                input->error = 0.f;
                return input;
            }
            end_time = Time::now();
            solution = input->solutions;
            delete (input->old_solutions);
            delete (input);
            return EOS;
        }

    };
}



float *jacobi_map(float **_coefficients, float *_terms, const ulong _size, const ulong _iterations,
                  const float _tolerance,
                  const ulong _nworkers) {

    //setting up global data structure
    coefficients = _coefficients;
    terms = _terms;
    tolerance = _tolerance;
    max_iterations = _iterations;
    size = _size;
    iteration = 0;
    nworkers = _nworkers;

    std::vector<std::unique_ptr<ff::ff_node>> workers;
    for (ulong i = 0; i < 1; ++i) {
        workers.push_back(ff::make_unique<mapWorker>());
    }

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

