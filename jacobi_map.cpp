//
// Created by marco on 9/2/17.
//

#include "jacobi_map.hpp"
#include "ff/map.hpp"

using namespace std;

namespace {

    typedef struct task_t {
        float *solutions;
        float *old_solutions;
        float error;
    } task_t;


    float **coefficients;

    float *terms;

    float tolerance;

    float error_computed;

    ulong max_iterations;
    ulong size;


    float *solution = nullptr;


    struct mapWorker : ff::ff_Map<task_t, task_t, float> {
        task_t *svc(task_t *input) {
            this->parallel_for(0, size, [&input](const ulong i) {
                float tmp = terms[i];
                for (ulong j = 0; j < size; ++j) {
                    if (i == j) continue;
                    tmp -= (input->old_solutions[j]) * (coefficients[i][j]);;
                }
                input->solutions[i] = tmp / (coefficients[i][i]);
            });

            auto reduce = [input](float &var, const ulong i) {
                var += abs(input->solutions[i] - input->old_solutions[i]);
            };
            auto reduce2 = [input](const ulong i, float &var) {
                var += abs(input->solutions[i] - input->old_solutions[i]);
            };

            this->parallel_reduce(input->error, 0.f, 0, size, reduce2, reduce);

            this->parallel_for(0, size, [&input](const ulong i) {
                input->old_solutions[i] = input->solutions[i];
            });

            ff_send_out(input);
            return GO_ON;

        }
    };

    ulong iteration = 0;


    struct receiver : ff::ff_monode_t<task_t> {
        task_t *svc(task_t *input) override {
            ff_send_out_to(input, 0);
            return GO_ON;
        }
    };


    struct generator : ff::ff_node_t<task_t> {
        task_t *svc(task_t *input) override {
            if (input == nullptr) {
                auto task = new task_t();
                task->solutions = new float[size]{(tolerance - tolerance)};
                task->old_solutions = new float[size]{(tolerance - tolerance)};
                return task;
            }
            ++iteration;
            input->error /= (float) size;
            if (input->error > tolerance && iteration < max_iterations) {
                input->error = 0.f;
                return input;
            }


            solution = input->solutions;
            delete (input->old_solutions);
            delete (input);
            return EOS;
        }

    };
}

#include <iostream>


float *jacobi_map(float **_coefficients, float *_terms, const ulong _size, const ulong _iterations,
                  const float _tolerance,
                  const ulong nworkers) {

    //setting up global data structure
    coefficients = _coefficients;
    terms = _terms;
    tolerance = _tolerance;
    max_iterations = _iterations;
    size = _size;

    std::vector<std::unique_ptr<ff::ff_node>> workers;
    for (ulong i = 0; i < nworkers; ++i) {
        workers.push_back(std::make_unique<mapWorker>());
    }

    ff::ff_Farm<task_t> farm(std::move(workers));
    generator myEmitter;
    receiver myReceiver;
    ff::ff_Pipe<task_t> pipe(myEmitter, farm, myReceiver);
    pipe.wrap_around();
    if (pipe.run_and_wait_end() < 0)
        std::cerr << "running pipe" << std::endl;
    std::cout << "iteration computed: " << iteration - 1 << " error: " << error_computed << std::endl;

    return solution;
}

