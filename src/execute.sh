#!/bin/bash

directory="./data"

rm -r -f $directory
mkdir -p $directory

cd $directory

iterations=50
workers=64
size=1024
tolerance=-1.f

algorithms=(thread)
#algorithms=(sequential thread fastflow omp)

sizes=(64 512 1024 2048 4096 16384)


echo "testing algorithms with different sizes"
for s in ${sizes[@]}; do
    for method in ${algorithms[@]}; do
        name=$method"_size.csv"
        ../main $method -t$tolerance -i$iterations -w$workers -s$s -p$name
    done
done

iter=(1 16 64 256 1024 2048 4096)

echo "testing algorithms with different iterations"
for i in ${iter[@]}; do
    for method in ${algorithms[@]}; do
        name=$method"_iterations.csv"
        ../main $method -t$tolerance -i$i -w$workers -s$size -p$name
    done
done

#unset algorithms
#algorithms=(thread fastflow omp)

echo "testing parallel algorithms with different workers"
for w in {1..128..1}; do
    for method in ${algorithms[@]}; do
        name=$method"_workers.csv"
        ../main $method -t$tolerance -i$iterations -w$w -s$size -p$name
    done
done

size=16384

echo "testing parallel algorithms with different workers and large size"
for w in {1..128..1}; do
    for method in ${algorithms[@]}; do
        name=$method"_workersxxl.csv"
        ../main $method -t$tolerance -i$iterations -w$w -s$size -p$name
    done
done
