#!/bin/bash

#SBATCH --job-name=hello-cuda
#SBATCH --output=hello-cuda-%j.out
#SBATCH --error=hello-cuda-%j.err

#SBATCH --partition training
#SBATCH --gres=gpu
#SBATCH --nodes 1
#SBATCH --time=01:00:00

#SBATCH --ntasks=1

module load cuda-11.2.1
module load gcc-6.5.0

/usr/local/cuda/bin/nvcc ../exp_force_main.cu -o ../output/ExForce

# for blocks in 540 1024 2048
# do
#     for stream_count in 2 4 8 
#     do
#         srun -N 1 ../output/ExForce ../new_fb_full.txt $blocks 1024 $stream_count 1 > stopwatch.txt
#     done 
# done

for threads in 256 512 1024
do
    for blocks in 32 64 128
    do
        for stream_count in 2 4 8 
        do
            srun -N 1 ../output/ExForce ../test_edgelist.txt $blocks $threads $stream_count 1 >> stopwatch.txt
        done 
    done
done
#../output/ExForce ../../test_edgelist.txt 1
#../output/ExForce ../../mini_graph.txt 1
#../output/ExForce ../new_fb_full.txt 540 1024 2 1


