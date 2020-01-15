from subprocess import call
import os

#for n in range(6,36):
#    call(["bin/sw_solve_big", str(n+1), str(10)], cwd="..")

problem_size_arr=[1000,2000,2000,2000,4000,4000,8000]
n_threads_arr = [2,4,6,8,16,18,20,22,24,25,26,28,32,40,48,50,52,58,64,72]
while True:
    for problem_size in problem_size_arr:
        if problem_size==8000:
            n_threads_arr = [8,16,20,24,32,40,50,64,72]
        for n_threads in n_threads_arr:
            print("*"*64)
            call(["./swps3", "-j", "%d"%n_threads, "matrices/blosum50.mat", 
                "test/query%d.fa"%problem_size, "test/db%d.fa"%problem_size], 
                    cwd="../benchmark/swps/swps3-20120913")
