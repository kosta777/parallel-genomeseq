from subprocess import call
import os

os.environ["OMP_PLACES"] = "cores"
os.environ["OMP_PROC_BIND"] = "close"
for n in range(36):
    call(["bin/sw_solve_big", str(n+1), str(10)], cwd="..")