from subprocess import call
import os

os.environ["OMP_PLACES"] = "cores"
os.environ["OMP_PROC_BIND"] = "close"
os.environ["OMP_DYNAMIC"] = "false"
for n in range(6,36):
    call(["bin/sw_solve_big", str(n+1), str(10)], cwd="..")
