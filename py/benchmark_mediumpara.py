from subprocess import call

for n in range(32):
    call(["bin/sw_solve_big", str(n+1)], cwd="..")