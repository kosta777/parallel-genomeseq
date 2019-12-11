from subprocess import call

for n in range(32):
    print("npiece =" + str(n+1))
    call(["bin/sw_solve_big", str(n+1)], cwd="..")