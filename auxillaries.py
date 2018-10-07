"""
Created on Tue Aug 21 16:16:26 2018
Author: Nicholas Irvin

Purpose: These auxillary files can be copy and pasted into the DB program for special tasks
"""








"""
Profiler.
The profiler details how much time each function takes to run in order to 
evaluate bottlenecks in speed. Pay attention to the tottime column.
Note that JIT cannot be used concurrently with the profiler.
"""
# Paste this code into the functions module. You might also have to paste in 
# the parameters from "Specify Parameters" from the run file as well as 
# print(generate_LU_table(spectr, LL, UL, nbg, ArrSize, tol, conc))
import cProfile, pstats, io
def profile(fnc): # 
         
    def inner(*args, **kwargs):
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval
    return inner
# Paste "@profile" right before the definition of the function you want to analyze.
# Paste the parameters into the function's file, then put call to print the function 
# with those parameters.





# Get the average time for n uses of the code. 
def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped 
wrapped = wrapper(optimize_stack,LL, UL, resolution, spectrum, concentration, control, print_file)
# The arguments here should be the name of the code's main function followed by 
# that function's arguments.
n=1
print(timeit.timeit(wrapped, number=n)/n)







"""
JIT
The JIT decorator speeds up array computation. Note that JIT cannot be used 
concurrently with the profiler.
"""

from numba import jit
# Paste "@jit" right before "def particle_flux" within the functions module







"""
Temperature sweep

Sweep through temepratures while outputing just maximum efficiency and optimal bandgap.
"""






"""
Silicon bandgap (eV) as a function of tempaerature (K)
"""

def SiBG(T_cell):
    return round(1.17 - 4.73E-4*T_cell**2/(T_cell + 636), 4)
# Source: https://books.google.com/books?id=Yg3SBQAAQBAJ&pg=PA118&lpg=PA118&dq=updated+varshni+parameters+for+silicon&source=bl&ots=xWQTc6dsU0&sig=5LKn_S2IBewzf25j7tFrcnmqPu8&hl=en&sa=X&ved=2ahUKEwiy6LOS79bdAhXDIDQIHSHUAEYQ6AEwBnoECAAQAQ#v=onepage&q=updated%20varshni%20parameters%20for%20silicon&f=false
