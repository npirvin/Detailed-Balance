# Detailed-Balance
Repository for generalized detailed-balance code. This package will be combined with ASU's photovoltaic package.

Edit the parameters in run.py. Pay attention to the 'upper_limits' and 'lower_limits' as well as 'spectrum', 'concentration', and 'control'.
The main file run.py runs stack.py, which builds the output as well as controls the way the cells are connected - i.e. in indepenent, series, or parallel connections.
The file stack.py uses the file single_cell.py to calculate the power of each cell.
The file auxillaires.py has some extra functions that you can paste in. For example, the profile can detail which functions take the most computation time.
