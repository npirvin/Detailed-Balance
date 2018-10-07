""" Created on Fri Aug 24 11:47:20 2018
Author: Nicholas Irvin

Outline: It is assumed that the cell's current is only generated and recombinated
for energies above the cell's band-gap, E1, and below the next cell's 
band gap, E2. The single_cell module calculates the power of the cell
for light between energies E1 and E2. This module stack.py samples this calculation
for various values of E1 and E2 for each cell. Finally, the optimal combination
of bandgaps can be determined. """

# Import libraries.
import numpy as np
import time
from scipy.optimize import minimize_scalar
from decimal import Decimal
import matplotlib.pyplot as plt
import single_cell



# Set constants.
f = 1/46200 # (sr)  (The solid angle of the sun from Earth)/Pi
Stefan_Boltzmann = 5.67051e-8 # (J/K)  Stefan-Boltzmann coefficient in SI standard units
T_sun = 5778 # (K)  Approximate temperature of the sun







def independent_cell_power(energies, lower_limits, upper_limits, nbg, resolution, T_cell, spectrum, concentration): 
    """ Purpose: Generate look-up table of maximum power 'power' and voltage 
    at max power 'vmp' at every combination of energies [E1, E2] specified by 
    the lower_limits and upper_limits vectors. The function stack can then use this array to combine
    the power generate from each cell.
    
    The outputs are a vector then two matrices: energies (sampled energy points
    of E1 and E2, vmp (voltage of max power), and power. Each array element is
    for a different point of [E1, E2].
    Used by the function optimize_stack. """
    
    array_size = len(energies)
    #  Initialize
    power = np.zeros((array_size, array_size))
    vmp = np.zeros((array_size, array_size))
    
    if spectrum == 0:
        lower_limits = lower_limits + [10]
        upper_limits = upper_limits + [10]
    else: 
        lower_limits = lower_limits + [4.428]
        upper_limits = upper_limits + [4.428]
    # Iterate through points of [E1, E2] while calculating the maximum power
    # and voltage at maximum power
    for ia in range(array_size-1):
        E1 = energies[ia]
        for ja in range(array_size):  # E2 goes to infinity but E1 doesn't.
            E2 = energies[ja]
            
            # Check if need calculation.
            skip = 0
            for k in range(nbg):
                if lower_limits[k] <= E1 <= upper_limits[k] and lower_limits[k+1] <= E2 <= upper_limits[k+1]:
                    skip = 1
            if skip == 1:
                [vmp[ia, ja], power[ia, ja]]  =  single_cell.max_power_pt(E1, E2, T_cell, spectrum, concentration)
                # Calculate maximum power and voltage at maximum power (vmp)
    return(vmp, power)






def series_stack_power( bandgaps, nbg, T_cell, spectrum, concentration):
    """ Purpose: Calculate total power for a stack of cell wired in series for a
    given set of bandgaps.
    
    Used by the function optimize_stack. """ 
    

    Jsc_array = [single_cell.find_Jsc(bandgaps[j], bandgaps[j+1], T_cell, spectrum, concentration) for j in range(nbg)]
    Jsc_low = min(Jsc_array)    
    # Create total power function that is to be maximizesd by varying voltage.    
    def power(current):
        volt = 0  # Initialize
        for m in range(nbg):  
            # Series connection, so add voltages and keep current equivalent.
            volt = volt + single_cell.find_voltage(
                    current, bandgaps[m], bandgaps[m+1], T_cell, spectrum, concentration, Jsc_array[m])  
            # This line uses the series assumption.
        return -current*volt 
    
    res = minimize_scalar(power, bounds=[0, Jsc_low], method = 'Bounded')
    return(abs(res.fun))  
    # Output max power.

    
    



def transmittance_series_stack_power( bandgaps, nbg, T_cell, spectrum, concentration):
    """ Purpose: Calculate total power for a stack of cell wired in series for a
    given set of bandgaps. Here we allow absortiance of the top cell to vary.
    
    Used by the function optimize_stack. """ 
    
    if spectrum == 0:
        infinity = 10
    else:
        infinity = 4.428
    
    power_max = 0
    for i in range(1001):
        a = i/1000
        
        # Find the stack's maximum current as the lowest Jsc.
        Jsc_array = [(single_cell.find_Jsc(bandgaps[0], bandgaps[1], T_cell, spectrum, concentration) + (1-a)*single_cell.find_Jsc(bandgaps[1], infinity, T_cell, spectrum, concentration)), 
        a*single_cell.find_Jsc(bandgaps[1], infinity, T_cell, spectrum, concentration)]
        Jsc_low = min(Jsc_array)
        
        # Create total power function that is to be maximized by varying voltage.    
        def power(current):
            volt = 0  # Initialize
            for m in range(nbg):  
                # Series connection, so add voltages and keep current equivalent.
                volt = volt + single_cell.find_voltage(
                        current, bandgaps[m], bandgaps[m+1], T_cell, spectrum, concentration, Jsc_array[m])  # This line uses the series assumption.
            return -current*volt 
        
        res = minimize_scalar(power, bounds=[0, Jsc_low], method = 'Bounded')
        power_max = max(power_max, abs(res.fun))
    return(power_max)  
    
    
    
    
    
    
def parallel_stack_power(bandgaps, nbg, T_cell, spectrum, concentration):
    """ Purpose: Calculate total power for a stack of cell wired in parallel for a
    given set of bandgaps.
    
    Used by the function optimize_stack. """ 
    
    # Construct each Jsc
    Jsc_array = [single_cell.find_Jsc(T_cell, bandgaps[j], bandgaps[j+1], T_cell, spectrum, concentration) for j in range(nbg)]

    Voc_array = [single_cell.find_voltage(T_cell, 0, bandgaps[j], bandgaps[j+1], T_cell, spectrum, concentration, Jsc_array[j]) for j in range(nbg)]

    power_max = 0  # Initialize
    
    
    # The voltage of each cell is primarily determined by E1 for each cell,
    # so the voltages are inherently mismatched. In order to match the voltages,
    # allow for each layer to be a number of cells connected in series.
    number_of_cells = np.ones(nbg, dtype='int')  
    number_of_cells[0] = 60  # Keep the bottom layer with 60 cells while varying 
                             # the number for the other layers.
    # For more info, see "Fundamental Analysis of Triple Layer Area Decoupled 
    # Photovoltaic Modules" by Saroj Pyakurel - page 11.
    
    
    for m in range(60**nbg): 
        
        # Create total power function that is to be maximized by varying voltage. 
        def power(volt):
            current = np.zeros(nbg)
            for k in range(nbg):
                current[k] = single_cell.find_current(volt/number_of_cells[k], bandgaps[k], bandgaps[k+1], T_cell, spectrum, concentration, Jsc_array[k])/number_of_cells[k]
            return -sum(current)*volt 
        
        res = minimize_scalar(power, bounds=[0, min(np.multiply(number_of_cells, Voc_array))], method='Bounded')
        power_max = max(power_max, abs(res.fun))
    
    
        # Move to the next number_of_cells configuration.
        if number_of_cells[nbg-1] == 60: 
             break  # Quit if it is the last configuration.
        layer = 1   # Otherwise, start at the second layer. 
                    # The first layer stays at 60 cells.
        while number_of_cells[layer] == 60 or (
                all(number_of_cells[i] <= number_of_cells[i+1] for i in range(len(number_of_cells)-1))):
            # number_of_cells of a layer is limited by 60 and the number of the next layer
            number_of_cells[layer] = 1 # return the bandgap to its lower limit...
            layer += 1  # and move to the next layer.
        number_of_cells[layer] += 1  
        # Increment the number of cells to find power at the optimal configuration.
        
    return abs(power_max)
    # Output max power.






def intermediate_band_power(bandgaps, nbg, T_cell, spectrum, concentration):    
    """ Purpose: Calculate total power for a intermediate-band cell. This is
    modeled by a bottom and middle cell in series followed by connection to the 
    top cell in parallel.
    
    Used by the function optimize_stack. """ 
    
    # The bottom cells' bandgaps must add up to the top cell bandgap.
    if bandgaps[0] + bandgaps[1] != bandgaps[2]:
        return 0
    
    Jsc_array = [single_cell.find_Jsc(bandgaps[j], bandgaps[j+1], T_cell, spectrum, concentration) for j in range(nbg)]
    Jsc_eff = min(Jsc_array[0], Jsc_array[1]) 
    # The current through the two in series is limited by the lowest Jsc.
        
    # Create total power function that is to be maximized by varying current 
    # through the intermediate band. 
    def power(current):
        current1 = current
        current2 = current  # Series condition
        voltage1 = single_cell.find_voltage(
                current, bandgaps[0], bandgaps[0+1], T_cell, spectrum, concentration, Jsc_array[0])
        voltage2 = single_cell.find_voltage(
                current, bandgaps[1], bandgaps[1+1], T_cell, spectrum, concentration, Jsc_array[1])
        voltage3 = voltage1 + voltage2  # Parallel condition
#        Voc_3 = single_cell.find_voltage(0, bandgaps[2], bandgaps[2+1], spectrum, concentration)
#        # Limit voltage by this Voc.
#        if voltage3 > Voc_3:
#            return 0
        current3 = single_cell.find_current(
                voltage3, bandgaps[2], bandgaps[2+1], T_cell, spectrum, concentration, Jsc_array[2]) 
        return -(current1*voltage1 + current2*voltage2 + current3*voltage3)
    
    return(abs(minimize_scalar(power, bounds=[0, Jsc_eff], method='Bounded').fun))  
    # Output max power.






"""
Under Construction
"""
def combination_stack(bandgaps, nbg, T_cell, spectrum, concentration, combination):
    power = 0
    for connection in combination:
        if connection == "Parallel":
            power = power +  intermediate_band_power(bandgaps, nbg, T_cell, spectrum, concentration)
        if connection == "Series":
            power = power +  intermediate_band_power(bandgaps, nbg, T_cell, spectrum, concentration)
    return power






    
def optimize_stack(lower_limits, upper_limits, resolution, T_cell, spectrum, concentration, control):
    """ Purpose: Maximize the total power such that, for each cell, E1 < E2, 
    and the E2 of a cell must be the E1 of the cell above it. 
    
    Used by module run. """
    
    start = time.time() # Start timing calculations
    
    nbg = len(upper_limits) # Number of bandgaps
    # Calculate incident sunlight power to be used in the efficiency ratio.
    if spectrum == 0: # Blackbody
        power_in = concentration*f*Stefan_Boltzmann*T_sun**4
        infinity = 10.  # Solar radiation is negligible above 10
    elif spectrum == 1:  # AM1.5G 
       power_in = concentration*1000.37
       infinity = 4.428  # AM1.5G and AM1.5D only reach energies of 4.428
    elif spectrum == 2:  # AM1.5D                
        power_in = concentration*900.14
        infinity = 4.428

    energies = np.arange(lower_limits[0], upper_limits[-1]+2*resolution, resolution)  
    # +2*resolution because arange doesn't include the endpoint, 
    # and we also want room for infinity.
    energies[-2] = upper_limits[-1]
    energies[-1] = infinity  # Set the last energy to cover the whole energy range
    
    # Create the list that stores the bandgap values that will be iteratively sampled.  
    bandgaps = lower_limits + [infinity]

    # Initialize Results.
    optimal_bandgaps = np.ones(nbg) # Optimum band gaps, 
    optimal_efficiency = 0.0; # Max efficiency
    file_data=[]

        
    # Calculates the voltage at max power and power generated for a single cell 
    # at each set [E1, E2]. E1 and E2 are sampled from vector energies
    if control == 0:  # Independent connection
        [vmp, power] = independent_cell_power(
            energies, lower_limits, upper_limits, nbg, resolution, T_cell, spectrum, concentration)
        # vmp is voltage at max power
        
                
    # For each bandgap configuration,
    for j in range(len(energies)**nbg):
        eff_all = 0.0  # Calculate total efficiency.
        
        if control == 0:  # For independent connection
            for i in range(nbg):  # For each cell,
                # find the elements in the power array that correspond to... 
                index_bottom = np.searchsorted(energies,bandgaps[i])
                # index_bottom = np.where(np.isclose(energies,bandgaps[i]))[0][0]
                # the cell's lower ..
                index_top = min(np.searchsorted(energies,bandgaps[i+1]), len(energies)-1)
                # and upper energy limits...
                eff_all  =  round(eff_all + 1e2*power[index_bottom,index_top]/power_in, 3) 
                # to find output power between those limits.
        else:
            if control == 1:  # For series connection                 
                power = transmittance_series_stack_power(bandgaps, nbg, T_cell, spectrum, concentration)
                # series_stack_power already gives total power of the stack
            if control == 2:  # For parallel connection                 
                power = parallel_stack_power(bandgaps, nbg, T_cell, spectrum, concentration)
                # series_stack_power already gives total power of the stack
            if control == 3:  # For parallel connection                 
                power = intermediate_band_power(bandgaps, nbg, T_cell, spectrum, concentration)
                # series_stack_power already gives total power of the stack
            eff_all = round(100*power/power_in, 3)
        
        # Save results
        file_data.append([bandgaps[0:nbg], eff_all])
        # Append vmp here to store voltage at max power            
        # Update max efficiency and optimal bandgaps
        if (eff_all >= optimal_efficiency): 
            optimal_efficiency = eff_all;
            for i in range(nbg):
                optimal_bandgaps[i] = bandgaps[i];
             
        
        # Move to the next bandgap configuration.
        if (nbg == 1 and bandgaps[nbg-1]>=upper_limits[nbg-1]) or (
                bandgaps[:-1] == upper_limits):
             break  # Quit if it is the last configuration.
        cell = 0 #  Otherwise, start at the bottom cell.
        while (bandgaps[cell] >= upper_limits[cell]) or (bandgaps[cell] >= bandgaps[cell+1]):
            # For each cell above its upper limits,
            bandgaps[cell] = lower_limits[cell] # return the bandgap to its lower limit...
            cell += 1  # and move to the next cell.
        bandgaps[cell] = round( bandgaps[cell] + resolution, 5 )  # Increment the bandgap...
                                              # And round to avoid numerical drifting
                         
                            
    # Write the file
    file = open('Data.txt', 'w')
    file.write('The maximum efficiency is ')
    file.write(str(optimal_efficiency))
    file.write('% with bandgaps of\n')
    file.write(str(optimal_bandgaps[0:nbg]))
    file.write(' ev. \n\n')
    
    file.write('Temperature (C) \t Bandgaps (eV)')
    for i in range(nbg):
        file.write('\t')
    file.write('\t')
    file.write('Efficiency (%)\n')
    
    for j in range(len(file_data)):
        file.write(str(T_cell-273.15))
        file.write('\t\t\t')
        for i in range(nbg):
            file.write(str(Decimal(str(file_data[j][0][i]))))  # Save bandgaps.
            file.write('\t\t')
        file.write(str(file_data[j][1]))  # Save efficiency.
        file.write('\n')
    file.close()
        
    
    # Plot bottom bandgap vs efficiency.
    plt.xlabel('Bottom-Cell Bandgap (eV)')
    plt.ylabel('Maximum Efficiency (%)')
    plt.legend()
    plt.plot([file_data[i][0][0] for i in range(int((upper_limits[0]-lower_limits[0])/resolution))], [file_data[i][1] for i in range(int((upper_limits[0]-lower_limits[0])/resolution))])
    
    
    # Prepare to print main results. 
    # One could also append vmp (voltage at max power) for the max efficiency.
    readout = ['\nThe maximum efficiency is '] #\n makes a new line.
    readout.append(str(optimal_efficiency))
    readout.append('% with bandgaps of\n')
    readout.append(optimal_bandgaps[:nbg])
    readout.append(' ev. \n\n')
    readout.append('Efficiencies at other band gaps have been written onto' 
               ' file Data.txt. \n')
    
    
    stop = time.time() # Stop the timer for the calculations.
    time_elapsed = stop - start
    readout.append('The calculations took ')
    readout.append(round(time_elapsed, 3))
    readout.append(' seconds.')

    return(readout)
    
