#!/bin/bash

# Loop through all directories matching the pattern example for HfNiSn IIS data
for dir in */; do
    # Remove trailing slash
    dir=${dir%/}
    
    echo "Processing directory: $dir"
    
    # Extract method and temperature from directory name
    if [[ $dir =~ ^(IIS|matthiessen|screenedPOP)_xmgrace_7methods_([0-9]+K|Temperature_Comparison)$ ]]; then
        method="${BASH_REMATCH[1]}"
        temp="${BASH_REMATCH[2]}"
        
        # Skip the Temperature_Comparison directories if they exist
        if [[ $temp == "Temperature_Comparison" ]]; then
            echo "Skipping comparison directory: $dir"
            continue
        fi
        
        # Construct the filename based on method
        case $method in
            IIS)
                filename="HfNiSn_IIS_${temp}_7methods_overall.dat"
                ;;
            matthiessen)
                filename="HfNiSn_matthiessen_${temp}_7methods_overall.dat"
                ;;
            screenedPOP)
                filename="HfNiSn_screenedPOP_${temp}_7methods_overall.dat"
                ;;
        esac
        
        # Check if the directory exists and file exists
        if [ -d "$dir" ] && [ -f "$dir/$filename" ]; then
            cd "$dir" || continue
            
            echo "  Extracting data from $filename"
            
            # Extract columns to tau files (skip comment lines)
            for i in {1..7}; do
                awk -v col=$((i+1)) '!/^#/ && NF >= (col) {print $1, $col}' "$filename" > "tau_M${i}.dat"
            done
            
            echo "  Done processing $dir"
            cd ..
        else
            echo "  Warning: File $dir/$filename not found"
        fi
    else
        echo "  Skipping: $dir (doesn't match expected pattern)"
    fi
    
    echo "---"
done

echo ""
echo "========================================="
echo "Processing total scattering with interpolation..."
echo "========================================="

# Define interpolation parameters
E_MIN=0.0
E_MAX=0.5
N_STEPS=500

# Get all unique temperatures from the directories
declare -A temperatures

for dir in */; do
    dir=${dir%/}
    if [[ $dir =~ ^(IIS|matthiessen|screenedPOP)_xmgrace_7methods_([0-9]+K)$ ]]; then
        temp="${BASH_REMATCH[2]}"
        temperatures[$temp]=1
    fi
done

# Process each temperature
for temp in "${!temperatures[@]}"; do
    echo ""
    echo "Processing temperature: $temp"
    
    # Create the total_scattering folder
    total_dir="total_scattering_${temp}"
    mkdir -p "$total_dir"
    
    # Define the paths to tau files for each method
    iis_dir="IIS_xmgrace_7methods_${temp}"
    matth_dir="matthiessen_xmgrace_7methods_${temp}"
    pop_dir="screenedPOP_xmgrace_7methods_${temp}"
    
    # Check if all required directories exist
    if [ ! -d "$iis_dir" ] || [ ! -d "$matth_dir" ] || [ ! -d "$pop_dir" ]; then
        echo "  Warning: Missing directories for temperature $temp"
        continue
    fi
    
    # Process each method (M1 to M7)
    for m in {1..7}; do
        iis_file="$iis_dir/tau_M${m}.dat"
        matth_file="$matth_dir/tau_M${m}.dat"
        pop_file="$pop_dir/tau_M${m}.dat"
        output_file="$total_dir/tau_total_M${m}.dat"
        interpolated_dir="$total_dir/interpolated"
        mkdir -p "$interpolated_dir"
        
        iis_interp="$interpolated_dir/tau_IIS_M${m}.dat"
        matth_interp="$interpolated_dir/tau_matth_M${m}.dat"
        pop_interp="$interpolated_dir/tau_POP_M${m}.dat"
        
        # Check if all three files exist
        if [ -f "$iis_file" ] && [ -f "$matth_file" ] && [ -f "$pop_file" ]; then
            echo "  Computing total scattering for M${m}"
            
            # Create a Python script for interpolation and combination
            python3 << EOF
import numpy as np
from scipy.interpolate import interp1d
import sys
import os

# Interpolation parameters
E_min = $E_MIN
E_max = $E_MAX
N_steps = $N_STEPS

# Create common energy grid
E_common = np.linspace(E_min, E_max, N_steps)

# Function to load data and interpolate
def load_and_interpolate(filename, E_common):
    try:
        data = np.loadtxt(filename)
        if data.ndim == 1:
            # Only one data point
            E = np.array([data[0]])
            tau = np.array([data[1]])
        else:
            E = data[:, 0]
            tau = data[:, 1]
        
        # Create interpolation function
        # Use linear interpolation with extrapolation
        f = interp1d(E, tau, kind='linear', bounds_error=False, 
                     fill_value=(tau[0], tau[-1]))
        tau_interp = f(E_common)
        return tau_interp
    except Exception as e:
        print(f"Error loading {filename}: {e}", file=sys.stderr)
        return np.zeros(N_steps)

# Load and interpolate all three mechanisms
tau_matth = load_and_interpolate("$matth_file", E_common)
tau_pop = load_and_interpolate("$pop_file", E_common)
tau_iis = load_and_interpolate("$iis_file", E_common)

# Save interpolated data
np.savetxt("$matth_interp", np.column_stack((E_common, tau_matth)), 
           fmt='%.6e', header='Energy(eV) Tau_matth(fs)')
np.savetxt("$pop_interp", np.column_stack((E_common, tau_pop)), 
           fmt='%.6e', header='Energy(eV) Tau_POP(fs)')
np.savetxt("$iis_interp", np.column_stack((E_common, tau_iis)), 
           fmt='%.6e', header='Energy(eV) Tau_IIS(fs)')

# Compute total scattering using Matthiessen's rule
# Avoid division by zero
tau_total = np.zeros(N_steps)
mask = (tau_matth > 0) & (tau_pop > 0) & (tau_iis > 0)
tau_total[mask] = 1.0 / (1.0/tau_matth[mask] + 1.0/tau_pop[mask] + 1.0/tau_iis[mask])

# Save total scattering
np.savetxt("$output_file", np.column_stack((E_common, tau_total)), 
           fmt='%.6e', header='Energy(eV) Tau_total(fs)')

print(f"    Interpolated {N_steps} points from {E_min} to {E_max} eV")
print(f"    Non-zero total tau points: {np.sum(mask)}")
print(f"    Max total tau: {np.max(tau_total):.4f} fs")
print(f"    Min total tau (non-zero): {np.min(tau_total[mask]):.6f} fs" if np.sum(mask) > 0 else "    All points are zero")
EOF
            
            if [ -s "$output_file" ]; then
                lines=$(wc -l < "$output_file")
                echo "    Created $output_file ($lines data lines on common grid)"
            else
                echo "    WARNING: Output file is empty!"
            fi
        else
            echo "  Warning: Missing files for M${m} at temperature $temp"
        fi
    done
    
    echo "  Completed total scattering for $temp"
done

echo ""
echo "Total scattering processing with interpolation complete!"
