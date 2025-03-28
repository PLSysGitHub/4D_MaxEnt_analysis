# Learning the dynamic organization of a replicating bacterial chromosome from time-course Hi-C data
## Analysis code

This repository contains code for analyzing chromosome configurations sampled from the 4D-MaxEnt model for replicating bacterial chromosomes. The same analyses are also performed for samples from mechanistic molecular dynamics simulations of loop-extruders and/or ori-pulling.

## Donwloading raw data

The raw data files with simulated polymer configurations can be downloaded from Zenodo, under the accession code [10.5281/zenodo.15099820](https://doi.org/10.5281/zenodo.15099820).

## Running the code

- This code is written in Julia, and was ran using Julia 1.11. 
- To install the required packages, start julia, and run the following commands:
    
    ```julia
    ]
    activate .
    instantiate
    ```

- After installing, you can run any of the included scripts, which do the following:

### Analysis for model construction

- Find_fork_positions_form_Hi-C -> Use Hi-C counts to infer fork positions

### Analysis of model statistics

- process_MaxEnt.jl -> Convert .txt files to .h5 format. Calculate and save statistics for 4D-MaxEnt model.

- process_MD.jl -> calculate statistics for bottom-up mechanistic simulations

- plot_results.jl -> plot all results and save plots

- check_MD_time_scales.jl -> check the time-scales of the MD simulations

- check_loop_size_distributions.jl -> check the loop size distributions of the MD simulations

- check_initial_loop_distributions.jl -> check how many generations needed before initial loop size distributions converged

- analyze_MaxEnt_local_extension.jl -> Average local extensions behind and ahead of forks, fit exponential decay

## Questions?

If you have any questions, feel free to get in touch with j.k.harju [at] vu.nl
