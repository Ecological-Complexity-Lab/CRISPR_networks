# UNDER DEVELOPMENT!!

This repository accompanies the paper Pilosof S, Alcalá-Corona SA, Wang T, Kim T, Maslov S. The network structure and eco-evolutionary dynamics of CRISPR-induced immune diversification. bioRxiv. 2020. Available: <https://www.biorxiv.org/content/10.1101/850800v1.abstract>. It contains the code and data for the simulation and empirical analysis included in the paper.

# Simulated data

## Simulations based on the Childs et. al. model

Here we implement the hybrid (deterministic/stochastic) model by Lauren M. Childs, et. al. in *"Multiscale model of CRISPR‐induced coevolutionary dynamics: diversification at the interface of Lamarck and Darwin"*. [*Evolution: International Journal of Organic Evolution 66, no. 7 (2012)*](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.2012.01595.x).


The C++ source code to make the simulations can be found [here](https://github.com/Ecological-Complexity-Lab/CRISPR_networks/blob/master/SourceCodeSimulatior/CRISPR_Model.cpp). 

### Compilation

The code is written in C++ (version 11). You must compile the code in your command line by using a proper C++11 compiler.
For example:
```
g++ -std=c++11 CRISPR_Model.cpp -o CRISPRsimulator
```

### Run

Run the executable compiled code in the command line, using the next parameters:
```
./CRISPRsimulator Dp mu runT Pts Sp seed 
```
Where:
* `Dp`: is the initial number of phage strains. Usually is set to 1 virus strain.
* `mu`: The prtospacer mutation rate per protospacer.
* `runT`: Specifies the running time-length as hours for the simulation. This is the maximum length of each simulation. The execution time will be different as the other parameters are varied. 
* `Pts`: the Number of Protospacers for the virus strains.
* `Sp`: Specifies the Number of Spacers for the host strains.
* `seed`: Set the seed for the random generator used for the stochastic part.

For example:
```
./CRISPRsimulator 1 1e-7 5000 15 10 12499
```

See [here](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.2012.01595.x) for more details about the description of the parameters. 


You can also set the seed of simulation, using a random number generated from your Unix/Linux System. This is for running different simulations with the same parameters.
```
./CRISPRsimulator 1 1e-7 5000 15 10 $RANDOM
```

#### Output files

Once the code has finished it will generate 7 output files, which will be ready for analysis which is detailed below.

### Analysis of the output files

In this repository we include the output of this simulator as input data, and the R code necessary to analyze these data. The simulator produces files named `mu1e-7_initialDiffDp1_S10P15_R-NNNN`, where NNNN is the seed number of a particular simulation. We have included the stochastic simulation outputs for our main simulation example (seed number 12499) in the folder `data\mu1e-7_initialDiffDp1_S10P15_R-12499`.

* `Bacteria-abundance.txt`: The density of bacteria at different times.
* `Phage-abundance.txt`: The density of viruses at different times.
* `Bacteria-TREE.txt`: Specifies the parents and children of bacteria. Converted to a nwk file in section `Trees` in file `simulations_analysis.R`.
* `data-bact.txt`: Spacer composition of bacteria strains at each time.
* `data-phage.txt`: Protospacer composition of virus strains at each time.

We also include a list with all the seeds for the 100 simulations we ran in file `seed_index.csv`. Running the simulator with these seeds will reproduce the exact data we use in the paper.

## Analysis of simulations

Output of the stochastic simulations is the input of the file `simulations_analysis.R`. We also include the file `host_spacer_simulation.Rmd` for specific analysis of modularity and phylogenetic distance in host-spacer and infection networks.

### Requirements

* The code is written in R and requires the packages specified in the code. Modularity analysis is performed with package [infomapecology](https://github.com/Ecological-Complexity-Lab/infomap_ecology_package). See instructions for installation there.
* The analysis is computationally intensive and we ran it on the Midway HPC cluster at the University of Chicago. Some adaptation of the R code will be necessary to adequate it to other HPC systems.
* Section Set up in file `simulations_analysis.R` include code that may need to be adapted and the R packages.

### Analysis breakdown

The R file `simulations_analysis.R` is divided into sections.

* **Set up**: The first line includes the parameter values used (Table is in the paper). Lines 11-24 are used when calling the file from an external sbatch file when running jobs on a HPC system.
* **Initialize**: Read simulation data and use it to define regimes. Also prepare the data frames.
* **Functions**: Define functions. The main function is `create_networks_hr`, which creates all the networks (Fig 2 in the paper).
* **bacteria / phage diversity**: Plot abundance profiles of hosts and viruses.
* **Generate networks**: Generate networks for each time.
* **Measures of diversity**: Calculate host, virus and spacer abundance and richness. Other indices of diversity can be included in this section.
* **Phage and bacteria diversification**: Calculate the persistence of hosts and viruses.
* **Trees**: Create the phylogenetic trees for hosts and viruses.
* **Modularity of infection networks**: Calculate the modularity of infection networks in time. If `make_plots` is TRUE then it will also produce a plot of the infection network per time step (good for exploratory analysis).
* **WNODF**: Calculate the WNDOF index for weighted nestedness in immunity networks at each time.
* **Spacer matches**: Calculate the proportion of spacer matches (the edge weights in immunity networks) at each time.
* **Extinctions**: Analyze order of virus extinction.
* **R0 for 0 and 1 matches**: This calculates the R0 and R1, which together form the R_pot (Figure 4)
* **Plot**: Produce a PDF with all the plots for convenience.

The following sections are also used in file `host_spacer_simulation.Rmd`. The analysis in these sections also serves for comparison to empirical data.

* **Significance of modularity of infection networks**: Calculate the significance of modularity at the end of each VDR by comparing to shuffled networks.
* **Phylogenetic signal in infection networks**: Is there a phylogenetic signal in infection networks at the end of each VDR?
* **Significance of modularity of host-spacer networks**: Aggregates host-spacer networks within each VDR and test if modularity is non-random compared to shuffled networks.
* **Phylogenetic signal in host-spacer modules**: Is there a phylogenetic signal in the host-spacer networks analyzed in the previous section?

# Empirical data

## Database

The data sets we use in the paper are stored in an SQL database `CRISPR_database_NEE.sqlite`.

## Analysis

Code is concentrated in file `empirical_data_analysis.R`. We also include `empirical_data_analysis.Rmd`, which sources the R file to create a markdown file for convenience. The main functions that perform the analysis are:

* **get_strain_spacers**: Obtains the spacer set of a virus or a host from the host-spacer/virus-protospacer matrix.
* **get_matches**: Match spacers and protospacers and produce an immunity network.
* **test_PD_modules**: the same function as in the simulated data, to test for phylogenetic signal.
* **calculate_ev_nestedness**: Calculate leading eigenvalue weighted nestedness (Staniczenko PPA, Kopp JC, Allesina S. The ghost of nestedness in ecological networks. Nat Commun. 2013;4: 1391.).
* **main**: A wrapper to read data, build and analyze immunity networks.

The analysis of modularity is also performed with infomapecology.


# Neutral model without explicit immunity

Here we implement a mean-field simplification of the model which does not explicitly include the diversity of hosts and viruses. Instead, we consider the model as a Lotka-Volterra model whit all the host strains considered as one single prey, and all the viruses strains are collapsed in one single predator. The immunity is taken as an average value of the matrix M which is artificial.

The Python code is just a simple ODE integrator of the LV mean-field simplified system, which uses different values of the average immunity matrix M.

The code found [here](https://github.com/Ecological-Complexity-Lab/CRISPR_networks/blob/master/MeanField-NeutralDynamics/MFDyn_AverageImmMatrix.py). 

## Run

```
python MFDyn_AverageImmMatrix.py M1 M2 M3 M4
```
Where M<sub>i</sub> are the artificial average values of the matrix M.

For example:
```
python MFDyn_AverageImmMatrix.py 0.35 0.7 0.85 0.9167
```
