# Disease-Transmission-Models-in-Networks-and-Their-Differences
This repository is built in Python, featuring six epidemiological models implemented on four different underlying networks. The epidemiological compartmental models used, are SI, SIS, SIR, SIRS, SEIR, SEIRS and the underlying networks are Erdős–Rényi, Watts-Strogatz, Barabási–Albert and finally a 2x2 lattice. Each simulation provides the results of individual stochastic runs as well as their averaged outcome across multiple simulations. These averages are then compared with the predictions of the corresponding deterministic (mean-field) ODE model. The simulations were carried out with the EoN (Epidemics on Networks) library, which employs the Gillespie stochastic algorithm. 

## Model Description

The models are implemented both:

* As a **mean-field ODE approximation**, which assumes homogeneous mixing of the population and describes the evolution of the proportions of susceptible and infected individuals over time.
* As a **stochastic network-based simulation** using the Gillespie algorithm, which accounts for discrete individuals and network structure, allowing for variability and local effects of infection spread.

#### SI
