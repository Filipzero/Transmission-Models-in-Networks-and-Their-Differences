# Disease-Transmission-Models-in-Networks-and-Their-Differences
This repository is built in Python, featuring six epidemiological models implemented on four different underlying networks. The epidemiological compartmental models used, are SI, SIS, SIR, SIRS, SEIR, SEIRS and the underlying networks are Erdős–Rényi, Watts-Strogatz, Barabási–Albert and finally a 2x2 lattice. Each simulation provides the results of individual stochastic runs as well as their averaged outcome across multiple simulations. These averages are then compared with the predictions of the corresponding deterministic (mean-field) ODE model. The simulations were carried out with the EoN (Epidemics on Networks) library, which employs the Gillespie stochastic algorithm and the networks were created with the NetworkX library. 

## Model Description

The models describe a population divided into compartments (depending on the epidemic model) and are implemented in two ways:

* As a **mean-field ODE approximation**, which assumes homogeneous mixing of the population and describes the evolution of the proportions of susceptible and infected individuals over time. The total population N is assumed constant, so S(t)+I(t)=N, or in terms of proportions, s(t)+i(t)=1.

* As a **stochastic network-based simulation** using the Gillespie algorithm, which accounts for discrete individuals and network structure, allowing for variability and local effects of infection spread, while keeping the total population fixed at N.

### S-I (Susceptible-Infected)

The SI model describes a population that is divided into two compartments: Susceptible(S) and Infected (I). Susceptible individuals can become infected and thus infectious through contact with another individual who has the disease and they remain infected forever.

The dynamics of this model can be described by the following mean-field ODEs:

$$ \frac{dS}{dt} = -\frac{β S I}{N} $$
$$ \frac{dI}{dt} = \frac{β S I}{N} $$
where:
* ds/dt 
