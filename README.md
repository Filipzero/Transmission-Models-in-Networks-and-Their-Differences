# Disease-Transmission-Models-in-Networks-and-Their-Differences
This repository is built in Python, featuring six epidemiological models implemented on four different underlying networks. The epidemiological compartmental models used, are SI, SIS, SIR, SIRS, SEIR, SEIRS and the underlying networks are Erdős–Rényi, Watts-Strogatz, Barabási–Albert and finally a 2x2 lattice. Each simulation provides the results of individual stochastic runs as well as their averaged outcome across multiple simulations. These averages are then compared with the predictions of the corresponding deterministic (mean-field) ODE model. The simulations were carried out with the EoN (Epidemics on Networks) library, which employs the Gillespie stochastic algorithm and the networks were created with the NetworkX library. 

## Model Description

The models describe a population divided into compartments (depending on the epidemic model) and are implemented in two ways:

* As a **mean-field ODE approximation**, which assumes homogeneous mixing of the population and describes the evolution of the proportions of susceptible and infected individuals over time. The total population N is assumed constant, so S(t)+I(t)=N, or in terms of proportions, s(t)+i(t)=1.

* As a **stochastic network-based simulation** using the Gillespie algorithm, which accounts for discrete individuals and network structure, allowing for variability and local effects of infection spread, while keeping the total population fixed at N.

### SI (Susceptible-Infected)

The SI model describes a population that is divided into two compartments: Susceptible(S) and Infected (I). Susceptible individuals can become infected and thus infectious through contact with another individual who has the disease and they remain infected forever.

The dynamics of this model can be described by the following mean-field ODEs:

$$ \frac{dS}{dt} = -\frac{β S I}{N} $$
$$ \frac{dI}{dt} = \frac{β S I}{N} $$

where:
* S(t) is the number of susceptible individuals at time t
* I(t) is the number of infected individuals at time t
* N is the total population, which remains constant (S+I=N)
* β is the infection rate, representing the probability per unit time that a susceptible individual becomes infected upon contact with an infected individual.  
<img width="1175" height="695" alt="image UCRTE3" src="https://github.com/user-attachments/assets/8107b70c-684a-4a32-b440-4ab6e92cb68c" />

$$ Transition\ from\ susceptible\ to\ infected. Full\ conversion\ occurs,\ as\ predicted\ by\ the\ ODEs,\ for\ N=10,000\ and\ β=0.3. $$

### SIS (Susceptible-Infected-Susceptible)

The SIS model describes a population that is divided into two compartments: Susceptible(S) and Infected (I). The difference between the SI model is that the infected individuals can recover and return to the susceptible state, allowing reinfection.

The dynamics here will change taking into account the reinfection:

$$ \frac{dS}{dt} = -\frac{β S I}{N} + γΙ $$
$$ \frac{dI}{dt} = \frac{β S I}{N} -γI $$

where:
* γ is the rate at which infected individuals return to the susceptible state.
<img width="1184" height="739" alt="image CVSOE3" src="https://github.com/user-attachments/assets/53b54594-1132-4e13-ac9d-d0a66a32b58e" />
text{Transition\ from\ susceptible\ to\ infected.\ The\ system\ reaches\ a\ stable\ balance\ between\ S-I\ with\ an\ infection\ rate\ β=0.3\ \and\ "recovery"\ rate\ γ=0.1.}






