# Disease Transmission Models in Networks and Their Differences
This repository is built in Python, featuring six epidemiological models implemented on four different underlying networks. The epidemiological compartmental models used, are SI, SIS, SIR, SIRS, SEIR, SEIRS and the underlying networks are Erdős–Rényi, Watts-Strogatz, Barabási–Albert and finally a 2x2 lattice. Each simulation provides the results of individual stochastic runs as well as their averaged outcome across multiple simulations. These averages are then compared with the predictions of the corresponding deterministic (mean-field) ODE model. The simulations were carried out with the EoN (Epidemics on Networks) library, which employs the Gillespie stochastic algorithm and the networks were created with the NetworkX library. 

1. [Requirements](#1-requirements)
2. [How to Run the Simulation](#2-how-to-run-the-simulation)
3. [Compartmental Model Description](#3-compartmental-model-description)
   - [SI](#si-susceptible-infected)
   - [SIS](#sis-susceptible-infected)
   - [SIR](#sir-susceptible-infected-recovered)
   - [SIRS](#sirs-susceptible-infected-recovered)
   - [SEIR](#seir-susceptible-exposed-infected-recovered)
   - [SEIRS](#seirs-suscpetible-exposed-infected-recovered)

## Requirements
This project requires Python 3.8+ and the libraries listed in requirements.txt.

## How to Run the Simulation
1. Clone or download the project repository.
2. Make sure the libraries in the requirements.txt file are installed
3. Run a script using an IDE or the terminal like so :
   python3 EoN_SI.py
4. Optional: Adjust simulation parameters in the script before running, such as:

    Network type (G)
    Network size (N)
    Infection rate (beta)
    Simulation duration (t_max)
    Number of simulations (num_simulations)
    Other related parameters

5. The scripts will display plots of the SI, SIS, SIR, SIRS, SEIR, SEIRS models and optionally save results to CSV if enabled.

## Compartmental Model Description

The following parameters are used across the models:

| Symbol | Description |
|--------|-------------|
| N      | Total population (constant) |
| S(t)   | Number of susceptible individuals at time t |
| I(t)   | Number of infected individuals at time t |
| R(t)   | Number of recovered individuals at time t (if applicable) |
| E(t)   | Number of exposed individuals at time t (if applicable) |
| β      | Infection rate: probability per unit time that a susceptible individual becomes infected upon contact with an infected individual |
| γ      | Recovery rate: rate at which infected individuals recover (**SIR, SIRS, SEIR, SEIRS**) | 
| σ      | Progression rate from exposed to infectious (SEIR, SEIRS) |
| ω      | Waning immunity rate (SIRS, SEIRS) |

The models describe a population divided into compartments (depending on the epidemic model) and are implemented in two ways:

* As a **mean-field ODE approximation**, which assumes homogeneous mixing of the population and describes the evolution of the proportions of susceptible and infected individuals over time. The total population N is assumed constant, so S(t)+I(t)=N, or in terms of proportions, s(t)+i(t)=1.

* As a **stochastic network-based simulation** using the Gillespie algorithm, which accounts for discrete individuals and network structure, allowing for variability and local effects of infection spread, while keeping the total population fixed at N.

### SI (Susceptible-Infected)

The SI model describes a population that is divided into two compartments: Susceptible(S) and Infected (I). Susceptible individuals can become infected and thus infectious through contact with another individual who has the disease and they remain infected forever.

The dynamics of this model can be described by the following mean-field ODEs:

$$ \frac{dS}{dt} = -\frac{β S I}{N} $$
$$ \frac{dI}{dt} = \frac{β S I}{N} $$
 
<img width="1175" height="695" alt="image UCRTE3" src="https://github.com/user-attachments/assets/8107b70c-684a-4a32-b440-4ab6e92cb68c" />

$$ Transition\ from\ susceptible\ to\ infected. Full\ conversion\ occurs,\ as\ predicted\ by\ the\ ODEs,\ for\ N=10,000\ and\ β=0.3. $$

### SIS (Susceptible-Infected)

The SIS model describes a population that is divided into two compartments: Susceptible(S) and Infected (I). The difference between the SI model is that the infected individuals can recover and return to the susceptible state, allowing reinfection.

The dynamics here will change taking into account the reinfection:

$$ \frac{dS}{dt} = -\frac{β S I}{N} + γΙ $$
$$ \frac{dI}{dt} = \frac{β S I}{N} -γI $$

where:
* γ is the rate at which infected individuals return to the susceptible state (**different from the other γ parameters in the other models**).
<img width="1184" height="739" alt="image CVSOE3" src="https://github.com/user-attachments/assets/53b54594-1132-4e13-ac9d-d0a66a32b58e" />

$$ Transition\ from\ susceptible\ to\ infected.\ The\ system\ reaches\ a\ stable\ balance\ between\ S-I\ with\ an\ infection\ rate\ β=0.3\ and\ "recovery"\ rate\ γ=0.1. $$

### SIR (Susceptible-Infected-Recovered)

The population is divided in three compartments this time, meaning that an individual can transistion through three stages: Susceptible(S), Infected(I) and Recovered(R).

The mean-field ODEs:

$$ \frac{dS}{dt} = -\frac{β S I}{N} $$
$$ \frac{dI}{dt} = \frac{β S I}{N} -γI $$
$$ \frac{dR}{dt} = γI $$

<img width="1171" height="703" alt="image SRKUE3" src="https://github.com/user-attachments/assets/b5633ce1-3dc6-4ba6-b131-374ffeebaaee" />

$$ Susceptibles\ become\ infected\ and\ then\ recover.\ With\ parameters\ β=0.3\ and\ γ=0.1, we\ can\ see\ that\ not\ the\ entire\ population\ becomes\ infected.\ Some\ remain\ susceptible\ and\ thats\ why\ the\ recovered\ curve\ does\ not\ reach\ 1. $$

### SIRS (Susceptible-Infected-Recovered)

The population remains divided in three compartments, but as in the SIS and SI models reinfection is possible. The immunity that they gained can now be lost over time.

The ODEs of this compartmental models are:

$$ \frac{dS}{dt} = -\frac{β S I}{N} + ωR $$
$$ \frac{dI}{dt} = \frac{β S I}{N} -γI $$
$$ \frac{dR}{dt} = γI -ωR $$

<img width="1180" height="732" alt="image 82ZOE3" src="https://github.com/user-attachments/assets/230d67ef-2479-4693-8c3e-69b4008e14d1" />

$$ We\ clearly\ see\ that\ the\ system\ reaches\ a\ stable\ balance\ between\ S-I-R (β=0.3,\ γ=0.1,\ ω=0.05\).\ Note\ that\ if\ we\ had\ different\ rates\ the\ system\ would\ still\ reach\ equilibrium\ but\ the\ peak\ values\ would\ differ. $$

### SEIR (Susceptible-Exposed-Infected-Recovered)

In this model a new compartment is introduced. It accounts for a latent period (E) before individuals become infectious. Individuals first enter an exposed state, during which they are infected but not yet infectious, before progressing to the infectious stage.

For the ODEs we get:

$$ \frac{dS}{dt} = -\frac{β S I}{N} $$
$$ \frac{dΕ}{dt} = \frac{β S I}{N} -σΕ $$
$$ \frac{dΙ}{dt} = σE - γI $$
$$ \frac{dR}{dt} = γI $$

<img width="1179" height="732" alt="image HTASE3" src="https://github.com/user-attachments/assets/bb4db37e-0d2b-4a47-a956-3eec4bb612a9" />

$$ Parameters\ given:\ β=0.3,\ σ=0.2,\ γ=0.1.\ The\ system's\ dynamics\ change,\ making\ the\ infectious\ state\ progress\ slower. $$

### SEIRS (Susceptible-Exposed-Infected-Recovered)

This model allows waning immunity, making reinfection possible, so its ODEs will be:

$$ \frac{dS}{dt} = -\frac{β S I}{N} +ωR $$
$$ \frac{dΕ}{dt} = \frac{β S I}{N} -σΕ $$
$$ \frac{dΙ}{dt} = σE - γI $$
$$ \frac{dR}{dt} = γI - ωR $$

<img width="1176" height="731" alt="image ZCOUE3" src="https://github.com/user-attachments/assets/312420a5-8627-4eff-b2e1-304fc38189e6" />

$$ Parameters\ given:\ β=0.3,\ σ=0.2,\ γ=0.1, ω=0.05.\ Similar\ final\ state\ to\ that\ of\ the\ other\ models\ with\ waning\ immunity. $$
