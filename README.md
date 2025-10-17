# Disease Transmission Models in Networks and Their Differences
This repository is built in Python, featuring six epidemiological models implemented on four different underlying networks. The epidemiological compartmental models used, are SI, SIS, SIR, SIRS, SEIR, SEIRS and the underlying networks are Erdős–Rényi, Watts-Strogatz, Barabási–Albert and finally a 2x2 lattice. Each simulation provides the results of individual stochastic runs as well as their averaged outcome across multiple simulations. These averages are then compared with the predictions of the corresponding deterministic (mean-field) ODE model, as well as the results obtained from different network topologies. The simulations were carried out with the EoN (Epidemics on Networks) library, which employs the Gillespie stochastic algorithm and the networks were created with the NetworkX library. 

## Table of contents
1. [Requirements](#requirements)
2. [How to Run the Simulation](#how-to-run-the-simulation)
3. [Compartmental Model Description](#compartmental-model-description)
   - [SI](#si-susceptible-infected)
   - [SIS](#sis-susceptible-infected)
   - [SIR](#sir-susceptible-infected-recovered)
   - [SIRS](#sirs-susceptible-infected-recovered)
   - [SEIR](#seir-susceptible-exposed-infected-recovered)
   - [SEIRS](#seirs-suscpetible-exposed-infected-recovered)
4. [Graph Topologies](#graph-topologies)

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

The models describe a population divided into X compartments (depending on the epidemic model) and are implemented in two ways:

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

$$ Parameters\ given:\ β=0.3,\ σ=0.2,\ γ=0.1,\ ω=0.05.\ Similar\ final\ state\ to\ that\ of\ the\ other\ models\ with\ waning\ immunity. $$

## Graph Topologies

1. Erdős–Rényi (Random) Network

   * Each pair of nodes is connected with a fixed probability p.
   * Useful for studying random interactions in large populations.

2. Watts–Strogatz (Small-World) Network

   * Nodes are initially connected in a regular lattice, with some edges randomly rewired.
   * Captures real-world networks with high clustering and short path lengths.

3. Barabási–Albert (Scale-Free) Network

   * Nodes are added sequentially with preferential attachment, resulting in hubs.
   * Reflects networks with heterogeneous connectivity, like social or contact networks.

4. Grid Network (Optional)

   * Nodes are placed on a 2D lattice with periodic boundary conditions (the edges “wrap around,” so corner nodes can interact across the boundary).

Each network type can influence the dynamics of disease spread, affecting infection speed and final outbreak size.

## Gillespie algorithm
The simulations for the disease spreading in networks were implemented using the [EoN](https://epidemicsonnetworks.readthedocs.io/en/latest/EoN.html) library. The algorithm that is used simulates the transitions between compartments by following a stochastic process for selecting the next event in time. More specifically, I will describe how the algorithm operates for the SIR model, with only minor differences for other models.

The algorithm begins by calculating the transition rates. Each susceptible node has a probability of becoming infected depending on how many infected neighbors it has. That is, if a node has one infected neighbor, its infection rate is **τ** (the transmission rate per edge), while if it has two infected neighbors, the rate becomes **2τ**, and so on. The probability that an infected node recovers is constant and independent of what happens around it. From these, we obtain the total rate **Τ**, which is the sum of all infection and recovery rates in the network.

The algorithm having computed the total rate, it must now decide when the next event will occur and what that event will be. The time until the next event is drawn from an exponential distribution with rate **T**. This means that if **T** is large, an event will occur soon, while if **Τ** is small, it will occur later. The following step for the algorithm is to select which event takes place, with a probability proportional to its rate. Since all recoveries occur at the same constant rate, the total recovery rate is **γI** and the probability that the next event is a recovery equals **γI/Τ**. Therefore, if a recovery occurs, we choose a random infected node to recover.

If not, an infection event is selected instead. In that case, it must determine which susceptible node becomes infected, with probabilities proportional to their individual infection rates. The complication here is that these rates differ among nodes, unlike the constant recovery rate, because each susceptible node’s risk depends on how many infected neighbors it has. To handle this, it draws a random number **r∈(0, T−γI)** and iterates over all susceptible nodes, subtracting their infection rate from **r** until it becomes negative. The node that causes this to happen is the one that becomes infected.

An alternative approach is to define a fixed maximum infection rate **$r^*$** greater than or equal to the actual maximum infection rate among all nodes **$r_u$** thereby increasing the total rate **Τ**. When a transmission is scheduled to occur, the algorithm randomly selects a susceptible node whose real infection rate is less than or equal to **$r^*$** . However, since it temporarily overestimated the infection probabilities, it generates a new random number **p∈(0, 1)**; if **p < ($r_u$ / $r^*$)** , the infection of node u is accepted. Otherwise, the time advances to the next event without a change.





















