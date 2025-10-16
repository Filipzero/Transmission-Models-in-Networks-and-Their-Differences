import numpy as np
import networkx as nx
import EoN
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd

# Parameters
N = 10000  # network size 
p = 0.01  # Probability of edge creation
m = 3 # nodes
#k = p * (N - 1)
beta = 0.3  # Infection rate
#τ = beta / k
t_max = 100  # Simulation duration
num_simulations = 10  # Number of simulations to average
report_times = np.linspace(0, t_max, 100)


#G = nx.grid_2d_graph(100, 100, periodic=True)
G = nx.fast_gnp_random_graph(N, p)
#G = nx.watts_strogatz_graph(N, 10, 0.1)  # N nodes, each connected to 10 neighbors, 10% rewiring - Small-World Network
#G = nx.barabasi_albert_graph(N, m)

avg_degree = sum(dict(G.degree()).values()) / N
τ = beta / avg_degree



# ODEs
def SI_ODE(t, y, beta):
    S, I = y
    dSdt = -beta * S * I
    dIdt = beta * S * I
    return [dSdt, dIdt]

# Initialize arrays for storing simulation results
S_means = np.zeros(len(report_times))
I_means = np.zeros(len(report_times))

# Initial conditions
S0, I0 = 0.99, 0.01

# Solve the ODE model
ode_solution = solve_ivp(
    SI_ODE,
    [0, t_max],
    [S0, I0],
    args=(beta,),
    t_eval=report_times) #Runge-Kutta επιλυση
t_ode = ode_solution.t
S_ode, I_ode = ode_solution.y

# Create transition graphs for the SI model
# Neighbor-induced transitions (I+S -> I+I)
J = nx.DiGraph()
J.add_edge(('I', 'S'), ('I', 'I'), rate=τ)

H = nx.DiGraph()  # No spontaneous transitions for the SI model



fig, ax = plt.subplots(figsize=(14, 8))


for sim in range(num_simulations):

    # Initialize node states
    IC = {}
    nodes = list(G.nodes())
    np.random.shuffle(nodes)
    
    n_infected = int(I0 * N)
    
    for i, node in enumerate(nodes):
        if i < n_infected:
            IC[node] = 'I'
        else:
            IC[node] = 'S'

    # Run the Gillespie simulation
    t, S, I = EoN.Gillespie_simple_contagion(
        G,
        H,  # Empty spontaneous transition graph
        J,  # Neighbor-induced transition graph
        IC,  # Initial conditions
        return_statuses=('S', 'I'),
        tmax=t_max
    )

    # Subsample the results
    S_subsampled = np.interp(report_times, t, S)
    I_subsampled = np.interp(report_times, t, I)
    
    ax.plot(report_times, S_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    ax.plot(report_times, I_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    
    # Add to means
    S_means += S_subsampled / N
    I_means += I_subsampled / N

# Average the simulation results
S_means /= num_simulations
I_means /= num_simulations

# Plot

ax.scatter(report_times, S_means, label="Simulation (S)", color="blue", linewidth=2, s=7)
ax.scatter(report_times, I_means, label="Simulation (I)", color="red", linewidth=2, s=7)

ax.plot(t_ode, S_ode, label="ODE (S)", color="navy", linewidth=2, alpha=0.6)
ax.plot(t_ode, I_ode, label="ODE (I)", color="darkred", linewidth=2, alpha=0.6)


ax.set_title("SI Model on Erdos-Renyi vs Mean-Field Approximation")
ax.set_xlabel("Time (t)", fontsize=14, labelpad=15, color="#333333")
ax.set_ylabel("Proportion of Population", fontsize=14, labelpad=15, color="#333333")
ax.tick_params(axis='both', which='major', labelsize=12, colors="#333333")
ax.legend(fontsize=14, loc="best", frameon=True, shadow=True)

# Customize grid and ticks
ax.grid(color='lightgray', linestyle='-', linewidth=1.5)
ax.set_axisbelow(True)

# Customize plot background
for spine in ax.spines.values():
    spine.set_edgecolor('#bbbbbb')
    
#fig.savefig("SIS_Model.png", dpi=300, bbox_inches="tight")  # Save as Png
plt.show()


# ******** |STORING DATA TO CSV| ********

#data = {
#        "t": report_times,
#        "S_mean": S_means,
#        "I_mean": I_means
#        }

#df_SI = pd.DataFrame(data)
#print(df)

#df_SI.to_csv('./SI.csv', sep=';', index=False) # if ';' doesn't work then ','
























