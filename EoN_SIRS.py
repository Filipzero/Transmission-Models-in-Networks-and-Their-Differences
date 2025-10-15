import numpy as np
import networkx as nx
import EoN
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd


# Parameters
N = 10000  # network size 
p = 0.01  # Probability of edge creation in Erdős–Rényi network
m = 3 # nodes
#k = p * (N-1)
beta = 0.3  # Infection rate
#τ = beta / k
gamma = 0.1  # Recovery rate
omega = 0.05  # Rate of immunity loss (R -> S transition)
t_max = 180  # Simulation duration
num_simulations = 10  # Number of simulations to average
report_times = np.linspace(0, t_max, 100)

#G = nx.grid_2d_graph(100, 100, periodic=True)
#G = nx.fast_gnp_random_graph(N, p)
#G = nx.watts_strogatz_graph(N, 10, 0.1)  # N nodes, each connected to 10 neighbors, 10% rewiring
G = nx.barabasi_albert_graph(N, m)

avg_degree = sum(dict(G.degree()).values()) / N
τ = beta / avg_degree



# ODEs
def SIRS_ODE(t, y, beta, gamma, omega):
    S, I, R = y
    dSdt = -beta * S * I + omega * R
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I - omega * R
    return [dSdt, dIdt, dRdt]

# Initialize arrays for storing simulation results
S_means = np.zeros(len(report_times))
I_means = np.zeros(len(report_times))
R_means = np.zeros(len(report_times))

# Initial conditions
S0, I0, R0 = 0.99, 0.01, 0.0

# Solve the ODE model
ode_solution = solve_ivp(
    SIRS_ODE, 
    [0, t_max], 
    [S0, I0, R0], 
    args=(beta, gamma, omega), 
    t_eval=report_times) #Runge-Kutta επιλυση
t_ode = ode_solution.t
S_ode, I_ode, R_ode = ode_solution.y

# Create transition graphs for SIRS model
# Spontaneous transitions (I->R, R->S)
H = nx.DiGraph()
H.add_node('S')  # Need to add isolated nodes explicitly
H.add_edge('I', 'R', rate=gamma)  # I -> R transition
H.add_edge('R', 'S', rate=omega)    # R -> S transition

# Neighbor-induced transitions (I+S->I+I)
J = nx.DiGraph()
J.add_edge(('I', 'S'), ('I', 'I'), rate=τ)

fig, ax = plt.subplots(figsize=(14, 8))

for sim in range(num_simulations):
    
    # Initialize node states
    IC = {}
    nodes = list(G.nodes())
    np.random.shuffle(nodes)
    
    n_infected = int(I0 * N)
    n_recovered = int(R0 * N)
    
    for i, node in enumerate(nodes):
        if i < n_infected:
            IC[node] = 'I'
        elif i < n_infected + n_recovered:
            IC[node] = 'R'
        else:
            IC[node] = 'S'


    t, S, I, R = EoN.Gillespie_simple_contagion(
        G, 
        H,  # spontaneous transition graph
        J,  # neighbor-induced transition graph
        IC, # initial conditions
        return_statuses=('S', 'I', 'R'),
        tmax=t_max
    )

    # Subsample the results
    S_subsampled = np.interp(report_times, t, S)
    I_subsampled = np.interp(report_times, t, I)
    R_subsampled = np.interp(report_times, t, R)
    
    ax.plot(report_times, S_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    ax.plot(report_times, I_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    ax.plot(report_times, R_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    
    # Add to means
    S_means += S_subsampled / N
    I_means += I_subsampled / N
    R_means += R_subsampled / N

# Average the simulation results
S_means /= num_simulations
I_means /= num_simulations
R_means /= num_simulations

# Plot
ax.scatter(report_times, S_means, label="Simulation (S)", color="blue", linewidth=2, s=7)
ax.scatter(report_times, I_means, label="Simulation (I)", color="red", linewidth=2, s=7)
ax.scatter(report_times, R_means, label="Simulation (R)", color="green", linewidth=2, s=7)

ax.plot(t_ode, S_ode, label="ODE (S)", color="navy", linewidth=2, alpha=0.6)
ax.plot(t_ode, I_ode, label="ODE (I)", color="darkred", linewidth=2, alpha=0.6)
ax.plot(t_ode, R_ode, label="ODE (R)", color="darkgreen", linewidth=2, alpha=0.6)


ax.set_title("SIRS Model on Barabasi-Albert vs Mean-Field Approximation")
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
    


peak_index_sim = np.argmax(I_means) # Index of max infection 
peak_time_sim = report_times[peak_index_sim] # Time at max infection 
peak_infection_sim = I_means[peak_index_sim] # Max infection proportion 
print(f"Simulation Peak Infection: {peak_infection_sim:.4f} at time {peak_time_sim:.2f}")
final_infection_sim = I_means[-1]
print(f"Final Simulation Infection Proportion: {final_infection_sim:.4f}")


peak_index_ode = np.argmax(I_ode) # Index of max infection 
peak_time_ode = t_ode[peak_index_ode] # Time at max infection 
peak_infection_ode = I_ode[peak_index_ode] # Max infection proportion 
print(f"ODE Peak Infection: {peak_infection_ode:.4f} at time {peak_time_ode:.2f}")



#fig.savefig("SIRS-barabasi_albert.png", dpi=300, bbox_inches="tight")  # Save as Png
plt.show()


# ******** |STORING DATA TO CSV| ********

#data = {
#        "t": report_times,
#        "S_mean": S_means,
#        "I_mean": I_means,
#        "R_mean": R_means
#        }

#df_SIRS = pd.DataFrame(data)
#print(df)

#df_SIRS.to_csv('./SIRS.csv', sep=';', index=False) # if ';' doesn't work then ','