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
sigma = 0.2  # Rate of progression from exposed to infectious
omega = 0.05  # Rate of immunity loss (S <- R transition)
t_max = 180 # Simulation duration
num_simulations = 10  # Number of simulations to average
report_times = np.linspace(0, t_max, 100)


#G = nx.grid_2d_graph(100, 100, periodic=True)
#G = nx.fast_gnp_random_graph(N, p)
#G = nx.watts_strogatz_graph(N, 10, 0.1)  # N nodes, each connected to 10 neighbors, 10% rewiring
G = nx.barabasi_albert_graph(N, m)


avg_degree = sum(dict(G.degree()).values()) / N
τ = beta / avg_degree


# ODEs 
def SEIRS_ODE(t, y, beta, sigma, gamma, omega):
    S, E, I, R = y
    dSdt = -beta * S * I + omega * R
    dEdt = beta * S * I - sigma * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I - omega * R
    return [dSdt, dEdt, dIdt, dRdt]

# Initialize arrays for storing simulation results
S_means = np.zeros(len(report_times))
E_means = np.zeros(len(report_times))
I_means = np.zeros(len(report_times))
R_means = np.zeros(len(report_times))

# Initial conditions
S0, E0, I0, R0 = 0.99, 0.01, 0.00, 0.0

# Solve the ODE model
ode_solution = solve_ivp(
    SEIRS_ODE, 
    [0, t_max], 
    [S0, E0, I0, R0], 
    args=(beta, sigma, gamma, omega), 
    t_eval=report_times)  #Runge-Kutta επιλυση
t_ode = ode_solution.t
S_ode, E_ode, I_ode, R_ode = ode_solution.y

# Create transition graphs for SEIRS model
# Spontaneous transitions (E->I, I->R, R->S)
H = nx.DiGraph()
H.add_node('S')  # Need to add isolated nodes explicitly
H.add_edge('E', 'I', rate=sigma)  # E -> I transition
H.add_edge('I', 'R', rate=gamma)  # I -> R transition
H.add_edge('R', 'S', rate=omega)    # R -> S transition

# Neighbor-induced transitions (I+S->I+E)
J = nx.DiGraph()
J.add_edge(('I', 'S'), ('I', 'E'), rate=τ)



fig, ax = plt.subplots(figsize=(14, 8))


for sim in range(num_simulations):
    
    # Initialize node states
    IC = {}
    nodes = list(G.nodes())
    np.random.shuffle(nodes)
    
    n_exposed = int(E0 * N)
    n_infected = int(I0 * N)
    n_recovered = int(R0 * N)
    
    for i, node in enumerate(nodes):
        if i < n_exposed:
            IC[node] = 'E'
        elif i < n_exposed + n_infected:
            IC[node] = 'I'
        elif i < n_exposed + n_infected + n_recovered:
            IC[node] = 'R'
        else:
            IC[node] = 'S'

    # Run the Gillespie simulation
    t, S, E, I, R = EoN.Gillespie_simple_contagion(
        G, 
        H,  # spontaneous transition graph
        J,  # neighbor-induced transition graph
        IC, # initial conditions
        return_statuses=('S', 'E', 'I', 'R'),
        tmax=t_max
    )

    # Subsample the results
    S_subsampled = np.interp(report_times, t, S)
    E_subsampled = np.interp(report_times, t, E)
    I_subsampled = np.interp(report_times, t, I)
    R_subsampled = np.interp(report_times, t, R)
    
    ax.plot(report_times, S_subsampled / N, color='gray', alpha=0.3, linewidth=2)
    ax.plot(report_times, E_subsampled / N, color='gray', alpha=0.3, linewidth=2)
    ax.plot(report_times, I_subsampled / N, color='gray', alpha=0.3, linewidth=2)
    ax.plot(report_times, R_subsampled / N, color='gray', alpha=0.3, linewidth=2)

    
    # Add to means
    S_means += S_subsampled / N
    E_means += E_subsampled / N
    I_means += I_subsampled / N
    R_means += R_subsampled / N

# Average the simulation results
S_means /= num_simulations
E_means /= num_simulations
I_means /= num_simulations
R_means /= num_simulations




# Plot

ax.scatter(report_times, S_means, label="Simulation (S)", color="blue", linewidth=2, s=7)
ax.scatter(report_times, E_means, label="Simulation (E)", color="orange", linewidth=2, s=7)
ax.scatter(report_times, I_means, label="Simulation (I)", color="red", linewidth=2, s=7)
ax.scatter(report_times, R_means, label="Simulation (R)", color="green", linewidth=2, s=7)

ax.plot(t_ode, S_ode, label="ODE (S)", color="navy", linewidth=2, alpha=0.6)
ax.plot(t_ode, E_ode, label="ODE (E)", color="darkorange", linewidth=2, alpha=0.6)
ax.plot(t_ode, I_ode, label="ODE (I)", color="darkred", linewidth=2, alpha=0.6)
ax.plot(t_ode, R_ode, label="ODE (R)", color="darkgreen", linewidth=2, alpha=0.6)


ax.set_title("SEIRS Model on Barabasi-Albert vs Mean-Field Approximation")
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


peak_index_sim = np.argmax(E_means) # Index of max infection 
peak_time_sim = report_times[peak_index_sim] # Time at max infection 
peak_exposed_sim = E_means[peak_index_sim] # Max infection proportion 
print(f"Simulation Peak Exposed: {peak_exposed_sim:.4f} at time {peak_time_sim:.2f}")




peak_index_ode = np.argmax(I_ode) # Index of max infection 
peak_time_ode = t_ode[peak_index_ode] # Time at max infection 
peak_infection_ode = I_ode[peak_index_ode] # Max infection proportion 
print(f"ODE Peak Infection: {peak_infection_ode:.4f} at time {peak_time_ode:.2f}")


peak_index_sim = np.argmax(E_ode) # Index of max infection 
peak_time_sim = report_times[peak_index_sim] # Time at max infection 
peak_exposed_sim = E_ode[peak_index_sim] # Max infection proportion 
print(f"ODE Peak Exposed: {peak_exposed_sim:.4f} at time {peak_time_sim:.2f}")
    


final_exposed_sim = E_means[-1]
print(f"Final simulation Exposed Proportion: {final_exposed_sim: .4f}")
final_infection_sim = I_means[-1]
print(f"Final Simulation Infection Proportion: {final_infection_sim:.4f}")
    
    
    
    
#fig.savefig("SEIRS-barabasi_albert.png", dpi=300, bbox_inches="tight")  # Save as Png
plt.show()


#data = {
#        "t": report_times,
#        "S_mean": S_means,
#        "E_mean": E_means,
#        "I_mean": I_means,
#        "R_mean": R_means
#        }

#df_SEIRS = pd.DataFrame(data)
#print(df)

#df_SEIRS.to_csv('./SEIRS.csv', sep=';', index=False) # if ';' doesn't work then ','