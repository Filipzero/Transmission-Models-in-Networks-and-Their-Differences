import numpy as np
import networkx as nx
import EoN
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd

# Parameters
N = 10000  # network size 
p = 0.01  # Probability of edge creation in an Erdős–Rényi network
m = 3  # BA model parameter
beta = 0.3  # Infection rate
gamma = 0.1  # Recovery rate
sigma = 0.2  # Exposed -> Infectious rate
t_max = 180
num_simulations = 10
report_times = np.linspace(0, t_max, 100)

# Network
G = nx.barabasi_albert_graph(N, m)
#G = nx.grid_2d_graph(100, 100, periodic=True)
#G = nx.fast_gnp_random_graph(N, p)
#G = nx.watts_strogatz_graph(N, 10, 0.1)  # N nodes, each connected to 10 neighbors, 10% rewiring

avg_degree = sum(dict(G.degree()).values()) / N
τ = beta / avg_degree

# SEIR ODEs
def SEIR_ODE(t, y, beta, sigma, gamma):
    S, E, I, R = y
    dSdt = -beta * S * I
    dEdt = beta * S * I - sigma * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I
    return [dSdt, dEdt, dIdt, dRdt]

# Initialize arrays
S_means = np.zeros(len(report_times))
E_means = np.zeros(len(report_times))
I_means = np.zeros(len(report_times))
R_means = np.zeros(len(report_times))

# Initial conditions
S0, E0, I0, R0 = 0.99, 0.01, 0.00, 0.0

# Solve ODE
ode_solution = solve_ivp(
    SEIR_ODE,
    [0, t_max],
    [S0, E0, I0, R0],
    args=(beta, sigma, gamma),
    t_eval=report_times
)
t_ode = ode_solution.t
S_ode, E_ode, I_ode, R_ode = ode_solution.y

# Transition graphs
H = nx.DiGraph()
H.add_node('S')
H.add_edge('E', 'I', rate=sigma)
H.add_edge('I', 'R', rate=gamma)

J = nx.DiGraph()
J.add_edge(('I', 'S'), ('I', 'E'), rate=τ)  # S->E infection

fig, ax = plt.subplots(figsize=(14, 8))

for sim in range(num_simulations):
    # Initialize node states
    IC = {}
    nodes = list(G.nodes())
    np.random.shuffle(nodes)
    
    n_infected = int(I0 * N)
    n_exposed = int(E0 * N)
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
    
    # Gillespie simulation
    t, S, E, I, R = EoN.Gillespie_simple_contagion(
        G,
        H,
        J,
        IC,
        return_statuses=('S', 'E', 'I', 'R'),
        tmax=t_max
    )

    # Use np.interp for subsampling
    S_subsampled = np.interp(report_times, t, S)
    E_subsampled = np.interp(report_times, t, E)
    I_subsampled = np.interp(report_times, t, I)
    R_subsampled = np.interp(report_times, t, R)

    ax.plot(report_times, S_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    ax.plot(report_times, E_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    ax.plot(report_times, I_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    ax.plot(report_times, R_subsampled / N, color='gray', alpha=0.4, linewidth=2)

    # Add to means
    S_means += S_subsampled / N
    E_means += E_subsampled / N
    I_means += I_subsampled / N
    R_means += R_subsampled / N

# Average
S_means /= num_simulations
E_means /= num_simulations
I_means /= num_simulations
R_means /= num_simulations

# Plot averages
ax.scatter(report_times, S_means, label="Simulation (S)", color="blue", s=7)
ax.scatter(report_times, E_means, label="Simulation (E)", color="orange", s=7)
ax.scatter(report_times, I_means, label="Simulation (I)", color="red", s=7)
ax.scatter(report_times, R_means, label="Simulation (R)", color="green", s=7)

ax.plot(t_ode, S_ode, label="ODE (S)", color="navy", linewidth=2, alpha=0.6)
ax.plot(t_ode, E_ode, label="ODE (E)", color="darkorange", linewidth=2, alpha=0.6)
ax.plot(t_ode, I_ode, label="ODE (I)", color="darkred", linewidth=2, alpha=0.6)
ax.plot(t_ode, R_ode, label="ODE (R)", color="darkgreen", linewidth=2, alpha=0.6)

ax.set_title("SEIR Model on Barabasi-Albert vs Mean-Field Approximation")
ax.set_xlabel("Time (t)", fontsize=14, labelpad=15, color="#333333")
ax.set_ylabel("Proportion of Population", fontsize=14, labelpad=15, color="#333333")
ax.tick_params(axis='both', which='major', labelsize=12, colors="#333333")
ax.legend(fontsize=14, loc="best", frameon=True, shadow=True)
ax.grid(color='lightgray', linestyle='-', linewidth=1.5)
ax.set_axisbelow(True)
for spine in ax.spines.values():
    spine.set_edgecolor('#bbbbbb')

# Peak infection
peak_index_sim = np.argmax(I_means)
peak_time_sim = report_times[peak_index_sim]
peak_infection_sim = I_means[peak_index_sim]
print(f"Simulation Peak Infection: {peak_infection_sim:.4f} at time {peak_time_sim:.2f}")
final_infection_sim = I_means[-1]
print(f"Final Simulation Infection Proportion: {final_infection_sim:.4f}")

peak_index_ode = np.argmax(I_ode)
peak_time_ode = t_ode[peak_index_ode]
peak_infection_ode = I_ode[peak_index_ode]
print(f"ODE Peak Infection: {peak_infection_ode:.4f} at time {peak_time_ode:.2f}")

plt.show()


# ******** |STORING DATA TO CSV| ********

#data = {
#        "t": report_times,
#        "S_mean": S_means,
#        "E_mean": E_means,
#        "I_mean": I_means,
#        "R_mean": R_means
#        }

#df_SEIR = pd.DataFrame(data)
#print(df)

#df_SEIR.to_csv('./SEIR.csv', sep=';', index=False) # if ';' doesn't work then ','