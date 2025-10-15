import numpy as np
import networkx as nx
import EoN
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd


# Parameters
N = 10000  # network size
m=3 # nodes
p = 0.01  # Probability of edge creation
#k = p * (N-1)
beta = 0.3  # Infection rate
gamma = 0.1  # Recovery rate
t_max = 200  # Simulation duration
num_simulations = 10  # Number of simulations to average
report_times = np.linspace(0, t_max, 100)

#G = nx.grid_2d_graph(100, 100, periodic=True)
G = nx.fast_gnp_random_graph(N, p)
#G = nx.watts_strogatz_graph(N, 10, 0.1)  # N nodes, each connected to 10 neighbors, 10% rewiring
#G = nx.barabasi_albert_graph(N, m)


avg_degree = sum(dict(G.degree()).values()) / N
τ = beta / avg_degree


#  ODEs 
def SIR_ODE(t, y, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]

# Solve the ODE model
S0, I0, R0 = 0.99, 0.01, 0.0  # Initial conditions
t_eval = report_times
ode_solution = solve_ivp(SIR_ODE, [0, t_max], [S0, I0, R0], args=(beta, gamma), t_eval=t_eval) #Runge-Kutta επιλυση
t_ode = ode_solution.t
S_ode, I_ode, R_ode = ode_solution.y


S_means = np.zeros(len(report_times))
I_means = np.zeros(len(report_times))
R_means = np.zeros(len(report_times))


fig, ax = plt.subplots(figsize=(14, 8))


for sim in range(num_simulations):
   
    node_list = list(G.nodes)
    # Initialize the infected nodes
    initial_infected = np.random.choice(len(node_list), size=int(I0 * N), replace=False)
    initial_infected = [node_list[i] for i in initial_infected]
    
    # Run the SIR simulation
    t, S, I, R = EoN.fast_SIR(G, tau=τ, gamma=gamma, initial_infecteds=initial_infected, tmax=t_max)
    
    
    S_subsampled, I_subsampled, R_subsampled = EoN.subsample(report_times, t, S, I, R)
    ax.plot(report_times, S_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    ax.plot(report_times, I_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    ax.plot(report_times, R_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    S_means += S_subsampled / N  # Normalize by population size
    I_means += I_subsampled / N
    R_means += R_subsampled / N

# Average the simulation results
S_means /= num_simulations
I_means /= num_simulations
R_means /= num_simulations

# Plot
ax.scatter(report_times, S_means, color='blue', label="Simulation (S)", s=7)
ax.scatter(report_times, I_means, color='red', label="Simulation (I)", s=7)
ax.scatter(report_times, R_means, color='green', label="Simulation (R)", s=7)
ax.plot(t_ode, S_ode, label="ODE (S)", color="navy", linewidth=2, alpha=0.6)
ax.plot(t_ode, I_ode, label="ODE (I)", color="darkred", linewidth=2, alpha=0.6)
ax.plot(t_ode, R_ode, label="ODE (R)", color="darkgreen", linewidth=2, alpha=0.6)

ax.set_title("SIR Model on Barabasi-Albert vs Mean-Field Approximation")
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



peak_index_ode = np.argmax(I_ode) # Index of max infection 
peak_time_ode = t_ode[peak_index_ode] # Time at max infection 
peak_infection_ode = I_ode[peak_index_ode] # Max infection proportion 
print(f"ODE Peak Infection: {peak_infection_ode:.4f} at time {peak_time_ode:.2f}")


#fig.savefig("SIR-barabasi_albert.png", dpi=300, bbox_inches="tight")  # Save as Png

plt.show()


# ******** |STORING DATA TO CSV| ********

#data = {
#        "t": report_times,
#        "S_mean": S_means,
#        "I_mean": I_means,
#        "R_mean": R_means
#        }

#df_SIR = pd.DataFrame(data)
#print(df)

#df_SIR.to_csv('./SIR.csv', sep=';', index=False) # if ';' doesn't work then ','