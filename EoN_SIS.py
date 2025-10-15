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
#k = p * (N-1)
beta = 0.3  # Infection rate
#τ = beta / k
gamma = 0.1  # loss immunity
t_max = 100  # Simulation duration
num_simulations = 10  # Number of simulations to average
report_times = np.linspace(0, t_max, 100)

#G = nx.grid_2d_graph(100, 100, periodic=True)
G = nx.fast_gnp_random_graph(N, p)
#G = nx.watts_strogatz_graph(N, 10, 0.1)  # N nodes, each connected to 10 neighbors, 10% rewiring
#G = nx.barabasi_albert_graph(N, m)


avg_degree = sum(dict(G.degree()).values()) / N
τ = beta / avg_degree


#  ODEs 
def SIS_ODE(t, y, beta, gamma):
    S, I = y
    dSdt = gamma * I - beta * S * I
    dIdt = beta * S * I - gamma * I
    return [dSdt, dIdt]

# Solve the ODE model
S0, I0 = 0.99, 0.01  # Initial conditions
t_eval = report_times
ode_solution = solve_ivp(SIS_ODE, [0, t_max], [S0, I0], args=(beta, gamma), t_eval=t_eval) #Runge-Kutta επιλυση
t_ode = ode_solution.t
S_ode, I_ode = ode_solution.y


S_means = np.zeros(len(report_times))
I_means = np.zeros(len(report_times))

fig, ax = plt.subplots(figsize=(14, 8))

for sim in range(num_simulations):
    
    # Initialize the infected 1% nodes
    node_list = list(G.nodes)
    initial_infected = np.random.choice(len(node_list), size=int(I0 * N), replace=False)
    initial_infected = [node_list[i] for i in initial_infected]
    
    # Run the SIS simulation
    t, S, I = EoN.fast_SIS(G, tau=τ, gamma=gamma, initial_infecteds=initial_infected, tmax=t_max)
    
    
    S_subsampled, I_subsampled = EoN.subsample(report_times, t, S, I)
    ax.plot(report_times, S_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    ax.plot(report_times, I_subsampled / N, color='gray', alpha=0.4, linewidth=2)
    S_means += S_subsampled / N  # Normalize by population size
    I_means += I_subsampled / N

# Average the simulation results
S_means /= num_simulations
I_means /= num_simulations




# Plot
ax.scatter(report_times, S_means, label="Simulation (S)", color="blue", linewidth=2, s=7)
ax.scatter(report_times, I_means, label="Simulation (I)", color="red", linewidth=2, s=7)

ax.plot(t_ode, S_ode, label="ODE (S)", color="navy", linewidth=2, alpha=0.6)
ax.plot(t_ode, I_ode, label="ODE (I)", color="darkred", linewidth=2, alpha=0.6)

ax.set_title("SIS Model on Erdos-Renyi vs Mean-Field Approximation")
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






#fig.savefig("SIS-BarabassiAlbert.png", dpi=300, bbox_inches="tight")  # Save as Png    
plt.show()


# ******** |STORING DATA TO CSV| ********

#data = {
#        "t": report_times,
#        "S_mean": S_means,
#        "I_mean": I_means
#        }

#df_SIS = pd.DataFrame(data)
#print(df)

#df_SIS.to_csv('./SIS.csv', sep=';', index=False) # if ';' doesn't work then ','
