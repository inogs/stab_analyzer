import numpy as np
import matplotlib.pyplot as plt
import lyapunovV as lv

# Parametri
a = 2  # parametro della mappa espansiva
sigma = 1e-7  # parametro di scala per il rumore
delta0 = 0.0001  # soglia minima
DELTA = 0.001  # massimo delta
B = 0.00001  # Perturbazione iniziale

# Funzione per calcolare l'errore tra le traiettorie
def f_map(x, t, a, b, T):
    """ Alterna tra mappa espansiva e contrattiva """
    if t % (2 * T) < T:
        return a * x % 1  # Mappa espansiva
    else:
        return b * x  # Mappa contrattiva

# Funzione per calcolare tau
def calculate_tau(x, y, delta0):
    tau_values = []  # Lista per raccogliere i valori di tau
    i = 0
    while i < len(x) - 1:
        tau = 0
        # Se la distanza tra x[i] e y[i] è minore o uguale a delta0, continua
        while i < len(x) - 1 and np.abs(x[i] - y[i]) <= delta0:
            tau += 1
            i += 1
        if tau > 0:
            tau_values.append(tau)  # Aggiungi il valore di tau valido
        else:
            i += 1  # Se non è valido, salta alla prossima coppia
    return tau_values

# Funzione per calcolare l'esponente di Lyapunov
def calculate_lyapunov(T, num_steps, a, b, sigma, delta0, DELTA, B):
    x_vals = [0.05]  # Condizione iniziale per la traiettoria originale
    y_vals = [x_vals[0] + B * x_vals[0]]  # Condizione iniziale per la traiettoria perturbata  
    for t in range(num_steps):
        # Evoluzione delle traiettorie
        x_new = f_map(x_vals[-1], t, a, b, T)
        y_new = f_map(y_vals[-1], t, a, b, T)

        # Aggiungi il rumore
        noise = np.random.uniform(-1/2, 1/2)  # Rumore uniforme
        y_new += sigma * noise  # Perturbazione della traiettoria

        x_vals.append(x_new)
        y_vals.append(y_new)

    # Calcola tau e l'esponente di Lyapunov
    tau_values = calculate_tau(x_vals, y_vals, delta0)
    if tau_values:
        mean_tau = np.mean(tau_values)
        lyap = (1 / mean_tau) * np.log(DELTA / delta0)
    else:
        lyap = None

    return lyap, x_vals, y_vals

# Parametri di esecuzione
num_steps = 70
#T_values = [5]
T_values = [2,5, 10, 20, 50, 80, 100, 200, 400, 600, 800, 1000]   # Ciclo su T (tempo per ciascuna mappa)

plt.figure(figsize=(10, 5))
# Calcolo di Lyapunov per b1 = 2/3
# Liste per i valori di Lyapunov e T
lyapunov_values_b1 = []  # Per b1 = 2/3
lyapunov_values_b2 = []  # Per b2 = 1/4
b1 = 2/3

for T in T_values:
    x_vals = [0.05] 
    for t in range(num_steps):
    # Evoluzione delle traiettorie
        x_new = f_map(x_vals[-1], t, a, b1, T)
        x_vals.append(x_new)
    x_vals = np.array(x_vals)
    lyap = lv.LYAP(x_vals)
    lyapunov_b1 = lyap.lyap_e_paladin(dt=1,delta0=delta0, Delta=DELTA, ndim=1, tau=1,ires=5)[-1]
    lyapunov_values_b1.append(lyapunov_b1)

# Calcolo di Lyapunov per b2 = 1/4
b2 = 1/4
for T in T_values:
    x_vals = [0.05] 
    for t in range(num_steps):
    # Evoluzione delle traiettorie
        x_new = f_map(x_vals[-1], t, a, b2, T)
        x_vals.append(x_new)
    x_vals = np.array(x_vals)
    lyapunov_b2 = lyap.lyap_e_paladin(dt=1,delta0=delta0, Delta=DELTA, ndim=1, tau=1, ires=5)[-1]
    lyapunov_values_b2.append(lyapunov_b2)

# Plot con quadrati per b1 = 2/3
plt.scatter(T_values, lyapunov_values_b1, marker='s', edgecolor='black', facecolors='none', label=f"b = 2/3")

# Plot con rombi per b2 = 1/4
plt.scatter(T_values, lyapunov_values_b2, marker='D', edgecolor='red', facecolors='none', label="b = 1/4")

# Impostazioni del grafico
plt.xscale('log')  # Scala logaritmica per l'asse x
plt.ylim(0, 0.4)  # Limiti dell'asse y
plt.xlim(1, 1000)  # Limiti dell'asse x

# Etichette e titolo
plt.xlabel('T (Time Interval)', fontsize=12)
plt.ylabel('Lyapunov Exponent', fontsize=12)
plt.title('Lyapunov Exponent vs T per due valori di b', fontsize=14)
plt.grid(True)
plt.savefig('prova.png')
