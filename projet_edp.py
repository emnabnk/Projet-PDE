import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import time


# Paramètres

D = 0.5  # Coefficient de diffusion
C = 0  #coefficient de convection
a=0
b=1
n = 500  # Nombre de points
t_final = 0.5  # Temps final

# Conditions aux limites
u_L= 0
u_R = 0  
delta_x = (b - a) / (n - 1)
gamma = 0.1
nu = 0.5
# Conditions de stabilité
delta_t1 = gamma * delta_x**2 / D
#delta_t2 = nu * delta_x / C
delta_t = delta_t1 #np.min(delta_t1,delta_t2)   


def u_exacte(x, t):
    u_exacte=np.sin(np.pi * x) * (1 + t)
    return u_exacte


def differences_finies (x, delta_x, delta_t, t_final, D, u_L, u_R):
    n = len(x)
    u = np.sin(np.pi * x)  # quand t=0, u(x,0) = sin(pi * x)
    u[0] = u_L
    u[-1] = u_R
    # Nombre de pas de temps
    n_t = int(t_final / delta_t) + 1
    # Matrice pour stocker les solutions à chaque instant t
    u_all = np.zeros((n_t, n))
    u_all[0, :] = u
    for t in range(1, n_t):
        u_new = np.copy(u)
        for i in range(1, n - 1):
            u_new[i] = u[i] + D * delta_t / delta_x**2 * (u[i-1] - 2*u[i] + u[i+1])
        u_new[0]=u_L 
        u_new[-1] = u_R  # Conditions aux limites
        u = u_new
        u_all[t, :] = u
    return u_all, delta_t * np.arange(n_t)


# Création de la grille spatiale
x = np.linspace(a, b, n)

# Calcul numérique par différences finies
u_numeric, t_values = differences_finies(x, delta_x, delta_t, t_final, D, u_L, u_R)

# Solution exacte
u_exact_all = np.array([u_exacte(x, t) for t in t_values])

# Calcul de l'erreur
#error = np.abs(u_numeric - u_exact_all)
error = np.linalg.norm(u_numeric - u_exact_all)/np.linalg.norm(u_exact_all)

print(error)
print(u_numeric)
print(u_exact_all)

# # Graphique : Évolution de u(x, t)
# plt.figure(figsize=(8, 5))
# plt.plot(u_numeric, t_values, label = 'u_numeric')

# plt.plot(u_exact_all, t_values,label='u_exact_all')
# plt.xlabel("x")
# plt.ylabel("t")
# plt.title("Évolution de u(x, t) (numérique et excate)")
# plt.show()

# Plot : Comparaison temporelle
plt.figure(figsize=(10, 6))
for i, t in enumerate([0, len(t_values)//4, len(t_values)//2, -1]):
    plt.plot(x, u_numeric[t, :], label=f'Numérique t={t_values[t]:.2f}')
    plt.plot(x, u_exact_all[t, :], '--', label=f'Exacte t={t_values[t]:.2f}')

plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.title('Comparaison entre solution numérique et exacte à différents instants')
plt.legend()
plt.grid()
plt.show()

# Plot : Comparaison spatiale
plt.figure(figsize=(10, 6))
for x_pos in [0, n//4, n//2, -1]:
    plt.plot(t_values, u_numeric[:, x_pos], label=f'Numérique x={x[x_pos]:.2f}')
    plt.plot(t_values, u_exact_all[:, x_pos], '--', label=f'Exacte x={x[x_pos]:.2f}')

plt.xlabel('t')
plt.ylabel('u(x, t)')
plt.title('Comparaison entre solution numérique et exacte à différents points spatiaux')
plt.legend()
plt.grid()
plt.show()

# # Graphique : Erreur entre la solution numérique et exacte
# plt.figure(figsize=(8, 5))
# plt.imshow(error, extent=[a, b, 0, t_final], origin='lower', aspect='auto', cmap='inferno')
# plt.colorbar(label="Erreur |u_num - u_exact|")
# plt.xlabel("x")
# plt.ylabel("t")
# plt.title("Erreur entre solution numérique et exacte")
# plt.show()