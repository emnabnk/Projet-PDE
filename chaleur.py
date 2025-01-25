import numpy as np
import matplotlib.pyplot as plt

######## ETUDE DE L'EQUATION DE LA CHALEUR (C = 0) : nsin(2pix) * exp(-4 pi**2 * D * t ) ##########
# PARAMETRES 
a = 0
b = 1
D = 0.5  # coefficient de diffusion
L = 1
Time = 0.1
NX = 100  # Discrétisation espace
NT = 1000  # Discrétisation temps  

u_L = 0
u_R = 0

delta_x = L / (NX - 1)
delta_t = Time / NT

C = delta_t * D / (delta_x**2)

# Initialisation
x = np.linspace(a, b, NX)
T = np.sin(2 * np.pi * x)
u = np.zeros((NX))
# Solution exacte
Texact = np.sin(2 * np.pi * x) * np.exp(-4 * np.pi**2 * D * Time)

for n in range(0, NT):
    for j in range(1, NX - 1):
        u[j] = C * (T[j-1] - 2 * T[j] + T[j+1])

    for j in range(1, NX - 1):
        T[j] += u[j]

    # Plot every 100 time steps
    if (n % 100 == 0):
        plotlabel = "t = %1.2f" % (n * delta_t)
        couleurs = plt.cm.plasma(n / NT)
        plt.plot(x, T, label=plotlabel, color=couleurs)

plt.xlabel('x', fontsize=26)
plt.ylabel('T', fontsize=26, rotation=0)
plt.title('Equation de la chaleur 1D')
plt.legend()
plt.show()

# Tracé de la solution exacte (après la boucle)
plt.plot(x, Texact, label="Solution exacte", color="black", linestyle="--")
plt.xlabel('x', fontsize=26)
plt.ylabel('T', fontsize=26, rotation=0)
plt.title('Solution exacte de l\'équation de la chaleur 1D')
plt.legend()
plt.show()

# Calcul de l'erreur
error = np.linalg.norm(T - Texact) / np.linalg.norm(Texact)
print(f"Erreur : {error}")
