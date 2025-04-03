import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

x = []
y = []
z = []

xi = []
yi = []
zi = []

xj = []
yj = []
zj = []
with open(sys.argv[1], "r") as TsvFile:
    next(TsvFile)
    for line in TsvFile:
        li = line.strip().split("\t")
        
        # Calcul avec précision décimale
        fraAliB = round(float(li[11]) * 100, 1)
        fraAliA = round(float(li[10]) * 100, 1)
        FracAli = round(((fraAliB + fraAliA)/2),1)
        ANI = round(float(li[13]),1)
        
        if fraAliA >= 80.0 and fraAliB >= 80.0 and ANI >= 50.0 and FracAli < 100 :
            x.append(ANI)
            y.append(FracAli)
            z.append(round(float(li[12]) * 100, 1))
        elif fraAliA >= 80.0 and fraAliB >= 80.0 and ANI >= 50.0 and FracAli == 100 :
            xi.append(ANI)
            yi.append(FracAli)
            zi.append(round(float(li[12]) * 100, 1))
        if fraAliA >= 80.0 and fraAliB >= 80.0 and ANI >= 50.0 and FracAli <= 100 :
            xj.append(ANI)
            yj.append(FracAli)
            zj.append(round(float(li[12]) * 100, 1))
data = pd.DataFrame({'ANI': x, 'FracAli': y, 'ratioLen': z})
contingency_table = pd.crosstab( data['FracAli'],data['ANI'])

data_i = pd.DataFrame({'ANI': xi, 'FracAli': yi, 'ratioLen': zi})
contingency_table_i = pd.crosstab(data_i['FracAli'],data_i['ANI'])

data_j = pd.DataFrame({'ANI': xj, 'FracAli': yj, 'ratioLen': zj})
contingency_table_j = pd.crosstab(data_j['FracAli'],data_j['ANI'])

print(contingency_table_j)


# Créer la figure
fig = plt.figure(figsize=(10, 12))

# Gridspec avec des ratios adaptés (rouge plus petite)
gs = fig.add_gridspec(2, 1, height_ratios=[0.01, 1])  # Rouge au-dessus

# Créer le premier subplot pour la heatmap rouge (en haut)
ax1 = fig.add_subplot(gs[0, 0])
heatmap1 = sns.heatmap(
    contingency_table_i,
    ax=ax1,
    cmap='Reds',
    cbar=False,
    xticklabels=False,
    yticklabels=False,
)
ax1.invert_yaxis()  # Inverser l'axe des y

# Créer le second subplot pour la heatmap bleue (en bas)
ax2 = fig.add_subplot(gs[1, 0])
heatmap2 = sns.heatmap(
    contingency_table,
    ax=ax2,
    cmap='Blues',
    cbar=False,
    xticklabels=True,
    yticklabels=True,
)

ax2.invert_yaxis()  # Inverser l'axe des y

# Ajouter des colorbars de même taille
cbar_ax = fig.add_axes([0.92, 0.05, 0.02, 0.4])  # Position pour la barre bleue
cbar_ax2 = fig.add_axes([0.92, 0.55, 0.02, 0.4])  # Position pour la barre rouge

plt.colorbar(heatmap2.collections[0], cax=cbar_ax, label="nb occurences (COV < 100%)")
plt.colorbar(heatmap1.collections[0], cax=cbar_ax2, label="nb occurences (COV = 100%)")

# Ajuster l'espacement pour qu'il n'y ait pas de séparation entre les heatmaps
plt.subplots_adjust(hspace=0)
#plt.savefig("heatmap_ANI_Frac_occurences.svg", format="svg")
# Afficher les graphiques
plt.show()




"""plt.figure(figsize=(10, 10))
ax = sns.heatmap(contingency_table, cmap='Blues',  square=True)
ax.invert_yaxis()
plt.title('number of occurrences')
plt.ylabel('FracAli')
plt.xlabel('ANI')
#plt.savefig("heatmap_ANI_Frac_occurences.svg", format="svg")
plt.show()

plt.figure(figsize=(10, 10))
ax = sns.heatmap(contingency_table_i, cmap='Reds',  square=True)
ax.invert_yaxis()
plt.title('number of occurrences')
plt.ylabel('FracAli')
plt.xlabel('ANI')
#plt.savefig("heatmap_ANI_Frac_occurences.svg", format="svg")
plt.show()
"""
df_sum = contingency_table.sum(axis=0)
df_sum_i = contingency_table_i.sum(axis=0)
df_sum_j = contingency_table_j.sum(axis=0)
"""
# Tracer la courbe bleu foncé
plt.figure(figsize=(10, 6))
plt.plot(df_sum.index, df_sum.values, marker='o')
plt.xlabel('ANI')
plt.ylabel('Somme des lignes')
plt.title('Somme des lignes pour chaque ANI')
plt.grid(True)
plt.show()

# Tracer la courbe en rouge
plt.figure(figsize=(10, 6))
plt.plot(df_sum_i.index, df_sum_i.values, marker='o')
plt.xlabel('ANI')
plt.ylabel('Somme des lignes')
plt.title('Somme des lignes pour chaque ANI')
plt.grid(True)
plt.show()
"""

# Tracer les courbes sur le même plot
plt.figure(figsize=(10, 6))
plt.plot(df_sum.index, df_sum.values, marker='', label='COV < 100%')
plt.plot(df_sum_i.index, df_sum_i.values, marker='', color='red', label='COV = 100%')
plt.plot(df_sum_j.index, df_sum_j.values, marker='', color='black', label='ALL')
# Ajouter des labels et un titre
plt.xlabel('ANI')
plt.ylabel('Nombre d\'occurences')
plt.title('')

# Ajouter une légende
plt.legend()

# Afficher la grille
plt.grid(True)
#plt.savefig("courbes_ANI_occurences.svg", format="svg")
# Afficher le plot
plt.show()


"""# Tracer la courbe
plt.figure(figsize=(10, 6))
plt.plot(df_sum_j.index, df_sum_j.values, marker='o')
plt.xlabel('ANI')
plt.ylabel('Somme des lignes')
plt.title('Somme des lignes pour chaque ANI')
plt.grid(True)
plt.show()"""
