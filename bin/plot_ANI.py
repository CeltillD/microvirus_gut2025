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
        if fraAliA >= 70.0 and fraAliB >= 70.0 and ANI >= 50.0 and FracAli <= 100 :
            xj.append(int(ANI))
            yj.append(int(fraAliA))
            zj.append(round(float(li[12]) * 100, 1))
data = pd.DataFrame({'ANI': x, 'FracAli': y, 'ratioLen': z})
contingency_table = pd.crosstab( data['FracAli'],data['ANI'])

data_i = pd.DataFrame({'ANI': xi, 'FracAli': yi, 'ratioLen': zi})
contingency_table_i = pd.crosstab(data_i['FracAli'],data_i['ANI'])

data_j = pd.DataFrame({'ANI': xj, 'FracAli': yj, 'ratioLen': zj})
contingency_table_j = pd.crosstab(data_j['FracAli'],data_j['ANI'])

print(contingency_table_j)


"""# Créer la figure
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
plt.show()"""




plt.figure(figsize=(10, 10))
ax = sns.heatmap(contingency_table_j, cmap='Blues',  square=True)
ax.invert_yaxis()
plt.title('number of occurrences')
plt.ylabel('FracAli')
plt.xlabel('ANI')
#plt.savefig("heatmap_clean_sans_ZUB.svg", format="svg")
plt.show()

# Définir les intervalles et les palettes de couleurs
intervals = [
    (50, 60, 'Reds'), 
    (60, 70, 'Oranges'), 
    (70, 80, 'Greens'), 
    (80, 90, 'Blues'), 
    (90, 99, 'Purples'), 
    (100, None, 'Greys') # None pour inclure uniquement la valeur unique de FracAli =100
]

# Créer le graphique
plt.figure(figsize=(10,6))

for start, end, cmap in intervals:
    # Filtrer les données pour l'intervalle donné
    if end is not None:
        filtered_data = contingency_table_j.loc[start:end-1]
        gradient = np.linspace(0.2, 1.0, len(filtered_data)) # Gradient clair à foncé
        for i, (index, row) in enumerate(filtered_data.iterrows()):
            color = plt.cm.get_cmap(cmap)(gradient[i])
            plt.plot(row.index.values.tolist(), row.values.tolist(), label=f"FracAli={index}", color=color)
    else:
        # Cas spécial pour FracAli=100
        row = contingency_table_j.loc[start]
        plt.plot(row.index.values.tolist(), row.values.tolist(), label="FracAli=100", color=plt.cm.get_cmap(cmap)(1.0))

# Ajouter des légendes et des labels
plt.title("Courbes avec dégradés de couleurs par intervalle")
plt.xlabel("ANI")
plt.ylabel("Valeurs")
plt.legend(loc='upper right', fontsize='small', ncol=2)
plt.grid(True)
plt.show()


intervals = [
    (50, 60, 'red'), 
    (60, 70, 'orange'), 
    (70, 80, 'green'), 
    (80, 90, 'blue'), 
    (90, None, 'purple'), # None pour inclure les valeurs jusqu'à FracAli=99
    (100, None, 'black') # FracAli=100
]

# Créer le graphique
plt.figure(figsize=(10,6))

# Tracer une courbe par intervalle
for start, end, color in intervals:
    if end is not None:
        # Filtrer les données pour l'intervalle donné
        filtered_data = contingency_table_j.loc[start:end-1]
        summed_values = filtered_data.sum(axis=0) # Somme des lignes pour chaque colonne ANI
        plt.plot(summed_values.index.values.tolist(), summed_values.values.tolist(), label=f"[{start}-{end}[", color=color)
    else:
        # Cas spécial pour FracAli=100 ou intervalle final
        if start == 100:
            row = contingency_table_j.loc[start]
            plt.plot(row.index.values.tolist(), row.values.tolist(), label="100", color=color)
        else:
            filtered_data = contingency_table_j.loc[start:]
            summed_values = filtered_data.sum(axis=0)
            plt.plot(summed_values.index.values.tolist(), summed_values.values.tolist(), label=f"[{start}-99]", color=color)

# Ajouter la courbe noire en pointillés pour la totalité des données
total_values = contingency_table_j.sum(axis=0) # Somme de toutes les lignes pour chaque colonne ANI
plt.plot(total_values.index.values.tolist(), total_values.values.tolist(), label="Total", color="black", linestyle="--")

# Ajouter des légendes et des labels
plt.title("Courbes par intervalle et total global")
plt.xlabel("ANI")
plt.ylabel("Valeurs")
plt.legend(loc='upper right', fontsize='small', ncol=2)
plt.grid(True)
#plt.savefig("courbes_clean_sans_ZUB.svg", format="svg")
plt.show()
"""
plt.figure(figsize=(10, 10))
ax = sns.heatmap(contingency_table_i, cmap='Reds',  square=True)
ax.invert_yaxis()
plt.title('number of occurrences')
plt.ylabel('FracAli')
plt.xlabel('ANI')
#plt.savefig("heatmap_ANI_Frac_occurences.svg", format="svg")
plt.show()

df_sum = contingency_table.sum(axis=0)
df_sum_i = contingency_table_i.sum(axis=0)
df_sum_j = contingency_table_j.sum(axis=0)

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

"""# Tracer les courbes sur le même plot
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
"""

"""# Tracer la courbe
plt.figure(figsize=(10, 6))
plt.plot(df_sum_j.index, df_sum_j.values, marker='o')
plt.xlabel('ANI')
plt.ylabel('Somme des lignes')
plt.title('Somme des lignes pour chaque ANI')
plt.grid(True)
plt.show()"""

intervals = [
    (70, 90, 'red'), 
    (90, 101, 'blue'), 
]

# Créer le graphique
plt.figure(figsize=(10,6))

# Tracer une courbe par intervalle
for start, end, color in intervals:
    if end is not None:
        # Filtrer les données pour l'intervalle donné
        filtered_data = contingency_table_j.loc[start:end-1]
        summed_values = filtered_data.sum(axis=0) # Somme des lignes pour chaque colonne ANI
        plt.plot(summed_values.index.values.tolist(), summed_values.values.tolist(), label=f"[{start}-{end}[", color=color)
    else:
        # Cas spécial pour FracAli=100 ou intervalle final
        if start == 100:
            row = contingency_table_j.loc[start]
            plt.plot(row.index.values.tolist(), row.values.tolist(), label="100", color=color)
        else:
            filtered_data = contingency_table_j.loc[start:]
            summed_values = filtered_data.sum(axis=0)
            plt.plot(summed_values.index.values.tolist(), summed_values.values.tolist(), label=f"[{start}-99]", color=color)

# Ajouter la courbe noire en pointillés pour la totalité des données
total_values = contingency_table_j.sum(axis=0) # Somme de toutes les lignes pour chaque colonne ANI
plt.plot(total_values.index.values.tolist(), total_values.values.tolist(), label="Total", color="black", linestyle="--")

# Ajouter des légendes et des labels
plt.title("Courbes par intervalle et total global")
plt.xlabel("ANI")
plt.ylabel("Valeurs")
plt.legend(loc='upper right', fontsize='small', ncol=2)
plt.grid(True)
plt.savefig("courbes_clean_avec_ZUB_2.svg", format="svg")
plt.show()
