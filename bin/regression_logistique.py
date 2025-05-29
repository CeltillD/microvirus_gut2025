#!/usr/bin/env python3

import sys
import pandas as pd 
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.linear_model import LogisticRegressionCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

def cluster_and_classify(df: pd.DataFrame, bin_vec: np.ndarray, K: int, T: int, condition_str: str = "CD"):
    """
    Objectifs :
    -> clusterisation hiérarchique des données par condition avec K clusters
    -> séléction aléatoire de T clusters pour le training 
    -> sous-compartimenter les données training/test par condition
    -> définition des données d'entrées (X) et cibles (Y)
    -> calibrage pour régularisation du modèle (L1)
    -> application de la régression logistique
    -> récupération des facteurs (genres) et leur influence
    -> sauvegarde des performances du modèle sur l'itération courante
    Sorties : valeur de régréssion (C) ; faux-positifs (FPR) ; vrais-positifs (TPR); performance (AUC) ; df des 10 meilleurs genres [facteur/influence/condition/inf_abs] 
    """
    idx_HC = np.where(bin_vec == 0)[0]                                                      # création d'index pour les conditions
    idx_case = np.where(bin_vec == 1)[0]

    cluster_HC = AgglomerativeClustering(n_clusters=K).fit_predict(df.iloc[idx_HC])         # application du clustering
    cluster_case = AgglomerativeClustering(n_clusters=K).fit_predict(df.iloc[idx_case])
    df = df.copy()                                                                          # préparation des données
    df['cluster'] = -1
    df.iloc[idx_HC, df.columns.get_loc('cluster')] = cluster_HC
    df.iloc[idx_case, df.columns.get_loc('cluster')] = cluster_case

    train_clusters_HC = np.random.choice(np.unique(cluster_HC), size=T, replace=False)      # sélection aléatoire des clusters d'entrainemet
    train_clusters_case = np.random.choice(np.unique(cluster_case), size=T, replace=False)

    train_HC = df[(bin_vec == 0) & (df['cluster'].isin(train_clusters_HC))]                 # sous-set des clusters de training HC
    test_HC = df[(bin_vec == 0) & (~df['cluster'].isin(train_clusters_HC))]                 # autres clusters HC pour test
    train_case = df[(bin_vec == 1) & (df['cluster'].isin(train_clusters_case))]
    test_case = df[(bin_vec == 1) & (~df['cluster'].isin(train_clusters_case))]
    train_df = pd.concat([train_HC, train_case])                                            # concaténation des données de test et training 
    test_df = pd.concat([test_HC, test_case])

    # données Y : cibles 
    train_y = train_df.index.astype(str).str.extract(r'(HC|CD|UC)', expand=False).map(lambda x: 0 if x == 'HC' else 1)
    test_y = test_df.index.astype(str).str.extract(r'(HC|CD|UC)', expand=False).map(lambda x: 0 if x == 'HC' else 1)
    # données X : entrées
    X_train = train_df.drop(columns=['cluster'])
    X_test = test_df.drop(columns=['cluster'])                                  
    
    # on cherche ici parmi ces valeurs c, celle qui nous renvoie les meilleurs résultats de test / minimise les pertes

    model = LogisticRegressionCV(   # validation croisée
    Cs=10,                          # test sur 10 valeurs de C
    cv=10,                          # 10-fold
    penalty='l1',                   # L1 car on considère que peu de genres associés
    solver='liblinear',             # solver classique
    scoring='neg_log_loss',         # se base sur le moins de pertes
    max_iter=10000                  # max d'itérations pour la convergence des résultats
    )
    model.fit(X_train, train_y)  
    print("\nC :", model.C_[0])     # affichage pour verification perso
    best_C = model.C_[0]        
    
    # éxecution de la régréssion linéaire 
    model = LogisticRegression(
    penalty='l1',
    solver='liblinear',   
    max_iter=10000,
    C=best_C
    )

    model.fit(X_train, train_y)

    y_scores_train = model.predict_proba(X_train)[:, 1]
    fpr_train, tpr_train, _ = roc_curve(train_y, y_scores_train)
    auc_train = auc(fpr_train, tpr_train)

    y_scores = model.predict_proba(X_test)[:, 1]
    fpr, tpr, _ = roc_curve(test_y, y_scores)
    auc_test = auc(fpr, tpr)                         # récupération des vecteurs faux positifs et vrais positifs et score AUC test

    ###
    coef = model.coef_[0]
    non_zero_mask = coef != 0                       
    influential_df = pd.DataFrame({                 
        'Facteur': X_train.columns[non_zero_mask],
        'Influence': coef[non_zero_mask],
        'Condition': [condition_str if w > 0 else 'HC' for w in coef[non_zero_mask]],
        'Abs': np.abs(coef[non_zero_mask])
    }).sort_values(by='Abs', ascending=False).drop(columns='Abs') # récupération des 10 facteurs les plus influents dans cette itération

    return best_C, fpr, tpr, auc_test, auc_train, influential_df


def run_pipeline(input_file: str, condition_str: str = "CD", K: int = 300, T: int = 210 ,Niter: int = 1):
    """
    Objectifs : 
    -> formater et séléctionner les données à comparer
    -> réaliser une boucle de N itérations de la régréssion logistique
    -> compiler et sauvegarder les informations résultants des itérations
    -> Affichage des performances du modèle (Courbes ROC)
    -> Affichage des influences différentielles des genres les plus réccurents et influents
    Sortie : -
    """

    glob_df = pd.read_csv(input_file, sep="\t", index_col=0)    # datas déja normalisées RPKM + CLR par virome
    t_df = glob_df.T                                            # avoir les prédicteurs en colonnes (genres) et viromes en lignes

    if condition_str == "IBD":                                  # IBD = UC + CD
        conditions = ['HC', 'CD', 'UC']
    else:
        conditions = ['HC', condition_str]

    filtered_df = t_df[t_df.index.astype(str).str.contains('|'.join(conditions))] # garde les viromes HC et (CD/UC)

    bin_vec = filtered_df.index.astype(str).str.extract(r'(HC|CD|UC)', expand=False)
    bin_vec = bin_vec.map(lambda x: 0 if x == 'HC' else 1).to_numpy()
    filtered_df['cluster'] = -1             # init colonnes cluster 
    results = []                            # liste les infos de chaque itération (AUC,C,K,T)
    influential_summary_list = []           # liste les 10 facteurs les plus influents par iteration (liste de dfs)

    plt.figure(figsize=(8, 6))              # plot courbes ROC

    ### Démarrage de N itérations de clustering et regression logistique :

    false_iter = 0

    for i in range(1, Niter+1):

        #####################################--- Barre progressive shell ---#####
        progress = int((i / Niter) * 100)                                       #
        bar = "[" + "=" * (progress - 1) + ">" + " " * (100 - progress) + "]"   #
        print(f"\r{bar} Iteration {i}/{Niter}", end="", flush=True)             #
        #########################################################################

        # appel de la fonction d'execution pour la regression logistique
        best_C, fpr, tpr, auc_test, auc_train, influential_df = cluster_and_classify(filtered_df, bin_vec, K, T, condition_str=condition_str)

        if auc_test != 0 and auc_train != 0 :
            # je sauvegarde les données des facteurs les plus influents sur l'itération courante
            top_factors_df = influential_df.copy()
            top_factors_df["Iteration"] = i
            influential_summary_list.append(
                top_factors_df
            )

            # meme chose pour les performances/infos de cette itération
            results.append({
                "Iteration": i,
                "condition": condition_str,
                "K": K,
                "T": T,
                "AUC_TEST": auc_test,
                "AUC_TRAIN": auc_train,
                "C": best_C
            })
            
            # tracé courbe ROC pour l'itération courante
            plt.plot(fpr, tpr, lw=1.5)

    # Courbes ROC global
    plt.plot([0, 1], [0, 1], 'k--', lw=1)
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.title(f"Courbes ROC : {Niter-false_iter} itérations - HC vs {condition_str}")
    plt.legend(loc="lower right", fontsize="small")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"ROC_{Niter-false_iter}iter_{condition_str}.svg", format="svg")
    plt.show()

    # concaténation des dpnnées des différentes itérations
    summary_df = pd.DataFrame(results)
    influential_summary_df = pd.concat(influential_summary_list, ignore_index=True)

    # seuil minimal d'occurence des facteurs de forte influence
    min_count = int(0.75 * Niter)
    # comptage des occurences
    factor_counts = influential_summary_df.groupby('Facteur')['Iteration'].nunique()
    # on garde les facteurs les plus réccurents (et influents) dans les N iterations
    frequent_factors = factor_counts[factor_counts >= min_count].index
    filtered_df = influential_summary_df[influential_summary_df['Facteur'].isin(frequent_factors)]

    # calcul d'un influence "moyenne" sur les N itérations 
    agg_df = (
        filtered_df
        .groupby("Facteur")["Influence"]
        .agg(['sum', 'std'])
        .reset_index()
        .assign(
            Influence_moyenne=lambda df: df["sum"] / Niter,
            Influence_std=lambda df: df["std"] / np.sqrt(Niter)  # erreur standard
        )
        .drop(columns="sum")
        .sort_values(by="Influence_moyenne", ascending=False)
    )

    mean_influence_df_sorted = agg_df.sort_values(by="Influence_moyenne")

    # Plot avec barres d'erreur
    plt.figure(figsize=(10, 8))
    colors = mean_influence_df_sorted["Influence_moyenne"].apply(lambda x: 'red' if x > 0 else 'blue')

    plt.barh(
        mean_influence_df_sorted["Facteur"],
        mean_influence_df_sorted["Influence_moyenne"],
        xerr=mean_influence_df_sorted["Influence_std"],
        color=colors,
        ecolor='black',     
        capsize=4           
    )

    plt.axvline(0, color='black', linewidth=1)
    plt.xlabel("Influence moyenne")
    plt.title(f"Influence différentielle des genres de microvirus (HC vs {condition_str})")
    plt.tight_layout()
    plt.text(-0.05, len(mean_influence_df_sorted), "HC", fontsize=12, color='blue', ha='center')
    plt.text(0.05, len(mean_influence_df_sorted), condition_str, fontsize=12, color='red', ha='center')
    plt.grid(True, axis='x', linestyle='--', alpha=0.5)
    #plt.savefig(f"FacteurInfluence_{Niter-false_iter}iter_{condition_str}.svg", format="svg")
    plt.show()

    # sauvegarde des résultats
    summary_df.to_csv(f"RegLineaire_{condition_str}-{Niter-false_iter}.tsv", sep="\t", index=False)
    influential_summary_df.to_csv(f"FacteurInfluence_{condition_str}-{Niter-false_iter}.tsv", sep="\t", index=False)
   

if __name__ == "__main__":
    if len(sys.argv) < 6 :
        print("Usage: ./regression_logistique.py < data_file.tsv > < Condition = UC/CD > < K = total_clusters > < T = training_clusters > < niter = iterations > < C = coeff-reg >")
        sys.exit(1)
    run_pipeline(sys.argv[1], condition_str=sys.argv[2], K=int(sys.argv[3]), T=int(sys.argv[4]), Niter=int(sys.argv[5]))
