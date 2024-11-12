#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 12:28:22 2024

@author: maria
"""

## Calcular distribucion de grados

import networkx as nx
import powerlaw
import pandas as pd
#import matplotlib.pyplot as plt

# 1. Crear un grafo de ejemplo (puedes reemplazar esto con tu grafo)

def cal_gamma(file):

    df = pd.read_csv(file, delimiter=',')
    G = nx.Graph()


    for index, row in df.iterrows():
        G.add_edge(row['interaction'], row['name'])

    # 2. Calcular la distribución de grados
    degrees = [degree for node, degree in G.degree()]

    # 3. Usar la biblioteca powerlaw para ajustar la ley de potencia
    fit = powerlaw.Fit(degrees, discrete=True)  # Ajuste de la ley de potencia

    # Obtener el valor de gamma
    gamma = fit.power_law.alpha
    print(file)
    print(f'El valor de gamma es: {gamma:.4f}')

# 4. Visualizar la distribución y el ajuste de la ley de potencia
# fig = fit.plot_pdf(color='b', linewidth=2, label='Datos')  # Plot de los datos de la distribución de grados
# fit.power_law.plot_pdf(color='r', linestyle='--', ax=fig, label=f'Ajuste de Ley de Potencia\n$\gamma = {gamma:.2f}$')
# plt.xlabel('Grado')
# plt.ylabel('P(k)')
# plt.legend()
# plt.title('Distribución de Grados y Ajuste de Ley de Potencia')
# plt.show()

files = [
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_MA01_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_MA02_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_MA03_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_MA04_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_MA05_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_MA06_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_MA07_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_MA08_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_MA09_softpower40_3_signed_8_090_allGenes.csv default edge.csv",

"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_RS01_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_RS02_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_RS03_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_RS04_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_RS05_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_RS06_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_RS07_softpower40_3_signed_8_090_allGenes.csv default edge.csv",
"/home/mauricio/tfm-viu/Tfm_R_scripts/scripts_python/WGCNA_RS08_softpower40_3_signed_8_090_allGenes.csv default edge.csv"
]


f = lambda x : cal_gamma(x)
[cal_gamma(file) for file in files ]
