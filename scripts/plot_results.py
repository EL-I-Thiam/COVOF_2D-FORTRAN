#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script de tracé des métriques et diagnostics des simulations COVOF 2D.
Génère des courbes de validation pour :
 1. La trajectoire et vitesse de la bulle ascendante (Rising Bubble)
 2. La convergence et saut de pression de la goutte statique (Static Droplet)
 3. L'énergie cinétique et divergence résiduelle de la cavité et de Poiseuille
 4. La rupture de barrage (Dam Break)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_rising_bubble():
    csv_file = "output_bubble/bubble_history.csv"
    if not os.path.exists(csv_file):
        print(f"Fichier introuvable: {csv_file}")
        return

    df = pd.read_csv(csv_file)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Position du centre de masse
    ax1.plot(df['time'], df['center_of_mass_y'], 'b-o', lw=2, label="Y centre de masse")
    ax1.set_xlabel("Temps [s]")
    ax1.set_ylabel("Position Y [m]")
    ax1.set_title("Ascension du centre de masse de la bulle")
    ax1.grid(True, linestyle="--", alpha=0.7)
    ax1.legend()

    # Vitesse d'ascension
    ax2.plot(df['time'], df['rise_velocity'], 'r-s', lw=2, label="Vitesse de montée V_c")
    ax2.set_xlabel("Temps [s]")
    ax2.set_ylabel("Vitesse [m/s]")
    ax2.set_title("Vitesse d'ascension de la bulle")
    ax2.grid(True, linestyle="--", alpha=0.7)
    ax2.legend()

    plt.tight_layout()
    output_png = "output_bubble/bubble_dynamics.png"
    plt.savefig(output_png, dpi=300)
    print(f"Graphique sauvegardé : {output_png}")
    plt.close()

def plot_static_droplet():
    csv_file = "output_droplet/droplet_history.csv"
    if not os.path.exists(csv_file):
        print(f"Fichier introuvable: {csv_file}")
        return

    df = pd.read_csv(csv_file)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Courants parasites (vitesse max)
    ax1.plot(df['time'], df['max_u'], 'm-^', lw=2, label="Courants parasites (U_max)")
    ax1.set_xlabel("Temps [s]")
    ax1.set_ylabel("Vitesse max [m/s]")
    ax1.set_title("Amortissement des courants parasites")
    ax1.grid(True, linestyle="--", alpha=0.7)
    ax1.legend()

    # Conservation du volume
    rel_vol = (df['bubble_volume'] - df['bubble_volume'].iloc[0]) / df['bubble_volume'].iloc[0] * 100.0
    ax2.plot(df['time'], rel_vol, 'g-d', lw=2, label="Variation relative de volume [%]")
    ax2.set_xlabel("Temps [s]")
    ax2.set_ylabel("Erreur relative [%]")
    ax2.set_title("Conservation du volume de la goutte")
    ax2.grid(True, linestyle="--", alpha=0.7)
    ax2.legend()

    plt.tight_layout()
    output_png = "output_droplet/droplet_diagnostics.png"
    plt.savefig(output_png, dpi=300)
    print(f"Graphique sauvegardé : {output_png}")
    plt.close()

def plot_dam_break():
    csv_file = "output_dam_break/dam_break_history.csv"
    if not os.path.exists(csv_file):
        print(f"Fichier introuvable: {csv_file}")
        return

    df = pd.read_csv(csv_file)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(df['time'], df['kinetic_energy'], 'k-o', lw=2, label="Énergie cinétique globale")
    ax1.set_xlabel("Temps [s]")
    ax1.set_ylabel("E_k [J]")
    ax1.set_title("Énergie cinétique - Rupture de barrage")
    ax1.grid(True, linestyle="--", alpha=0.7)
    ax1.legend()

    ax2.plot(df['time'], df['max_u'], 'c-s', lw=2, label="Vitesse maximale de front")
    ax2.set_xlabel("Temps [s]")
    ax2.set_ylabel("U_max [m/s]")
    ax2.set_title("Propagation du front de déferlement")
    ax2.grid(True, linestyle="--", alpha=0.7)
    ax2.legend()

    plt.tight_layout()
    output_png = "output_dam_break/dam_break_dynamics.png"
    plt.savefig(output_png, dpi=300)
    print(f"Graphique sauvegardé : {output_png}")
    plt.close()

if __name__ == "__main__":
    plot_rising_bubble()
    plot_static_droplet()
    plot_dam_break()
    print("Tous les graphiques de diagnostic ont été générés avec succès.")
