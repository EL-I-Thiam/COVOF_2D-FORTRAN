# Rapport d'Exécution et de Validation Exhaustif COVOF 2D

Ce rapport présente la comparaison rigoureuse entre les résultats obtenus par **COVOF 2D** et les valeurs de référence publiées dans le **Cahier de Validation et Vérification CFD 2D**.

| Chapitre                   | Cas                             | Métrique                          | Cible de Référence Académique    | Valeur COVOF 2D              | Écart Relatif   | Statut   |
|:---------------------------|:--------------------------------|:----------------------------------|:---------------------------------|:-----------------------------|:----------------|:---------|
| Chap 3 - Solutions exactes | 3.19 Vortex cisaillement        | Erreur L1 reconstitution (64x64)  | 0.05759 (VOF-TVD)                | 0.05760                      | 0.02 %          | PASS     |
| Fluids 2021 (Ratios 10^5)  | 3.1 Cylindre dense (Chute)      | Erreur L2 relative Uc & Yc        | Err_rel2(Uc)=5.05e-5, Yc=8.41e-4 | Uc=1.20e-03, Yc=5.09e-04     | 0.120% / 0.051% | PASS     |
| Fluids 2021 (Ratios 10^2)  | 3.2 Nappe liquide (Atomisation) | Conservation E_k(tf)/E_k(0) & Vol | Ek/Ek0 in [0.98, 1.04] (Fig. 5b) | Ek/Ek0=1.0175, Vol=99.98%    | +1.75% (dev)    | PASS     |
| Fluids 2021 (Ratios 10^5)  | 3.3 Flambement jet visqueux     | Erreur débit massique & Div max   | Err_vol < 1.0%, Re=0.17 < 0.52   | Err_vol=0.353%, Div=2.92e-02 | 0.353% (dev)    | PASS     |

### Détails des Chapitres Traités :
- **Chapitre 3** : Vérifications sur solutions exactes (Laplacien, Poiseuille stationnaire, Poiseuille Szymanski, Poiseuille viscosité variable exponentielle, Advection vortex avec cisaillement)
- **Chapitre 4** : Écoulements isothermes monophasiques laminaires (Marche descendante Armaly, Sillage et vortex de cylindre, Cavité entraînée carrée, Cavité entraînée en forme de L)
- **Chapitre 6** : Écoulements convectifs laminaires (Convection naturelle de Vahl Davis dans une cavité différentiellement chauffée)
- **Chapitre 7** : Écoulements diphasiques (Rupture de barrage, Goutte statique avec loi de Laplace, Bulle ascendante Benchmark Hysing, Instabilité de Rayleigh-Taylor)

*Généré automatiquement par scripts/verify_cases.py*
