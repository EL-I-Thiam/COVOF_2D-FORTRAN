# Synthèse Comparative des Méthodes VOF — Cas 3.19 (Vortex Cisaillé)

| Méthode                      |   Erreur L1 | Conservation Masse   | Gain vs Baseline   | Statut    |
|:-----------------------------|------------:|:---------------------|:-------------------|:----------|
| VOF-TVD Superbee (Référence) |     0.0569  | 100.00%              | 1.00x              | RÉFÉRENCE |
| CIAM + PLIC (Géométrique)    |     0.00623 | 104.76%              | 9.24x              | CONFORME  |
| CIAM + SLIC (Géométrique)    |     0.00623 | 104.76%              | 9.24x              | CONFORME  |
| CIAM + ELVIRA (Géométrique)  |     0.00623 | 104.76%              | 9.24x              | CONFORME  |
| Weymouth & Yue + ELVIRA      |     0.01197 | 100.00%              | 4.81x              | CONFORME  |
| Weymouth & Yue + PLIC        |     0.01197 | 100.00%              | 4.81x              | CONFORME  |
| Weymouth & Yue + SLIC        |     0.01197 | 100.00%              | 4.81x              | CONFORME  |
| Weymouth & Yue + LVIRA       |     0.01197 | 100.00%              | 4.81x              | CONFORME  |
| CIAM + LVIRA                 |     0.00623 | 104.76%              | 9.24x              | CONFORME  |
| VOF-TVD HRIC (NVD)           |     0.04043 | 100.00%              | 1.42x              | CONFORME  |
| VOF-TVD MUSCL                |     0.06424 | 100.00%              | 0.90x              | CONFORME  |
| VOF-TVD Van Leer             |     0.06701 | 100.00%              | 0.86x              | CONFORME  |

*Comparaison entre schémas algébriques TVD, géométriques CIAM et multidimensionnels Weymouth & Yue.*
