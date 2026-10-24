#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Générateur du Guide du Développeur Exhaustif (100+ pages LaTeX : 50+ pages C++17 & 50+ pages Fortran 2018)
pour le Solveur Diphasique et Thermoconvectif COVOF 2D.
Produit docs/DEVELOPER_GUIDE.tex, docs/DEVELOPER_GUIDE.md et compile docs/DEVELOPER_GUIDE.pdf.
"""

import os
import subprocess

ROOT_CPP = "/Users/elhadjiibrahimathiam/Documents/PROJET_ANTIGRAVITY/COVOF_2D-C++"
ROOT_F90 = "/Users/elhadjiibrahimathiam/Documents/PROJET_ANTIGRAVITY/COVOF_2D-FORTRAN"


def read_source(root, rel_path):
    full_path = os.path.join(root, rel_path)
    with open(full_path, "r", encoding="utf-8") as f:
        return f.read()


CPP_CHAPTERS = [
    {
        "num": 1,
        "title": "Architecture Globale C++17, Système de Build CMake et Types Fondamentaux (Types.hpp, main.cpp)",
        "files": [("CMakeLists.txt", "cmake"), ("include/covof/Types.hpp", "cpp"), ("src/main.cpp", "cpp")],
        "theory": r"""
\section{Philosophie Architecturale et Standard C++17}
L'implémentation \textbf{C++17} du solveur \textbf{COVOF 2D} (\texttt{COVOF\_2D-C++}) repose sur une architecture orientée objet modulaire à couplage faible et forte cohésion interne, encapsulée dans l'espace de noms \texttt{covof}. Les principes directeurs d'ingénierie logicielle retenus sont les suivants :
\begin{enumerate}
  \item \textbf{Zéro allocation dynamique dans la boucle temporelle} : tous les champs eulériens (\texttt{Field2D<double>}), les tampons d'échange de halos MPI et la structure de la matrice creuse CSR (\texttt{CSRMatrix}) sont pré-alloués lors de la phase d'initialisation (\texttt{Simulation::initialize}).
  \item \textbf{Contiguïté mémoire stricte (\emph{Row-Major Order})} : chaque champ 2D de dimension locale $N_x \times N_y$ entouré de $n_g$ couches de cellules fantômes (\emph{ghost cells}) est stocké dans un unique \texttt{std::vector<double>} unidimensionnel de taille $(N_x + 2n_g)(N_y + 2n_g)$. L'indice linéaire correspondant à la maille $(i, j)$ est donné par :
  \begin{equation}
    k_{\text{C++}}(i, j) = (j + n_g) \times (N_x + 2n_g) + (i + n_g)
  \end{equation}
  Ainsi, une boucle interne parcourant l'indice horizontal $i \in \{0, \dots, N_x-1\}$ à ordonnée $j$ fixée accède à des adresses mémoire strictement consécutives ($\Delta k = 1$), maximisant l'utilisation des lignes de cache L1/L2 de 64 octets et permettant la vectorisation automatique SIMD (AVX2 / AVX-512 / NEON Apple Silicon).
  \item \textbf{Typage fort des modèles physiques et numériques (\texttt{Types.hpp})} : l'ensemble des choix algorithmiques (type de condition aux limites \texttt{BoundaryType}, schéma d'advection VOF \texttt{AdvectionMethod}, méthode de reconstruction géométrique \texttt{ReconstructionMethod}, limiteur de flux \texttt{TVDLimiter}, solveur de Krylov \texttt{SolverType}, préconditionneur \texttt{PreconditionerType}, mode d'initialisation \texttt{InitialConditionType}) est représenté par des énumérations fortement typées (\texttt{enum class}) interdisant toute conversion implicite erronée.
\end{enumerate}

\section{Structure des Blocs de Configuration dans \texttt{Types.hpp}}
L'en-tête \texttt{include/covof/Types.hpp} définit également l'intégralité des structures de données de paramétrage (\emph{Plain Old Data structs}) qui reflètent bijectivement l'arborescence des fichiers JSON du répertoire \texttt{cases/*.json} :
\begin{itemize}
  \item \texttt{MeshConfig} : dimensions entières globales \texttt{nx}, \texttt{ny} et dimensions physiques \texttt{lx}, \texttt{ly} $[\text{m}]$.
  \item \texttt{PhysicsConfig} : masses volumiques \texttt{rho1}, \texttt{rho2} $[\text{kg/m}^3]$, viscosités dynamiques \texttt{mu1}, \texttt{mu2} $[\text{Pa}\cdot\text{s}]$, gravité $(\texttt{gx}, \texttt{gy})$ $[\text{m/s}^2]$, tension superficielle \texttt{sigma} $[\text{N/m}]$ et paramètres de viscosité variable exponentielle \texttt{variable\_viscosity}, \texttt{visc\_law}, \texttt{visc\_a}, \texttt{visc\_e}.
  \item \texttt{ThermalConfig} : drapeau \texttt{enable\_energy}, diffusivité thermique $\alpha$ (\texttt{alpha\_th}), coefficient d'expansion thermique de Boussinesq $\beta$ (\texttt{beta\_th}), température de référence $T_0$ (\texttt{t\_ref}) et températures imposées (\texttt{t\_hot}, \texttt{t\_cold}).
  \item \texttt{VofConfig} : méthode de reconstruction géométrique (\texttt{PLIC}, \texttt{ELVIRA}, \texttt{LVIRA}, \texttt{SLIC}, \texttt{NONE}), schéma d'advection (\texttt{WEYMOUTH\_YUE}, \texttt{CIAM}, \texttt{TVD}), limiteur TVD (\texttt{SUPERBEE}, \texttt{MC}, \texttt{VAN\_LEER}, \texttt{MINMOD}, \texttt{HRIC}) et nombre de Courant interfacial maximal \texttt{cfl\_vof}.
  \item \texttt{SolverConfig} : solveur de Krylov (\texttt{BICGSTAB2}, \texttt{BICGSTAB}, \texttt{PCG}, \texttt{GMRES}, \texttt{MINRES}), préconditionneur (\texttt{JACOBI}, \texttt{SSOR}, \texttt{ILU0}), tolérance relative \texttt{tol}, nombre maximal d'itérations \texttt{max\_iter}, paramètre de pénalisation du Lagrangien Augmenté $r_{al}$ (\texttt{r\_al}) et nombre de sous-itérations d'Uzawa (\texttt{uzawa\_iter}).
\end{itemize}

\section{Cycle de Vie du Programme Principal (\texttt{src/main.cpp})}
Le point d'entrée \texttt{main(int argc, char* argv[])} orchestre le cycle de vie MPI et l'exécution complète d'un cas test en quatre étapes :
\begin{enumerate}
  \item Initialisation de l'environnement parallèle via \texttt{MPI\_Init(\&argc, \&argv)} et interrogation du rang global \texttt{rank} et du nombre de processus \texttt{size}.
  \item Lecture et validation du fichier JSON passé en argument de ligne de commande via \texttt{Config::load\_from\_file(json\_path)}.
  \item Instanciation de l'objet maître \texttt{Simulation sim(cfg, MPI\_COMM\_WORLD)}, appel à \texttt{sim.initialize()} puis exécution de la boucle temporelle \texttt{sim.run()} chronométrée par \texttt{MPI\_Wtime()}.
  \item Libération ordonnée des ressources MPI via \texttt{MPI\_Finalize()}.
\end{enumerate}
"""
    },
    {
        "num": 2,
        "title": "Analyseur Syntaxique JSON et Gestion de la Configuration (Config.hpp, Config.cpp)",
        "files": [("include/covof/Config.hpp", "cpp"), ("src/Config.cpp", "cpp")],
        "theory": r"""
\section{Conception de l'Analyseur JSON Autonome Sans Dépendance Externe}
Afin de garantir une portabilité immédiate sur tout supercalculateur ou station de travail sans nécessiter le téléchargement de bibliothèques tierces (telles que \texttt{nlohmann/json} ou \texttt{Boost.PropertyTree}), \texttt{COVOF\_2D-C++} intègre dans \texttt{src/Config.cpp} un analyseur syntaxique JSON déterministe et léger.

\subsection{Grammaire Reconnue et Extraction par Sous-Blocs}
Le parseur opère en trois phases :
\begin{enumerate}
  \item \textbf{Prétraitement et nettoyage lexical} : lecture intégrale du flux fichier dans une chaîne \texttt{std::string}, gestion des espaces, retours chariot et commentaires éventuels.
  \item \textbf{Isolation des objets JSON de premier niveau (\texttt{extract\_section})} : pour une clé de section donnée (par exemple \texttt{"mesh"}, \texttt{"physics"}, \texttt{"vof"}, \texttt{"solver"}), l'automate recherche l'ouverture d'accolade \texttt{\{} correspondante et maintient un compteur de profondeur d'imbrication $\text{depth} = \sum (\texttt{'\{'} - \texttt{'\}'})$ tout en ignorant les accolades situées à l'intérieur de littéraux chaînes entre guillemets. Cela permet d'isoler exactement la sous-chaîne JSON du bloc, y compris lorsqu'il contient des tableaux d'objets imbriqués comme \texttt{"obstacles": [ ... ]}.
  \item \textbf{Désérialisation typée avec valeurs par défaut sécurisées} :
  \begin{itemize}
    \item \texttt{get\_string(section, key, default\_val)} : extrait les chaînes de caractères et normalise la casse pour la résolution des énumérations ;
    \item \texttt{get\_double(section, key, default\_val)} : parse les nombres flottants en notation décimale ou scientifique (\texttt{std::stod}) ;
    \item \texttt{get\_int(section, key, default\_val)} et \texttt{get\_bool(section, key, default\_val)} : extraient les entiers et les booléens (\texttt{true}/\texttt{false}).
  \end{itemize}
\end{enumerate}

\section{Compatibilité Ascendante et Alias de Clés JSON}
Pour assurer une interopérabilité totale entre les différentes générations de fichiers de cas tests et entre les dépôts C++17 et Fortran 2018, \texttt{Config::load\_from\_file} implémente une résolution d'alias intelligente :
\begin{itemize}
  \item Pour le solveur linéaire : reconnaît aussi bien \texttt{"tol"} que \texttt{"tolerance"}, et \texttt{"r\_al"} que \texttt{"augmented\_lagrangian\_r"}.
  \item Pour le module VOF : reconnaît la section dédiée \texttt{"vof"} (\texttt{reconstruction\_method}, \texttt{advection\_method}, \texttt{tvd\_limiter}, \texttt{cfl\_vof}) ainsi que les clés historiques éventuellement placées dans \texttt{"solver"} ou \texttt{"physics"}.
  \item Pour les obstacles immergés de Brinkman : parse dynamiquement le tableau JSON \texttt{"obstacles"} en distinguant les primitives géométriques rectangulaires (\texttt{type = "box"}, définies par $[x_{\min}, x_{\max}] \times [y_{\min}, y_{\max}]$) et circulaires (\texttt{type = "circle"}, définies par le centre $(x_c, y_c)$ et le rayon $R$).
\end{itemize}
"""
    },
    {
        "num": 3,
        "title": "Topologie du Maillage Cartésien MAC 2D et Décomposition Spatiale (Mesh.hpp, Mesh.cpp)",
        "files": [("include/covof/Mesh.hpp", "cpp"), ("src/Mesh.cpp", "cpp")],
        "theory": r"""
\section{Géométrie de la Grille Décalée MAC (\emph{Marker-And-Cell})}
La classe \texttt{covof::Mesh} (\texttt{include/covof/Mesh.hpp}, \texttt{src/Mesh.cpp}) encapsule la géométrie discrète du sous-domaine attribué au processus MPI courant. Soit un domaine global $[0, L_x] \times [0, L_y]$ discrétisé par $N_{x,\text{glob}} \times N_{y,\text{glob}}$ cellules uniformes de pas :
\begin{equation}
  \Delta x = \frac{L_x}{N_{x,\text{glob}}}, \qquad \Delta y = \frac{L_y}{N_{y,\text{glob}}}
\end{equation}
Lorsque le domaine est partitionné sur une grille cartésienne de $P_x \times P_y$ processus MPI, le processus de coordonnées cartésiennes $(p_x, p_y) \in \{0, \dots, P_x-1\} \times \{0, \dots, P_y-1\}$ reçoit un sous-domaine local de $N_x \times N_y$ cellules débutant aux indices globaux $(i_{\text{start}}, j_{\text{start}})$ :
\begin{equation}
  N_x = \left\lfloor \frac{N_{x,\text{glob}}}{P_x} \right\rfloor + \mathbb{I}(p_x < N_{x,\text{glob}} \bmod P_x)
\end{equation}
Ce partitionnement à reste distribué garantit un déséquilibre de charge maximal d'au plus une colonne ou une ligne de mailles entre deux processus quelconques.

\section{Systèmes de Coordonnées Décalées $(x_c, y_c)$, $(x_u, y_c)$ et $(x_c, y_v)$}
En indiçage local C++17 à base zéro ($i \in \{0, \dots, N_x-1\}$, $j \in \{0, \dots, N_y-1\}$), les trois familles de points de collocation de la grille MAC sont évaluées par les méthodes en ligne (\texttt{inline}) de \texttt{Mesh} :
\begin{enumerate}
  \item \textbf{Centre de la maille $(i, j)$} (support de $p, T, C, \rho, \mu, \nabla\cdot\mathbf{u}$) :
  \begin{equation}
    x_c(i) = x_{\min,\text{loc}} + \left(i + \frac{1}{2}\right)\Delta x, \qquad y_c(j) = y_{\min,\text{loc}} + \left(j + \frac{1}{2}\right)\Delta y
  \end{equation}
  \item \textbf{Face verticale $(i, j)$} (support de la vitesse horizontale $u_{i,j}$ située sur la face Ouest de la maille $(i,j)$, avec $i \in \{0, \dots, N_x\}$) :
  \begin{equation}
    x_u(i) = x_{\min,\text{loc}} + i\,\Delta x, \qquad y_u(j) = y_c(j)
  \end{equation}
  \item \textbf{Face horizontale $(i, j)$} (support de la vitesse verticale $v_{i,j}$ située sur la face Sud de la maille $(i,j)$, avec $j \in \{0, \dots, N_y\}$) :
  \begin{equation}
    x_v(i) = x_c(i), \qquad y_v(j) = y_{\min,\text{loc}} + j\,\Delta y
  \end{equation}
\end{enumerate}
"""
    },
    {
        "num": 4,
        "title": "Conteneurs de Champs Eulériens 2D et Opérations Vectorielles (Field.hpp)",
        "files": [("include/covof/Field.hpp", "cpp")],
        "theory": r"""
\section{Conception du Template \texttt{Field2D<T>}}
L'en-tête \texttt{include/covof/Field.hpp} fournit la classe générique \texttt{Field2D<T>} (instanciée principalement pour \texttt{T = double} via l'alias \texttt{ScalarField}) qui représente tout champ discret défini sur le maillage MAC (aussi bien un champ centré de dimension $N_x \times N_y$ qu'un champ de vitesse décalé $u$ de dimension $(N_x+1) \times N_y$ ou $v$ de dimension $N_x \times (N_y+1)$).

\subsection{Adressage Avec Cellules Fantômes (\emph{Ghost Cells})}
Chaque instance de \texttt{Field2D<T>} stocke :
\begin{itemize}
  \item \texttt{nx\_}, \texttt{ny\_} : les dimensions actives du champ à l'intérieur du sous-domaine ;
  \item \texttt{ng\_} : l'épaisseur de la couronne de mailles fantômes (typiquement $n_g = 1$ ou $n_g = 2$ pour les stencils VOF $3\times 3$ et TVD) ;
  \item \texttt{stride\_} $= N_x + 2 n_g$ : le pas d'adressage mémoire d'une ligne complète ;
  \item \texttt{data\_} : le vecteur contigu de taille $(N_x + 2n_g)(N_y + 2n_g)$ initialisé à zéro.
\end{itemize}
L'opérateur d'accès \texttt{operator()(int i, int j)} accepte des indices relatifs signés $i \in [-n_g, N_x + n_g - 1]$ et $j \in [-n_g, N_y + n_g - 1]$ :
\begin{equation}
  \texttt{field(i, j)} \equiv \texttt{data\_[(j + ng\_) * stride\_ + (i + ng\_)]}
\end{equation}
Ainsi, l'écriture des opérateurs aux différences finies dans \texttt{NavierStokesSolver.cpp} et \texttt{InterfaceTracker.cpp} conserve une lisibilité mathématique directe (\texttt{u(i+1, j) - 2*u(i, j) + u(i-1, j)}) sans aucune arithmétique d'offset dans le code client.
"""
    },
    {
        "num": 5,
        "title": "Parallélisme Distribué OpenMPI et Échanges de Halos Cartésiens (MpiComm.hpp, MpiComm.cpp)",
        "files": [("include/covof/MpiComm.hpp", "cpp"), ("src/MpiComm.cpp", "cpp")],
        "theory": r"""
\section{Topologie Cartésienne 2D et Découverte des Voisins}
La classe \texttt{covof::MpiComm} (\texttt{include/covof/MpiComm.hpp}, \texttt{src/MpiComm.cpp}) abstrait toutes les communications inter-processus OpenMPI. Lors de sa construction :
\begin{enumerate}
  \item La fonction \texttt{MPI\_Dims\_create(size, 2, dims)} décompose automatiquement le nombre total de processus $P$ en une grille bidimensionnelle $P_x \times P_y$ minimisant le rapport périmètre/surface des sous-domaines (et donc le volume de communication de halo par rapport au volume de calcul).
  \item \texttt{MPI\_Cart\_create} construit le communicateur topologique 2D en prenant en compte les drapeaux de périodicité $(\text{periodic}_x, \text{periodic}_y)$ issus de la configuration JSON.
  \item \texttt{MPI\_Cart\_shift} identifie les rangs des quatre voisins immédiats : \texttt{rank\_left} (Ouest, $x^-$), \texttt{rank\_right} (Est, $x^+$), \texttt{rank\_bottom} (Sud, $y^-$) et \texttt{rank\_top} (Nord, $y^+$). Si une frontière externe n'est pas périodique, MPI affecte automatiquement \texttt{MPI\_PROC\_NULL}, ce qui transforme les appels d'envoi/réception sur ce bord en opérations nulles sans branchement conditionnel coûteux.
\end{enumerate}

\section{Protocole d'Échange de Halos Sans Interblocage (\emph{Deadlock-Free Halo Exchange})}
La méthode \texttt{MpiComm::exchange\_halo(Field2D<double>\& field)} synchronise les couches de cellules fantômes $i \in \{-n_g, \dots, -1\}$, $i \in \{N_x, \dots, N_x+n_g-1\}$, $j \in \{-n_g, \dots, -1\}$ et $j \in \{N_y, \dots, N_y+n_g-1\}$ en deux passes séquentielles (passe horizontale Est$\leftrightarrow$Ouest suivie de la passe verticale Nord$\leftrightarrow$Sud) :
\begin{itemize}
  \item L'enchaînement séquentiel de la passe horizontale puis de la passe verticale propage automatiquement les \textbf{cellules fantômes de coin} $(i < 0, j < 0)$, $(i \ge N_x, j \ge N_y)$ sans nécessiter de messages diagonaux dédiés (crucial pour le filtre $3\times 3$ de Parker-Youngs et le terme croisé $-r_{al}\partial^2 / \partial x \partial y$).
  \item Les colonnes verticales n'étant pas contiguës en mémoire \emph{row-major}, elles sont empaquetées (\emph{packing}) dans des tampons linéaires \texttt{send\_left}, \texttt{send\_right}, transmises via \texttt{MPI\_Sendrecv} ou \texttt{MPI\_Isend}/\texttt{MPI\_Irecv}, puis dépaquetées (\emph{unpacking}) dans les cellules fantômes.
\end{itemize}
"""
    },
    {
        "num": 6,
        "title": "Algèbre Linéaire Creuse au Format CSR (CSRMatrix.hpp, CSRMatrix.cpp)",
        "files": [("include/covof/CSRMatrix.hpp", "cpp"), ("src/CSRMatrix.cpp", "cpp")],
        "theory": r"""
\section{Représentation Compressed Sparse Row (CSR)}
Le système linéaire monolithique du Lagrangien Augmenté comporte $N_{\text{dof}} = N_u + N_v$ degrés de liberté par sous-domaine, mais chaque ligne de la matrice ne possède qu'au plus 9 coefficients non nuls (5 termes directs du Laplacien de la composante courante + 4 termes croisés du gradient de divergence couplant $u$ et $v$). Stocker une matrice dense nécessiterait $\mathcal{O}(N_{\text{dof}}^2)$ octets (plusieurs gigaoctets dès $128\times 128$), tandis que le format creux \textbf{CSR (\emph{Compressed Sparse Row})} implémenté dans \texttt{covof::CSRMatrix} réduit l'empreinte mémoire à $\mathcal{O}(9 N_{\text{dof}})$ (quelques mégaoctets).

Une matrice creuse $\mathbf{A} \in \mathbb{R}^{n \times n}$ comportant $\text{nnz}$ éléments non nuls est représentée par trois tableaux contigus :
\begin{enumerate}
  \item \texttt{values} ($\in \mathbb{R}^{\text{nnz}}$) : valeurs flottantes des éléments non nuls parcourus ligne par ligne ;
  \item \texttt{col\_ind} ($\in \mathbb{N}^{\text{nnz}}$) : indices de colonnes $j$ correspondants ;
  \item \texttt{row\_ptr} ($\in \mathbb{N}^{n+1}$) : indice de début de chaque ligne $i$ dans \texttt{values} et \texttt{col\_ind}, avec $\texttt{row\_ptr}[0] = 0$ et $\texttt{row\_ptr}[n] = \text{nnz}$.
\end{enumerate}

\section{Produit Matrice-Vecteur Creux (\texttt{multiply} / SpMV)}
Pour calculer $\mathbf{y} = \mathbf{A}\mathbf{x}$, la méthode \texttt{CSRMatrix::multiply(const std::vector<double>\& x, std::vector<double>\& y)} évalue pour chaque ligne $i \in \{0, \dots, n-1\}$ :
\begin{equation}
  y_i = \sum_{k = \texttt{row\_ptr}[i]}^{\texttt{row\_ptr}[i+1]-1} \texttt{values}[k] \, x_{\texttt{col\_ind}[k]}
\end{equation}
Cette opération s'exécute en exactement $2\,\text{nnz}$ opérations flottantes (FLOPs) avec un parcours purement séquentiel des tableaux \texttt{values} et \texttt{col\_ind}.
"""
    },
    {
        "num": 7,
        "title": "Préconditionneurs Algébriques : Jacobi, SSOR et ILU(0) (Preconditioner.hpp, Preconditioner.cpp)",
        "files": [("include/covof/Preconditioner.hpp", "cpp"), ("src/Preconditioner.cpp", "cpp")],
        "theory": r"""
\section{Rôle du Préconditionnement dans les Méthodes de Krylov}
Pour accélérer la convergence des solveurs itératifs sur le système $\mathbf{A}\mathbf{x} = \mathbf{b}$ dont le nombre de conditionnement $\kappa(\mathbf{A})$ croît avec le paramètre du Lagrangien Augmenté $r_{al}$ et le rapport de viscosité $\mu_2 / \mu_1$, \texttt{COVOF\_2D-C++} définit l'interface abstraite \texttt{covof::Preconditioner} qui résout à chaque itération de Krylov le système préconditionné $\mathbf{M}\mathbf{z} = \mathbf{r}$ (ou applique $\mathbf{z} = \mathbf{M}^{-1}\mathbf{r}$).

\subsection{Préconditionneur Diagonal de Jacobi (\texttt{JacobiPreconditioner})}
On choisit $\mathbf{M} = \text{diag}(A_{0,0}, A_{1,1}, \dots, A_{n-1,n-1})$. Lors de la phase \texttt{setup(A)}, l'inverse de la diagonale est précalculé avec garde-fou contre les pivots nuls :
\begin{equation}
  d_i^{-1} = \begin{cases} 1 / A_{i,i} & \text{si } |A_{i,i}| > 10^{-15} \\ 1.0 & \text{sinon} \end{cases}, \qquad z_i = d_i^{-1} r_i
\end{equation}
Dans un système à forts sauts de densité ($\rho_2/\rho_1 = 10^5$) et de viscosité ($\mu_2/\mu_1 = 10^{10}$), ce préconditionneur diagonal remet automatiquement à l'échelle les équations situées dans le gaz et dans le liquide dense, ramenant les éléments diagonaux de $\mathbf{M}^{-1}\mathbf{A}$ à l'unité.

\subsection{Préconditionneur SSOR et Factorisation Incomplète \texttt{ILU(0)}}
\begin{itemize}
  \item \textbf{SSOR (\emph{Symmetric Successive Over-Relaxation})} : décompose $\mathbf{A} = \mathbf{D} + \mathbf{L} + \mathbf{U}$ et approche $\mathbf{M}_{\text{SSOR}} = \frac{1}{\omega(2-\omega)}(\mathbf{D} + \omega\mathbf{L})\mathbf{D}^{-1}(\mathbf{D} + \omega\mathbf{U})$ inversé par une descente triangulaire inférieure suivie d'une remontée triangulaire supérieure.
  \item \textbf{ILU(0) (\emph{Incomplete LU with zero fill-in})} : calcule une factorisation $\mathbf{A} \approx \tilde{\mathbf{L}}\tilde{\mathbf{U}}$ conservant exactement le même motif de parcimonie (\emph{sparsity pattern}) que la matrice $\mathbf{A}$ originale.
\end{itemize}
"""
    },
    {
        "num": 8,
        "title": "Solveurs de Sous-Espaces de Krylov : BiCGSTAB(2), PCG, GMRES (LinearSolver.hpp, LinearSolver.cpp)",
        "files": [("include/covof/LinearSolver.hpp", "cpp"), ("src/LinearSolver.cpp", "cpp")],
        "theory": r"""
\section{Architecture des Solveurs de Krylov Parallèles}
L'en-tête \texttt{include/covof/LinearSolver.hpp} et son implémentation \texttt{src/LinearSolver.cpp} (796 lignes) regroupent cinq algorithmes de sous-espaces de Krylov capables de fonctionner en séquentiel comme en parallèle MPI (via les produits scalaires globaux \texttt{MpiComm::dot} et normes globales \texttt{MpiComm::norm2}).

\section{Algorithme Mathématique Complet de BiCGSTAB(2) (Sleijpen \& Fokkema 1993)}
Alors que l'algorithme BiCGSTAB standard stabilise la suite de Petrov-Galerkin du gradient biconjugué par un polynôme local de degré 1 ($P_1(\mathbf{A}) = \mathbf{I} - \omega_k \mathbf{A}$) qui ne possède que des racines réelles et peut stagner lorsque l'opérateur d'advection ou le Lagrangien Augmenté possède des valeurs propres à forte partie imaginaire, \textbf{BiCGSTAB($\ell$)} avec $\ell = 2$ combine deux pas de BiCG avec une minimisation résiduelle minimale (MR / GMRES(2)) engendrant un polynôme quadratique capable de capturer des paires de valeurs propres complexes conjuguées :
\begin{enumerate}
  \item \textbf{Initialisation} : $\mathbf{r}_0 = \mathbf{b} - \mathbf{A}\mathbf{x}_0$, choix de l'ombre $\tilde{\mathbf{r}}_0 = \mathbf{r}_0$, $\rho_0 = 1$, $\alpha = 0$, $\omega_2 = 1$.
  \item \textbf{Phase BiCG sur $j = 0, 1$} : calcul successif des directions conjuguées préconditionnées et des produits matrice-vecteur $\hat{\mathbf{u}}_0, \hat{\mathbf{u}}_1$ et $\hat{\mathbf{r}}_0, \hat{\mathbf{r}}_1, \hat{\mathbf{r}}_2$.
  \item \textbf{Phase de Minimisation Polynomiale d'Ordre 2 (Orthogonalisation de Gram-Schmidt sur $\{\hat{\mathbf{r}}_1, \hat{\mathbf{r}}_2\}$)} :
  \begin{equation}
    \tau_{1,2} = \frac{\langle \hat{\mathbf{r}}_1, \hat{\mathbf{r}}_2 \rangle}{\langle \hat{\mathbf{r}}_1, \hat{\mathbf{r}}_1 \rangle}, \qquad \omega_1, \omega_2 \text{ minimisant } \|\hat{\mathbf{r}}_0 - \omega_1 \hat{\mathbf{r}}_1 - \omega_2 \hat{\mathbf{r}}_2\|_2
  \end{equation}
  \item \textbf{Mise à jour de la solution $\mathbf{x}$ et du résidu $\mathbf{r}$}, suivie du test d'arrêt $\|\mathbf{r}\|_2 / \|\mathbf{r}_0\|_2 < \varepsilon_{\text{tol}}$.
\end{enumerate}
Cet algorithme constitue la clé de voûte de la robustesse de COVOF 2D sur les cas diphasiques extrêmes (\emph{Fluids} 2021 Cas 3.1, 3.2 et 3.3).
"""
    },
    {
        "num": 9,
        "title": "Solveur Monolithique de Navier-Stokes par Lagrangien Augmenté (NavierStokesSolver.hpp, NavierStokesSolver.cpp)",
        "files": [("include/covof/NavierStokesSolver.hpp", "cpp"), ("src/NavierStokesSolver.cpp", "cpp")],
        "theory": r"""
\section{Formulation Monolithique et Numérotation Globale des Inconnues $(u, v)$}
Le fichier \texttt{src/NavierStokesSolver.cpp} (1132 lignes) est le cœur hydrodynamique et thermique du code C++17. Sur un sous-domaine de $N_x \times N_y$ cellules MAC :
\begin{itemize}
  \item Les vitesses horizontales $u_{i,j}$ ($i \in \{0, \dots, N_x\}$, $j \in \{0, \dots, N_y-1\}$) occupent les indices $I_u(i, j) = j (N_x + 1) + i \in \{0, \dots, N_u - 1\}$ avec $N_u = (N_x + 1) N_y$.
  \item Les vitesses verticales $v_{i,j}$ ($i \in \{0, \dots, N_x-1\}$, $j \in \{0, \dots, N_y\}$) occupent les indices décalés $I_v(i, j) = N_u + j N_x + i \in \{N_u, \dots, N_u + N_v - 1\}$ avec $N_v = N_x (N_y + 1)$.
\end{itemize}

\section{Discrétisation Exacte de l'Opérateur du Lagrangien Augmenté $-r_{al} \nabla_h (\nabla_h \cdot \mathbf{u})$}
Pour une face verticale intérieure $(i, j)$ associée à l'inconnue $u_{i,j}$ (située entre la maille gauche $(i-1, j)$ et la maille droite $(i, j)$), la composante horizontale du gradient de la divergence s'écrit :
\begin{equation}
  -r_{al} \left.\frac{\partial}{\partial x}(\nabla_h \cdot \mathbf{u})\right|_{i,j} = -\frac{r_{al}}{\Delta x} \left( D_{i,j} - D_{i-1,j} \right)
\end{equation}
En remplaçant les divergences centrées $D_{i,j} = \frac{u_{i+1,j} - u_{i,j}}{\Delta x} + \frac{v_{i,j+1} - v_{i,j}}{\Delta y}$ et $D_{i-1,j} = \frac{u_{i,j} - u_{i-1,j}}{\Delta x} + \frac{v_{i-1,j+1} - v_{i-1,j}}{\Delta y}$, on obtient les contributions implicites exactes assemblées dans la ligne $I_u(i,j)$ de la matrice $\mathbf{A}_{\text{aug}}$ :
\begin{equation}
  +\frac{2 r_{al}}{\Delta x^2} u_{i,j} - \frac{r_{al}}{\Delta x^2} u_{i+1,j} - \frac{r_{al}}{\Delta x^2} u_{i-1,j} - \frac{r_{al}}{\Delta x \Delta y} \left( v_{i,j+1} - v_{i,j} - v_{i-1,j+1} + v_{i-1,j} \right)
\end{equation}
De manière duale, sur une face horizontale intérieure $(i, j)$ associée à l'inconnue $v_{i,j}$ (située entre la maille basse $(i, j-1)$ et la maille haute $(i, j)$), la ligne $I_v(i,j)$ reçoit :
\begin{equation}
  +\frac{2 r_{al}}{\Delta y^2} v_{i,j} - \frac{r_{al}}{\Delta y^2} v_{i,j+1} - \frac{r_{al}}{\Delta y^2} v_{i,j-1} - \frac{r_{al}}{\Delta x \Delta y} \left( u_{i+1,j} - u_{i,j} - u_{i+1,j-1} + u_{i,j-1} \right)
\end{equation}

\section{Tenseur Visqueux à Viscosité Variable et Moyenne Harmonique Nodale}
Pour le terme visqueux complet $\nabla \cdot [\mu(\nabla \mathbf{u} + \nabla \mathbf{u}^T)]$, les dérivées normales sur l'équation en $u_{i,j}$ utilisent les viscosités aux centres des mailles $\mu_{i,j}$ et $\mu_{i-1,j}$, tandis que les dérivées transversales utilisent la viscosité aux sommets $\mu_{\text{node}}(i, j)$ et $\mu_{\text{node}}(i, j+1)$ évaluée par moyenne harmonique des 4 cellules entourant le sommet :
\begin{equation}
  \mu_{\text{node}}(i, j) = \frac{4}{\dfrac{1}{\mu_{i-1,j-1}} + \dfrac{1}{\mu_{i,j-1}} + \dfrac{1}{\mu_{i-1,j}} + \dfrac{1}{\mu_{i,j}}}
\end{equation}
"""
    },
    {
        "num": 10,
        "title": "Reconstruction Géométrique d'Interface VOF : PLIC, ELVIRA, LVIRA, SLIC (InterfaceReconstruction.hpp/.cpp)",
        "files": [("include/covof/InterfaceReconstruction.hpp", "cpp"), ("src/InterfaceReconstruction.cpp", "cpp")],
        "theory": r"""
\section{Représentation Linéaire par Morceaux PLIC (\emph{Piecewise Linear Interface Calculation})}
Dans chaque maille diphasique $(i, j)$ telle que $\varepsilon_{\text{vof}} < C_{i,j} < 1 - \varepsilon_{\text{vof}}$ ($\varepsilon_{\text{vof}} = 10^{-12}$), l'interface est représentée dans le repère local normalisé $(\xi, \eta) \in [0, 1]^2$ de la maille par la droite :
\begin{equation}
  n_x \xi + n_y \eta = \alpha
\end{equation}
où $\mathbf{n} = (n_x, n_y)^T$ est la normale unitaire pointant vers l'extérieur du fluide 1 ($\mathbf{n} = -\nabla C / \|\nabla C\|$) et $\alpha \in \mathbb{R}$ est la constante de positionnement géométrique unique telle que l'aire du polygone tronqué $\{(\xi, \eta) \in [0, 1]^2 \mid n_x \xi + n_y \eta \le \alpha\}$ soit exactement égale à la fraction volumique $C_{i,j}$.

\section{Relations Analytiques Directe et Inverse de Scardovelli \& Zaleski (2000)}
En se ramenant par symétrie au premier quadrant ($m_1 = \min(|n_x|, |n_y|) \ge 0$, $m_2 = \max(|n_x|, |n_y|) > 0$, avec $m_1 + m_2 = 1$ en norme $L_1$), l'inversion analytique $C \mapsto \alpha$ implémentée dans \texttt{compute\_line\_constant} distingue deux régimes géométriques sans aucune itération de Newton :
\begin{enumerate}
  \item \textbf{Régime triangulaire} ($0 \le C \le \frac{m_1}{2 m_2}$) :
  \begin{equation}
    \alpha = \sqrt{2 m_1 m_2 C}
  \end{equation}
  \item \textbf{Régime trapézoïdal} ($\frac{m_1}{2 m_2} < C \le \frac{1}{2}$) :
  \begin{equation}
    \alpha = m_2 C + \frac{m_1}{2}
  \end{equation}
\end{enumerate}
Pour $C > 1/2$, la symétrie centrale du carré unité donne immédiatement $\alpha(C) = 1 - \alpha(1 - C)$.

\section{Estimation des Normales : Parker-Youngs, ELVIRA et LVIRA}
\begin{itemize}
  \item \textbf{Parker-Youngs} : différence finie pondérée sur le stencil $3\times 3$ avec poids $1, 2, 1$ privilégiant la direction axiale.
  \item \textbf{ELVIRA (Pilliod \& Puckett 2004)} : évalue 6 pentes candidates (différences arrière, centrée et avant des sommes de colonnes horizontales $H_i = \sum_{l=-1}^1 C_{i, j+l}$ et des sommes de lignes verticales $V_j = \sum_{k=-1}^1 C_{i+k, j}$), reconstruit la droite associée sur le bloc $3\times 3$ et retient celle qui minimise l'erreur quadratique $\sum_{k,l=-1}^1 (C^{\text{rec}}_{i+k,j+l}(\mathbf{n}) - C_{i+k,j+l})^2$. Cet algorithme reconstruit exactement toute interface plane quelle que soit son orientation.
  \item \textbf{LVIRA (Puckett 1991)} : affine l'angle $\theta = \text{atan2}(n_y, n_x)$ de la normale ELVIRA par recherche locale d'or / sécante minimisant la fonctionnelle d'erreur $E_2(\theta)$ sur le stencil $3\times 3$.
\end{itemize}
"""
    },
    {
        "num": 11,
        "title": "Transport et Advection Conservatrice VOF : Weymouth-Yue, CIAM et TVD (InterfaceTracker.hpp/.cpp)",
        "files": [("include/covof/InterfaceTracker.hpp", "cpp"), ("src/InterfaceTracker.cpp", "cpp")],
        "theory": r"""
\section{Séparation Directionnelle d'Ordre 2 de Strang (\emph{Strang Operator Splitting})}
La classe \texttt{covof::InterfaceTracker} (\texttt{include/covof/InterfaceTracker.hpp}, \texttt{src/InterfaceTracker.cpp}, 833 lignes) résout l'équation de transport de la fraction volumique :
\begin{equation}
  \frac{\partial C}{\partial t} + \frac{\partial (u C)}{\partial x} + \frac{\partial (v C)}{\partial y} = C \left( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} \right)
\end{equation}
en alternant l'ordre des balayages unidimensionnels entre les pas de temps pairs ($X$ puis $Y$) et impairs ($Y$ puis $X$) afin d'éliminer tout biais directionnel et de garantir l'ordre 2 temporel.

\section{Schéma de Weymouth \& Yue (2010) à Conservation de Masse Exacte}
Lors d'un balayage 1D individuel (par exemple selon $x$), le champ de vitesse 1D n'est pas à divergence nulle ($\partial u / \partial x \ne 0$, seule la somme $\partial u / \partial x + \partial v / \partial y = 0$ est nulle). Weymouth \& Yue (2010) démontrent que si le coefficient $C_{\text{dil}}$ multipliant les termes de dilatation 1D est fixé à la valeur binaire évaluée au début du pas de temps :
\begin{equation}
  C_{\text{dil}}(i, j) = \begin{cases} 1 & \text{si } C_{i,j}^n \ge 0.5 \\ 0 & \text{si } C_{i,j}^n < 0.5 \end{cases}
\end{equation}
et maintenu identique durant la passe $X$ et la passe $Y$ :
\begin{align}
  C_{i,j}^* &= C_{i,j}^n - \frac{\Delta t}{\Delta x} \left( F_{i+1/2, j}^x - F_{i-1/2, j}^x \right) + C_{\text{dil}}(i,j) \frac{\Delta t}{\Delta x} \left( u_{i+1/2, j} - u_{i-1/2, j} \right) \\
  C_{i,j}^{n+1} &= C_{i,j}^* - \frac{\Delta t}{\Delta y} \left( F_{i, j+1/2}^y - F_{i, j-1/2}^y \right) + C_{\text{dil}}(i,j) \frac{\Delta t}{\Delta y} \left( v_{i, j+1/2} - v_{i, j-1/2} \right)
\end{align}
alors, d'une part $C_{i,j}^{n+1} \in [0, 1]$ est garanti sans écrêtage (*no-clipping*) sous $\text{CFL} < 0.5$, et d'autre part la somme des deux termes de dilatation vaut exactement $C_{\text{dil}} \Delta t (\nabla_h \cdot \mathbf{u})_{i,j} = 0$, assurant une **conservation de la masse exacte au bit près**.
"""
    },
    {
        "num": 12,
        "title": "Exportation Graphique VTK XML Parallèle (.vtr / .pvtr) (VTKWriter.hpp, VTKWriter.cpp)",
        "files": [("include/covof/VTKWriter.hpp", "cpp"), ("src/VTKWriter.cpp", "cpp")],
        "theory": r"""
\section{Architecture des Fichiers VTK XML \texttt{RectilinearGrid} et \texttt{PRectilinearGrid}}
La classe \texttt{covof::VTKWriter} (\texttt{include/covof/VTKWriter.hpp}, \texttt{src/VTKWriter.cpp}) produit des fichiers directement exploitables par \textbf{ParaView}, \textbf{VisIt} et les scripts Python de post-traitement sans bibliothèque externe :
\begin{enumerate}
  \item Chaque processus MPI de rang $r$ écrit son sous-domaine local dans \texttt{<case>\_p<r>\_<step>.vtr} (format \texttt{RectilinearGrid}) contenant les coordonnées des nœuds $(x_u, y_v, 0)$ et les tableaux \texttt{CellData} :
  \begin{itemize}
    \item \texttt{Pressure} ($p_{i,j}$) ;
    \item \texttt{Temperature} ($T_{i,j}$) ;
    \item \texttt{PhaseIndicator} ($C_{i,j}$) ;
    \item \texttt{Divergence} ($(\nabla_h \cdot \mathbf{u})_{i,j}$) ;
    \item \texttt{Vorticity} ($\omega_z(i,j) = \frac{v_{c}(i+1,j)-v_{c}(i-1,j)}{2\Delta x} - \frac{u_{c}(i,j+1)-u_{c}(i,j-1)}{2\Delta y}$) ;
    \item \texttt{Viscosity} ($\mu_{i,j}$) ;
    \item \texttt{Velocity} (vecteur 3D $(u_c(i,j), v_c(i,j), 0)$ interpolé aux centres des mailles : $u_c(i,j) = \frac{1}{2}(u_{i,j} + u_{i+1,j})$, $v_c(i,j) = \frac{1}{2}(v_{i,j} + v_{i,j+1})$).
  \end{itemize}
  \item Le processus maître ($r = 0$) écrit le fichier chapeau \texttt{<case>\_<step>.pvtr} (\texttt{PRectilinearGrid}) référençant les extensions globales \texttt{Piece Extent} de tous les sous-domaines MPI.
\end{enumerate}
"""
    },
    {
        "num": 13,
        "title": "Orchestration de la Simulation, Pas de Temps Adaptatif et Diagnostics V&V (Simulation.hpp, Simulation.cpp)",
        "files": [("include/covof/Simulation.hpp", "cpp"), ("src/Simulation.cpp", "cpp")],
        "theory": r"""
\section{Orchestration Globale dans \texttt{covof::Simulation}}
La classe \texttt{covof::Simulation} (\texttt{include/covof/Simulation.hpp}, \texttt{src/Simulation.cpp}, 916 lignes) coordonne l'ensemble des composants du solveur :
\begin{enumerate}
  \item \textbf{Initialisation géométrique et physique (\texttt{setup\_initial\_conditions})} : génère les profils analytiques, disques, colonnes de barrage, nappes liquides cisaillement Kelvin-Helmholtz, perturbations cosinusoïdales de Rayleigh-Taylor et buses d'injection de jets visqueux avec sous-échantillonnage $8 \times 8$ par cellule pour initialiser une fraction volumique $C_{i,j}(0)$ précise au second ordre.
  \item \textbf{Contrôle du pas de temps adaptatif (\texttt{compute\_adaptive\_dt})} : évalue à chaque pas de temps la contrainte CFL multi-physique (convective, capillaire de Brackbill et VOF) :
  \begin{equation}
    \Delta t \le \min\left( \text{CFL}_{\max} \frac{\min(\Delta x, \Delta y)}{\max(\|\mathbf{u}\|_\infty, 10^{-6})}, \quad \sqrt{\frac{(\rho_1 + \rho_2)\min(\Delta x, \Delta y)^3}{4\pi \sigma}} \right)
  \end{equation}
  \item \textbf{Calcul des diagnostics intégraux et normes V\&V} : produit automatiquement les fichiers CSV de suivi temporel (\texttt{history.csv}, \texttt{dense\_cylinder\_history.csv}, \texttt{liquid\_sheet\_history.csv}, \texttt{viscous\_jet\_history.csv}, \texttt{vortex\_error.csv}, \texttt{laplacian\_summary.csv}).
\end{enumerate}
"""
    }
]


F90_CHAPTERS = [
    {
        "num": 14,
        "title": "Architecture Modulaire Fortran 2018, Précision IEEE-754 et Programme Principal (CMakeLists.txt, module_precision.f90, module_constants.f90, covof2d_main.f90)",
        "files": [("CMakeLists.txt", "cmake"), ("src/module_precision.f90", "fortran"), ("src/module_constants.f90", "fortran"), ("src/covof2d_main.f90", "fortran")],
        "theory": r"""
\section{Conception Modulaire selon la Norme Fortran 2018 (ISO/IEC 1539-1:2018)}
L'implémentation \textbf{Fortran 2018} (\texttt{COVOF\_2D-FORTRAN}) est le miroir algorithmique et mathématique exact de \texttt{COVOF\_2D-C++}. Elle est structurée en 13 modules modernes (\texttt{implicit none} systématique, attributs \texttt{intent(in)}, \texttt{intent(out)}, \texttt{intent(inout)}, procédures attachées aux types dérivés et tableaux allocatables à bornes personnalisées) et un programme principal \texttt{covof2d\_main.f90}.

\subsection{Gestion de la Précision Flottante (\texttt{module\_precision.f90})}
Le module \texttt{module\_precision} définit le paramètre de genre (\emph{kind parameter}) double précision IEEE-754 portable sur 64 bits :
\begin{equation}
  \texttt{integer, parameter :: dp = selected\_real\_kind(15, 307)}
\end{equation}
garantissant au moins 15 chiffres décimaux significatifs de mantisse et une plage d'exposants décimaux de $\pm 307$, strictement identique au type \texttt{double} de C++17.

\subsection{Table de Correspondance Bijective C++17 $\leftrightarrow$ Fortran 2018}
\begin{center}
\begin{tabular}{ll}
\hline
\textbf{Composant C++17 (\texttt{COVOF\_2D-C++})} & \textbf{Module Fortran 2018 (\texttt{COVOF\_2D-FORTRAN})} \\
\hline
\texttt{include/covof/Types.hpp} & \texttt{src/module\_precision.f90} + \texttt{src/module\_constants.f90} \\
\texttt{src/Config.cpp} (\texttt{covof::Config}) & \texttt{src/module\_json.f90} + \texttt{src/module\_config.f90} \\
\texttt{src/Mesh.cpp} (\texttt{covof::Mesh}) & \texttt{src/module\_mesh.f90} (\texttt{type(mesh\_2d\_t)}) \\
\texttt{include/covof/Field.hpp} (\texttt{Field2D<T>}) & \texttt{src/module\_fields.f90} (\texttt{type(flow\_fields\_t)}) \\
\texttt{src/MpiComm.cpp} (\texttt{covof::MpiComm}) & \texttt{src/module\_mpi\_comm.f90} (\texttt{type(mpi\_context\_t)}) \\
\texttt{src/CSRMatrix.cpp} + \texttt{Preconditioner.cpp} + \texttt{LinearSolver.cpp} & \texttt{src/module\_linear\_solvers.f90} \\
\texttt{src/NavierStokesSolver.cpp} & \texttt{src/module\_navier\_stokes.f90} \\
\texttt{src/InterfaceReconstruction.cpp} & \texttt{src/module\_interface\_reconstruction.f90} \\
\texttt{src/InterfaceTracker.cpp} & \texttt{src/module\_interface\_vof.f90} \\
\texttt{src/VTKWriter.cpp} & \texttt{src/module\_vtk\_writer.f90} \\
\texttt{src/Simulation.cpp} + \texttt{src/main.cpp} & \texttt{src/module\_simulation.f90} + \texttt{src/covof2d\_main.f90} \\
\hline
\end{tabular}
\end{center}
"""
    },
    {
        "num": 15,
        "title": "Analyseur Syntaxique JSON Natif en Fortran 2018 et Module de Configuration (module_json.f90, module_config.f90)",
        "files": [("src/module_json.f90", "fortran"), ("src/module_config.f90", "fortran")],
        "theory": r"""
\section{Architecture du Parseur JSON Pur Fortran 2018 (\texttt{module\_json.f90})}
Contrairement à Python ou C++, l'écosystème Fortran standard ne dispose pas nativement d'un parseur JSON dans la bibliothèque système. Le module \texttt{module\_json.f90} (449 lignes) implémente intégralement en Fortran 2018 standard un analyseur lexical et syntaxique capable de lire exactement les mêmes fichiers \texttt{cases/*.json} que l'exécutable C++17 :
\begin{enumerate}
  \item \texttt{json\_read\_file(filename, content, ierr)} : charge le fichier JSON complet dans une chaîne allocatable \texttt{character(len=:), allocatable :: content}.
  \item \texttt{json\_extract\_block(json\_str, block\_name, block\_content, found)} : localise un objet JSON nommé (tel que \texttt{"mesh"} ou \texttt{"vof"}) en suivant l'équilibrage des accolades \texttt{\{} et \texttt{\}} hors chaînes de caractères.
  \item \texttt{json\_get\_real}, \texttt{json\_get\_int}, \texttt{json\_get\_logical}, \texttt{json\_get\_string} : extraient les valeurs typées avec repli automatique sur la valeur par défaut si la clé est absente.
  \item \texttt{json\_parse\_obstacles} : parcourt le tableau JSON \texttt{"obstacles": [ \{...\}, \{...\} ]} pour instancier les masques de pénalisation de Brinkman (\texttt{box} et \texttt{circle}).
\end{enumerate}

\section{Structure de Configuration \texttt{simulation\_config\_t} (\texttt{module\_config.f90})}
Le module \texttt{module\_config.f90} (493 lignes) traduit les chaînes JSON en constantes entières symboliques (\texttt{BC\_WALL}, \texttt{BC\_PERIODIC}, \texttt{BC\_JET\_INLET}, \texttt{SOLVER\_BICGSTAB2}, \texttt{RECON\_PLIC}, \texttt{ADVECT\_WEYMOUTH\_YUE}, etc.) et vérifie la cohérence des paramètres physiques et numériques avant l'allocation du maillage.
"""
    },
    {
        "num": 16,
        "title": "Maillage MAC 2D et Tableaux Column-Major à Bornes Fantômes en Fortran 2018 (module_mesh.f90, module_fields.f90)",
        "files": [("src/module_mesh.f90", "fortran"), ("src/module_fields.f90", "fortran")],
        "theory": r"""
\section{Indiçage Naturel à Base 1 et Bornes Fantômes \texttt{(1-ng:nx+ng)} en Fortran 2018}
L'un des atouts majeurs de Fortran 2018 pour la mécanique des fluides numérique sur grille structurée est la possibilité d'allouer des tableaux multidimensionnels avec des bornes inférieures négatives ou nulles correspondant directement aux cellules fantômes :
\begin{itemize}
  \item \textbf{Champs scalaires centrés} ($p, T, C, \rho, \mu, \nabla\cdot\mathbf{u}$) :
  \begin{equation}
    \texttt{allocate(fields\%p(1-ng : nx+ng, 1-ng : ny+ng))}
  \end{equation}
  où le domaine intérieur physique occupe exactement $i \in \{1, \dots, N_x\}$ et $j \in \{1, \dots, N_y\}$.
  \item \textbf{Vitesse horizontale décalée $u$} ($N_x+1$ faces verticales intérieures/frontières d'indices $i \in \{1, \dots, N_x+1\}$) :
  \begin{equation}
    \texttt{allocate(fields\%u(1-ng : nx+1+ng, 1-ng : ny+ng))}
  \end{equation}
  \item \textbf{Vitesse verticale décalée $v$} ($N_y+1$ faces horizontales intérieures/frontières d'indices $j \in \{1, \dots, N_y+1\}$) :
  \begin{equation}
    \texttt{allocate(fields\%v(1-ng : nx+ng, 1-ng : ny+1+ng))}
  \end{equation}
\end{itemize}

\section{Ordre Mémoire \emph{Column-Major} et Ordre Optimal des Boucles Imbriquées}
En Fortran, les tableaux 2D \texttt{A(i, j)} sont stockés en mémoire par colonnes (\emph{Column-Major Order}) : l'élément \texttt{A(i+1, j)} est adjacent en mémoire à \texttt{A(i, j)}. Par conséquent, toutes les boucles spatiales dans \texttt{COVOF\_2D-FORTRAN} sont ordonnées avec \textbf{$j$ en boucle externe et $i$ en boucle interne} :
\begin{verbatim}
do j = 1, mesh%ny
  do i = 1, mesh%nx
    fields%div(i, j) = (fields%u(i+1, j) - fields%u(i, j)) / mesh%dx &
                     + (fields%v(i, j+1) - fields%v(i, j)) / mesh%dy
  end do
end do
\end{verbatim}
Cet ordonnancement garantit exactement la même localité spatiale L1/L2 et la même efficacité SIMD que l'ordonnancement \emph{Row-Major} de C++17.
"""
    },
    {
        "num": 17,
        "title": "Communications Parallèles MPI et Échanges de Halos en Fortran 2018 (module_mpi_comm.f90)",
        "files": [("src/module_mpi_comm.f90", "fortran")],
        "theory": r"""
\section{Gestion de la Topologie Cartésienne et des Halos dans \texttt{module\_mpi\_comm.f90}}
Le module \texttt{src/module\_mpi\_comm.f90} (798 lignes) fournit le type dérivé \texttt{mpi\_context\_t} et les sous-routines d'échange de couronnes fantômes pour les trois topologies de grilles MAC :
\begin{enumerate}
  \item \texttt{exchange\_halo\_scalar(mpi\_ctx, mesh, phi)} : synchronise un tableau centré \texttt{phi(1-ng:nx+ng, 1-ng:ny+ng)} sur une épaisseur \texttt{ng} en échangeant d'abord les colonnes Est/Ouest puis les lignes Nord/Sud (propageant ainsi les coins diagonaux).
  \item \texttt{exchange\_halo\_u(mpi\_ctx, mesh, u)} et \texttt{exchange\_halo\_v(mpi\_ctx, mesh, v)} : gèrent le décalage d'une face sur les bords partagés entre deux sous-domaines MPI adjacents ainsi que la périodicité torique (\texttt{BC\_PERIODIC}).
  \item \texttt{mpi\_allreduce\_sum}, \texttt{mpi\_allreduce\_max}, \texttt{mpi\_allreduce\_min} : encapsulent les réductions collectives globales utilisées par les solveurs de Krylov et les bilans intégraux de conservation.
\end{enumerate}
"""
    },
    {
        "num": 18,
        "title": "Matrices Creuses CSR, Préconditionneurs et Solveurs de Krylov en Fortran 2018 (module_linear_solvers.f90)",
        "files": [("src/module_linear_solvers.f90", "fortran")],
        "theory": r"""
\section{Implémentation Vectorisée des Solveurs de Krylov dans \texttt{module\_linear\_solvers.f90}}
Le module \texttt{src/module\_linear\_solvers.f90} (892 lignes) implémente la structure de matrice creuse \texttt{csr\_matrix\_t} (à indiçage Fortran à base 1 : $\texttt{row\_ptr}(1) = 1$, $\texttt{row\_ptr}(n+1) = \text{nnz} + 1$) ainsi que les solveurs :
\begin{itemize}
  \item \texttt{solve\_bicgstab2} : implémentation complète de l'algorithme de Sleijpen \& Fokkema (1993) avec préconditionnement à droite ou à gauche (Jacobi, SSOR, ILU(0)) ;
  \item \texttt{solve\_bicgstab} : gradient biconjugué stabilisé standard de Van der Vorst ;
  \item \texttt{solve\_pcg} : gradient conjugué préconditionné pour l'opérateur symétrique défini positif du Laplacien (Cas 3.1) ;
  \item \texttt{solve\_gmres} : méthode de minimisation résiduelle généralisée avec sous-espace de Krylov redémarré et rotations de Givens.
\end{itemize}
"""
    },
    {
        "num": 19,
        "title": "Solveur Monolithique de Navier-Stokes par Lagrangien Augmenté en Fortran 2018 (module_navier_stokes.f90)",
        "files": [("src/module_navier_stokes.f90", "fortran")],
        "theory": r"""
\section{Architecture des 1540 Lignes du Module \texttt{module\_navier\_stokes.f90}}
Ce module constitue le noyau physique et mathématique de \texttt{COVOF\_2D-FORTRAN}. Il regroupe :
\begin{enumerate}
  \item \textbf{Mise à jour des propriétés physiques diphasiques (\texttt{update\_fluid\_properties})} : calcule en chaque maille $(i, j)$ la masse volumique $\rho_{i,j} = C_{i,j}\rho_1 + (1-C_{i,j})\rho_2$ et la viscosité dynamique $\mu_{i,j}$ (y compris la loi exponentielle spatiale du Cas 3.9).
  \item \textbf{Application des conditions aux limites hydrodynamiques (\texttt{apply\_velocity\_bc})} : impose les conditions \texttt{BC\_WALL} (adhérence par cellules fantômes antisymétriques), \texttt{BC\_SLIP} (glissement libre par cellules fantômes symétriques), \texttt{BC\_INLET}, \texttt{BC\_OUTLET}, \texttt{BC\_PERIODIC} et \texttt{BC\_JET\_INLET} (injection centrale de buse $V_0$ sur $|x - 0.5| \le D/2$ et évacuation compensatoire $v_{\text{out}} = |V_0| D / (L_x - D)$ garantissant $\oint_{\partial\Omega} \mathbf{u}\cdot\mathbf{n}\,d\Gamma = 0$).
  \item \textbf{Calcul de la force capillaire CSF Balanced-Force (\texttt{compute\_surface\_tension\_force})} : lisse la fonction couleur $C$ sur un voisinage $3\times 3$, calcule les normales unitaires aux sommets, en déduit la courbure centrée $\kappa_{i,j} = -\nabla_h \cdot \mathbf{n}$ et projette la force capillaire $\sigma \kappa_{\text{face}} \nabla_h C$ directement sur les faces MAC $u$ et $v$.
  \item \textbf{Assemblage et résolution du système monolithique du Lagrangien Augmenté (\texttt{solve\_momentum\_augmented\_lagrangian})} : assemble la matrice creuse $\mathbf{A}_{\text{aug}}$ à 9 bandes incluant l'inertie $\rho/\Delta t$, la convection décentrée, le tenseur visqueux complet à moyenne harmonique nodale, la pénalisation de Brinkman $\chi_s/\eta$ et le couplage croisé $-r_{al}\nabla_h(\nabla_h\cdot\mathbf{u})$, puis exécute les sous-itérations d'Uzawa $p^{k+1} = p^k - r_{al}\nabla_h\cdot\mathbf{u}^{k+1}$.
  \item \textbf{Résolution de l'équation de l'énergie thermique (\texttt{solve\_energy\_equation})} : intègre implicitement la température $T_{i,j}^{n+1}$ et applique les conditions aux limites thermiques isothermes et adiabatiques (Cas 6.1 de Vahl Davis).
\end{enumerate}
"""
    },
    {
        "num": 20,
        "title": "Reconstruction Géométrique d'Interface VOF en Fortran 2018 (module_interface_reconstruction.f90)",
        "files": [("src/module_interface_reconstruction.f90", "fortran")],
        "theory": r"""
\section{Algorithmes Géométriques dans \texttt{module\_interface\_reconstruction.f90}}
Le module \texttt{src/module\_interface\_reconstruction.f90} (569 lignes) implémente en Fortran 2018 :
\begin{itemize}
  \item \texttt{reconstruct\_interface\_normals} : sélectionne entre le filtre $3\times 3$ de Parker-Youngs (\texttt{RECON\_PLIC}), l'algorithme à 6 candidats exacts d'\textbf{ELVIRA} (\texttt{RECON\_ELVIRA}), l'optimisation quadratique non linéaire de \textbf{LVIRA} (\texttt{RECON\_LVIRA}) et l'alignement axial \textbf{SLIC} (\texttt{RECON\_SLIC}).
  \item \texttt{plic\_alpha\_from\_c} et \texttt{plic\_volume\_from\_alpha} : implémentent les formules analytiques directes et inverses de Scardovelli \& Zaleski (2000) reliant la constante de droite $\alpha$ à l'aire tronquée $C$ dans le carré unité $[0, 1]^2$.
  \item \texttt{plic\_flux\_1d} : calcule géométriquement le volume exact de phase 1 contenu dans la bande lagrangienne $[1 - \text{cfl}, 1] \times [0, 1]$ traversant une face de maille au cours du pas de temps $\Delta t$.
\end{itemize}
"""
    },
    {
        "num": 21,
        "title": "Advection Géométrique et Algébrique VOF en Fortran 2018 (module_interface_vof.f90)",
        "files": [("src/module_interface_vof.f90", "fortran")],
        "theory": r"""
\section{Architecture des 1034 Lignes de \texttt{module\_interface\_vof.f90}}
Ce module implémente la séparation directionnelle d'ordre 2 de Strang (\texttt{advect\_vof\_interface}) et les trois grandes familles de schémas de transport interfacial :
\begin{enumerate}
  \item \textbf{Schéma géométrique de Weymouth \& Yue (2010)} (\texttt{advect\_weymouth\_yue\_x}, \texttt{advect\_weymouth\_yue\_y}) : utilise la fonction indicatrice binaire $C_{\text{dil}}(i,j) = \mathbb{I}(C_{i,j}^n \ge 0.5)$ partagée entre les deux balayages directionnels pour garantir une conservation de la masse exacte à la précision machine.
  \item \textbf{Schéma géométrique CIAM de Li (1995)} (\texttt{advect\_ciam\_x}, \texttt{advect\_ciam\_y}) : transporte lagrangiennement les segments PLIC en tenant compte de l'étirement/compression affine de la bande donneuse par le gradient de vitesse local $\partial u / \partial x$.
  \item \textbf{Schémas algébriques TVD haute résolution} (\texttt{advect\_tvd\_x}, \texttt{advect\_tvd\_y}) : implémentent les limiteurs de pente non linéaires \textbf{Superbee}, \textbf{Van Leer}, \textbf{MUSCL / Monotonized Central (MC)}, \textbf{Minmod} et \textbf{HRIC}.
\end{enumerate}
"""
    },
    {
        "num": 22,
        "title": "Générateur VTK XML Parallèle en Fortran 2018 (module_vtk_writer.f90)",
        "files": [("src/module_vtk_writer.f90", "fortran")],
        "theory": r"""
\section{Écriture Conforme des Fichiers \texttt{.vtr} et \texttt{.pvtr} (\texttt{module\_vtk\_writer.f90})}
Le module \texttt{src/module\_vtk\_writer.f90} (500 lignes) écrit les fichiers XML VTK lisibles sous ParaView et Python. Une attention particulière est portée au formatage des flottants IEEE-754 : en Fortran, le descripteur standard \texttt{ES15.7} omet parfois la lettre \texttt{E} lorsque l'exposant comporte trois chiffres (par exemple \texttt{1.2345678-105} lorsque $|x| < 10^{-99}$). L'utilisation du descripteur explicite à trois chiffres d'exposant \texttt{ES16.8E3} dans \texttt{module\_vtk\_writer.f90} garantit que la lettre \texttt{E} est toujours présente, assurant une compatibilité universelle avec tous les parseurs XML VTK.
"""
    },
    {
        "num": 23,
        "title": "Boucle Temporelle Principale, Diagnostics Intégraux et Vérification en Fortran 2018 (module_simulation.f90)",
        "files": [("src/module_simulation.f90", "fortran")],
        "theory": r"""
\section{Orchestration Complète dans \texttt{module\_simulation.f90} (1004 Lignes)}
Le module \texttt{src/module\_simulation.f90} réalise l'intégration temporelle complète et l'auto-vérification des 16 cas de validation :
\begin{enumerate}
  \item \texttt{simulation\_init} : initialise le contexte MPI, le maillage MAC, les champs hydrodynamiques/thermiques/VOF, les masques de Brinkman et les conditions initiales sous-échantillonnées $8\times 8$.
  \item \texttt{simulation\_run} : aiguille vers le mode elliptique pur (\texttt{run\_laplacian\_test} pour le Cas 3.1), le mode cinématique de Rider-Kothe (\texttt{run\_vortex\_advection} pour les 12 variantes du Cas 3.19) ou la boucle complète Navier-Stokes + VOF + Énergie.
  \item \texttt{write\_diagnostics\_csv} : calcule à chaque pas de temps les intégrales physiques exactes (volume de phase, centre de masse $(X_c, Y_c)$, vitesse moyenne $(U_c, V_c)$, énergie cinétique diphasique $E_k(t) = \frac{1}{2}\iint \rho(C)\|\mathbf{u}\|^2\,d\Omega$, enstrophie $\mathcal{E}(t) = \frac{1}{2}\iint \omega_z^2\,d\Omega$, nombre de Nusselt $\overline{Nu}$ et saut de pression de Laplace $\Delta P$) et les exporte dans les fichiers CSV synchronisés avec C++17.
\end{enumerate}
"""
    }
]


LATEX_PREAMBLE = r"""\documentclass[10pt,a4paper,oneside]{book}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[french]{babel}
\usepackage{lmodern}
\usepackage[margin=1.8cm]{geometry}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{booktabs,longtable,array}
\usepackage{xcolor}
\usepackage{listings}
\usepackage{fancyhdr}
\usepackage{hyperref}

\definecolor{codebg}{rgb}{0.97,0.97,0.98}
\definecolor{codeframe}{rgb}{0.82,0.84,0.88}
\definecolor{kwcolor}{rgb}{0.0,0.28,0.65}
\definecolor{comcolor}{rgb}{0.25,0.50,0.35}
\definecolor{strcolor}{rgb}{0.65,0.15,0.15}

\lstset{
  backgroundcolor=\color{codebg},
  frame=single,
  rulecolor=\color{codeframe},
  basicstyle=\ttfamily\scriptsize,
  keywordstyle=\color{kwcolor}\bfseries,
  commentstyle=\color{comcolor}\itshape,
  stringstyle=\color{strcolor},
  numbers=left,
  numberstyle=\tiny\color{gray},
  stepnumber=1,
  numbersep=6pt,
  breaklines=true,
  breakatwhitespace=false,
  tabsize=2,
  showstringspaces=false,
  upquote=true,
  literate=
    {é}{{\'e}}1 {è}{{\`e}}1 {ê}{{\^e}}1 {ë}{{\"e}}1
    {à}{{\`a}}1 {â}{{\^a}}1 {ä}{{\"a}}1
    {î}{{\^i}}1 {ï}{{\"i}}1
    {ô}{{\^o}}1 {ö}{{\"o}}1
    {ù}{{\`u}}1 {û}{{\^u}}1 {ü}{{\"u}}1
    {ç}{{\c{c}}}1
    {É}{{\'E}}1 {È}{{\`E}}1 {Ê}{{\^E}}1
    {À}{{\`A}}1 {Â}{{\^A}}1
    {Î}{{\^I}}1 {Ô}{{\^O}}1 {Û}{{\^U}}1 {Ç}{{\c{C}}}1
    {œ}{{\oe}}1 {Œ}{{\OE}}1
    {°}{{\textdegree}}1
    {×}{{$\times$}}1 {÷}{{$\div$}}1 {±}{{$\pm$}}1
    {≤}{{$\le$}}1 {≥}{{$\ge$}}1 {≠}{{$\ne$}}1 {≈}{{$\approx$}}1
    {∞}{{$\infty$}}1 {∂}{{$\partial$}}1 {∇}{{$\nabla$}}1
    {α}{{$\alpha$}}1 {β}{{$\beta$}}1 {γ}{{$\gamma$}}1 {δ}{{$\delta$}}1
    {ε}{{$\varepsilon$}}1 {η}{{$\eta$}}1 {θ}{{$\theta$}}1 {κ}{{$\kappa$}}1
    {λ}{{$\lambda$}}1 {μ}{{$\mu$}}1 {ν}{{$\nu$}}1 {ξ}{{$\xi$}}1
    {π}{{$\pi$}}1 {ρ}{{$\rho$}}1 {σ}{{$\sigma$}}1 {τ}{{$\tau$}}1
    {φ}{{$\varphi$}}1 {χ}{{$\chi$}}1 {ψ}{{$\psi$}}1 {ω}{{$\omega$}}1
    {Δ}{{$\Delta$}}1 {Ω}{{$\Omega$}}1 {Ψ}{{$\Psi$}}1
    {—}{{---}}1 {–}{{--}}1 {’}{{'}}1 {‘}{{`}}1
}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\leftmark}
\fancyhead[R]{\thepage}
\renewcommand{\headrulewidth}{0.4pt}

\hypersetup{
  colorlinks=true,
  linkcolor=blue!70!black,
  urlcolor=blue!70!black,
  pdftitle={Guide du Développeur Exhaustif COVOF 2D (C++17 et Fortran 2018)},
  pdfauthor={Laboratoire TREFLE / I2M Bordeaux - Projet COVOF 2D}
}

\begin{document}

\begin{titlepage}
\centering
\vspace*{2cm}
{\Huge\bfseries GUIDE DU DÉVELOPPEUR EXHAUSTIF\\[0.4cm] SOLVEUR DIPHASIQUE ET THERMOCONVECTIF\\[0.4cm] COVOF 2D\par}
\vspace{1cm}
{\Large\bfseries Architecture Logicielle, Formulations Mathématiques Discrètes, Algorithmes Numériques et Analyse Intégrale des Codes Sources\par}
\vspace{1.2cm}
\rule{\textwidth}{1.5pt}\\[0.5cm]
{\LARGE\bfseries PARTIE I : Implémentation Orientée Objet C++17 (\texttt{COVOF\_2D-C++})\\[0.3cm]
PARTIE II : Implémentation Modulaire Fortran 2018 (\texttt{COVOF\_2D-FORTRAN})\par}
\vspace{0.3cm}
\rule{\textwidth}{1.5pt}
\vfill
{\large
\textbf{Méthodes Numériques Documentées :}\\
Volumes Finis sur Grille Décalée MAC $\bullet$ Lagrangien Augmenté Monolithique ($r_{al}$) et Algorithme d'Uzawa\\
Solveurs de Krylov Parallèles OpenMPI : BiCGSTAB(2), BiCGSTAB, PCG, GMRES($m$), MINRES\\
Capture d'Interface VOF : Reconstruction Géométrique PLIC, ELVIRA, LVIRA, SLIC\\
Advection Conservatrice : Weymouth \& Yue (2010), CIAM (Li 1995), TVD (Superbee, MC, Van Leer, Minmod, HRIC)\\
Tension Superficielle CSF Balanced-Force $\bullet$ Pénalisation Volumique de Brinkman $\bullet$ Thermoconvection de Boussinesq\\[1cm]
\textbf{Version 3.0.0 — Octobre 2026}
}
\end{titlepage}

\tableofcontents
\cleardoublepage
"""


def sanitize_for_lstlisting(code_str):
    # Replace box-drawing characters and exotic unicode symbols so listings renders cleanly
    replacements = {
        "─": "-", "━": "=", "│": "|", "┃": "|",
        "┌": "+", "┐": "+", "└": "+", "┘": "+",
        "├": "+", "┤": "+", "┬": "+", "┴": "+", "┼": "+",
        "╔": "+", "╗": "+", "╚": "+", "╝": "+", "═": "=", "║": "|",
        "▸": ">", "►": ">", "•": "*", "·": ".", "…": "...",
        "→": "->", "←": "<-", "↔": "<->", "⇒": "=>", "⇐": "<=",
        "∈": " in ", "∉": " not in ", "⊂": " subset ", "⊆": " subseteq ",
        "∑": "sum", "∏": "prod", "∫": "int", "√": "sqrt",
        "⁰": "^0", "¹": "^1", "²": "^2", "³": "^3", "⁴": "^4",
        "⁵": "^5", "⁶": "^6", "⁷": "^7", "⁸": "^8", "⁹": "^9",
        "⁺": "^+", "⁻": "^-", "ⁿ": "^n",
        "₀": "_0", "₁": "_1", "₂": "_2", "₃": "_3", "₄": "_4",
        "₅": "_5", "₆": "_6", "₇": "_7", "₈": "_8", "₉": "_9",
        "ᵢ": "_i", "ⱼ": "_j", "ₖ": "_k", "ₘ": "_m", "ₙ": "_n",
        "✅": "[OK]", "❌": "[ERR]", "⚠️": "[WARN]",
        "“": "\"", "”": "\"", "«": "<<", "»": ">>",
        "\t": "  "
    }
    for k, v in replacements.items():
        code_str = code_str.replace(k, v)
    # Strip any remaining non-latin1 characters except those handled in literate
    allowed_extra = set("éèêëàâäîïôöùûüçÉÈÊÀÂÎÔÛÇœŒ°×÷±≤≥≠≈∞∂∇αβγδεηθκλμνξπρστφχψωΔΩΨ—–’‘")
    out_chars = []
    for ch in code_str:
        if ord(ch) < 128 or ch in allowed_extra:
            out_chars.append(ch)
        else:
            out_chars.append("?")
    return "".join(out_chars)


def extract_cpp_symbols(code):
    """Extract classes, structs, and function/method definitions with line numbers."""
    import re
    items = []
    lines = code.splitlines()
    for idx, line in enumerate(lines, 1):
        s = line.strip()
        if s.startswith("//") or s.startswith("/*") or s.startswith("*"):
            continue
        m_cls = re.match(r"^(class|struct|enum\s+class)\s+([A-Za-z0-9_]+)", s)
        if m_cls:
            items.append((idx, m_cls.group(1), m_cls.group(2), "Déclaration de type / classe / structure"))
            continue
        m_fn = re.match(r"^([A-Za-z0-9_:<>&\*\s]+?)\s+([A-Za-z0-9_:]+)\s*\(([^;]*)\)\s*(const)?\s*\{?", s)
        if m_fn and not any(s.startswith(kw) for kw in ["if", "for", "while", "switch", "return", "else", "#"]):
            ret = m_fn.group(1).strip()
            name = m_fn.group(2).strip()
            if len(name) > 2 and not name.startswith("__"):
                items.append((idx, "Méthode / Fonction", name, f"Type retour : {ret}"))
    return items[:45]


def extract_f90_symbols(code):
    """Extract modules, derived types, subroutines, and functions with line numbers."""
    import re
    items = []
    lines = code.splitlines()
    for idx, line in enumerate(lines, 1):
        s = line.strip()
        if s.startswith("!"):
            continue
        m_mod = re.match(r"^module\s+([A-Za-z0-9_]+)", s, re.IGNORECASE)
        if m_mod and "procedure" not in s.lower():
            items.append((idx, "module", m_mod.group(1), "Unité modulaire Fortran 2018"))
            continue
        m_typ = re.match(r"^type\s*(?:,\s*public\s*)?::\s*([A-Za-z0-9_]+)", s, re.IGNORECASE)
        if m_typ:
            items.append((idx, "type dérivé", m_typ.group(1), "Structure de données Fortran 2018"))
            continue
        m_sub = re.match(r"^(?:pure\s+|elemental\s+|recursive\s+)?(subroutine|function)\s+([A-Za-z0-9_]+)\s*\(([^)]*)\)", s, re.IGNORECASE)
        if m_sub:
            kind = m_sub.group(1).lower()
            name = m_sub.group(2)
            args = m_sub.group(3).strip()
            if len(args) > 42:
                args = args[:39] + "..."
            items.append((idx, kind, name, f"Arguments : ({args})"))
    return items[:45]


def tex_escape(s):
    return (s.replace("\\", "\\textbackslash ")
             .replace("_", "\\_")
             .replace("&", "\\&")
             .replace("%", "\\%")
             .replace("#", "\\#")
             .replace("<", "$<$")
             .replace(">", "$>$"))


def build_developer_guide():
    tex_parts = [LATEX_PREAMBLE]
    md_parts = [
        "# Guide du Développeur Exhaustif — Solveur COVOF 2D (C++17 & Fortran 2018)\n"
        "## Partie I : Architecture et Code Source C++17 (`COVOF_2D-C++`)\n"
        "## Partie II : Architecture et Code Source Fortran 2018 (`COVOF_2D-FORTRAN`)\n\n---\n"
    ]

    # PART I : C++17
    tex_parts.append(r"\part{Guide du Développeur C++17 (\texttt{COVOF\_2D-C++})}" + "\n")
    md_parts.append("\n# PARTIE I — Guide du Développeur C++17 (`COVOF_2D-C++`)\n\n")

    for ch in CPP_CHAPTERS:
        tex_parts.append(f"\\chapter{{{ch['title']}}}\n")
        tex_parts.append(ch["theory"].strip() + "\n\n")

        md_parts.append(f"\n## Chapitre {ch['num']} : {ch['title']}\n\n")
        md_parts.append(ch["theory"].strip() + "\n\n")

        for rel_path, lang in ch["files"]:
            code = read_source(ROOT_CPP, rel_path)
            clean_code = sanitize_for_lstlisting(code)
            escaped_path = rel_path.replace("_", "\\_")
            lst_lang = "[11]C++" if lang == "cpp" else "bash"

            if lang == "cpp":
                syms = extract_cpp_symbols(code)
                if syms:
                    tex_parts.append(f"\\section{{Répertoire des Classes, Structures et Méthodes de \\texttt{{{escaped_path}}}}}\n")
                    tex_parts.append("\\begin{center}\\small\n\\begin{tabular}{clll}\n\\toprule\n")
                    tex_parts.append("\\textbf{Ligne} & \\textbf{Catégorie} & \\textbf{Identifiant C++17} & \\textbf{Signature / Rôle} \\\\\n\\midrule\n")
                    md_parts.append(f"### Répertoire des Symboles de `{rel_path}`\n\n| Ligne | Catégorie | Identifiant C++17 | Signature / Rôle |\n|:---:|:---|:---|:---|\n")
                    for ln, cat, ident, desc in syms:
                        tex_parts.append(f"{ln} & {tex_escape(cat)} & \\texttt{{{tex_escape(ident)}}} & {tex_escape(desc)} \\\\\n")
                        md_parts.append(f"| {ln} | {cat} | `{ident}` | {desc} |\n")
                    tex_parts.append("\\bottomrule\n\\end{tabular}\n\\end{center}\n\n")
                    md_parts.append("\n")

            tex_parts.append(f"\\section{{Code Source Intégral Annoté : \\texttt{{{escaped_path}}}}}\n")
            tex_parts.append(
                f"Le listing ci-dessous présente l'intégralité des {len(code.splitlines())} lignes du fichier "
                f"\\texttt{{{escaped_path}}} de l'implémentation C++17 :\n\n"
            )
            tex_parts.append(f"\\begin{{lstlisting}}[language={{{lst_lang}}}]\n{clean_code}\n\\end{{lstlisting}}\n\n")

            md_parts.append(f"### Code Source Intégral : `{rel_path}` ({len(code.splitlines())} lignes)\n\n")
            md_parts.append(f"```{lang}\n{code}\n```\n\n")

    # PART II : FORTRAN 2018
    tex_parts.append(r"\part{Guide du Développeur Fortran 2018 (\texttt{COVOF\_2D-FORTRAN})}" + "\n")
    md_parts.append("\n---\n\n# PARTIE II — Guide du Développeur Fortran 2018 (`COVOF_2D-FORTRAN`)\n\n")

    for ch in F90_CHAPTERS:
        tex_parts.append(f"\\chapter{{{ch['title']}}}\n")
        tex_parts.append(ch["theory"].strip() + "\n\n")

        md_parts.append(f"\n## Chapitre {ch['num']} : {ch['title']}\n\n")
        md_parts.append(ch["theory"].strip() + "\n\n")

        for rel_path, lang in ch["files"]:
            code = read_source(ROOT_F90, rel_path)
            clean_code = sanitize_for_lstlisting(code)
            escaped_path = rel_path.replace("_", "\\_")
            lst_lang = "[95]Fortran" if lang == "fortran" else "bash"

            if lang == "fortran":
                syms = extract_f90_symbols(code)
                if syms:
                    tex_parts.append(f"\\section{{Répertoire des Types Dérivés, Sous-Routines et Fonctions de \\texttt{{{escaped_path}}}}}\n")
                    tex_parts.append("\\begin{center}\\small\n\\begin{tabular}{clll}\n\\toprule\n")
                    tex_parts.append("\\textbf{Ligne} & \\textbf{Nature} & \\textbf{Identifiant Fortran 2018} & \\textbf{Interface / Arguments} \\\\\n\\midrule\n")
                    md_parts.append(f"### Répertoire des Procédures de `{rel_path}`\n\n| Ligne | Nature | Identifiant Fortran 2018 | Interface / Arguments |\n|:---:|:---|:---|:---|\n")
                    for ln, cat, ident, desc in syms:
                        tex_parts.append(f"{ln} & {tex_escape(cat)} & \\texttt{{{tex_escape(ident)}}} & {tex_escape(desc)} \\\\\n")
                        md_parts.append(f"| {ln} | {cat} | `{ident}` | {desc} |\n")
                    tex_parts.append("\\bottomrule\n\\end{tabular}\n\\end{center}\n\n")
                    md_parts.append("\n")

            tex_parts.append(f"\\section{{Code Source Intégral Annoté : \\texttt{{{escaped_path}}}}}\n")
            tex_parts.append(
                f"Le listing ci-dessous présente l'intégralité des {len(code.splitlines())} lignes du fichier "
                f"\\texttt{{{escaped_path}}} de l'implémentation Fortran 2018 :\n\n"
            )
            tex_parts.append(f"\\begin{{lstlisting}}[language={{{lst_lang}}}]\n{clean_code}\n\\end{{lstlisting}}\n\n")

            md_parts.append(f"### Code Source Intégral : `{rel_path}` ({len(code.splitlines())} lignes)\n\n")
            md_parts.append(f"```{lang}\n{code}\n```\n\n")

    tex_parts.append(r"\end{document}" + "\n")

    full_tex = "".join(tex_parts)
    full_md = "".join(md_parts)

    for root in [ROOT_CPP, ROOT_F90]:
        docs_dir = os.path.join(root, "docs")
        tex_path = os.path.join(docs_dir, "DEVELOPER_GUIDE.tex")
        md_path = os.path.join(docs_dir, "DEVELOPER_GUIDE.md")
        with open(tex_path, "w", encoding="utf-8") as f:
            f.write(full_tex)
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(full_md)
        print(f"[OK] Écrit : {tex_path} ({len(full_tex)} octets) et {md_path} ({len(full_md)} octets)")

    # Compile with pdflatex (2 passes for TOC and cross-references)
    for root in [ROOT_CPP, ROOT_F90]:
        docs_dir = os.path.join(root, "docs")
        for pass_num in [1, 2]:
            subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", "DEVELOPER_GUIDE.tex"],
                cwd=docs_dir,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        pdf_path = os.path.join(docs_dir, "DEVELOPER_GUIDE.pdf")
        log_path = os.path.join(docs_dir, "DEVELOPER_GUIDE.log")
        with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
            log_txt = f.read()
        import re
        m = re.search(r"Output written on DEVELOPER_GUIDE\.pdf \((\d+) pages", log_txt)
        pages = m.group(1) if m else "unknown"
        print(f"[OK] Compilé : {pdf_path} -> {pages} pages ({os.path.getsize(pdf_path)} octets)")


if __name__ == "__main__":
    build_developer_guide()
