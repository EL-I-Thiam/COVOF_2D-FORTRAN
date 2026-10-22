!===============================================================================
! Module : module_constants
! Description : Définition des constantes physiques fondamentales, paramètres
!               numériques seuils, énumérations de conditions aux limites, types
!               de solveurs et schémas convectifs pour le code COVOF 2D.
!===============================================================================
module module_constants
  ! Importation du paramètre de précision double standardisé
  use module_precision, only: dp
  ! Obligation stricte de typage explicite pour toutes les constantes
  implicit none
  ! Encapsulation : seules les constantes marquées 'public' sont exportées
  private

  ! --- 1. Constantes Mathématiques et Physiques Fondamentales ---
  ! Valeur exacte du nombre Pi (35 décimales significatives en précision double)
  real(dp), parameter, public :: PI = 3.141592653589793238462643383279502884197_dp
  ! Seuil numérique infinitésimal pour prévenir les divisions par zéro dans les gradients
  real(dp), parameter, public :: EPS_SMALL = 1.0e-30_dp
  ! Perméabilité de Darcy par défaut dans les obstacles solides pénalisés (loi de Brinkman)
  real(dp), parameter, public :: PERMEABILITY_SOLID_DEFAULT = 1.0e-40_dp

  ! --- 2. Conditions aux Limites Hydrodynamiques (Vitesse) ---
  ! Vitesse imposée (Adhérence u=0, paroi mobile de cavité, ou profil d'entrée de Poiseuille)
  integer, parameter, public :: BC_DIRICHLET = 1
  ! Sortie libre convective ou gradient normal nul (du/dn = 0, dv/dn = 0)
  integer, parameter, public :: BC_NEUMANN   = 2
  ! Raccordement spatial périodique (face Est reliée à la face Ouest, Nord à Sud)
  integer, parameter, public :: BC_PERIODIC  = 3
  ! Glissement parfait sans frottement tangentiel (vitesse normale nulle, cisaillement nul)
  integer, parameter, public :: BC_FREESLIP  = 4

  ! --- 3. Conditions aux Limites Thermiques (Énergie) ---
  ! Température pariétale imposée (Paroi isotherme chaude ou froide de Vahl Davis)
  integer, parameter, public :: THERMAL_DIRICHLET = 1
  ! Flux thermique imposé ou paroi adiabatique sans transfert de chaleur (dT/dn = 0)
  integer, parameter, public :: THERMAL_NEUMANN   = 2

  ! --- 4. Formulations du Couplage Vitesse-Pression Navier-Stokes ---
  ! Méthode du Lagrangien Augmenté avec pénalité de divergence et actualisation d'Uzawa
  integer, parameter, public :: SOLVER_AUGMENTED_LAGRANGIAN = 1
  ! Formulation couplée monolithique directe avec intégration du gradient de pression
  integer, parameter, public :: SOLVER_MONOLITHIC           = 2

  ! --- 5. Algorithmes de Résolution Linéaire de Krylov (Matrix-Free) ---
  ! Biconjugué Gradient Stabilisé d'ordre 2 (Sleijpen & van der Vorst) pour systèmes asymétriques
  integer, parameter, public :: LINEAR_BICGSTAB2 = 1
  ! Biconjugué Gradient Stabilisé standard (van der Vorst 1992)
  integer, parameter, public :: LINEAR_BICGSTAB  = 2
  ! Gradient Conjugué Préconditionné pour opérateurs symétriques définis positifs (Poisson)
  integer, parameter, public :: LINEAR_CG        = 3
  ! Generalized Minimal Residual (Saad & Schultz 1986) avec orthogonalisation de Krylov
  integer, parameter, public :: LINEAR_GMRES     = 4
  ! Solveur direct bande ou factorisation LU locale
  integer, parameter, public :: LINEAR_DIRECT    = 5

  ! --- 6. Préconditionneurs Linéaires Parallèles ---
  ! Aucun préconditionnement (résolution brute de l'opérateur non préconditionné)
  integer, parameter, public :: PRECOND_NONE   = 0
  ! Préconditionneur diagonal par blocs de Jacobi (inversion de la diagonale locale)
  integer, parameter, public :: PRECOND_JACOBI = 1
  ! Symmetric Successive Over-Relaxation (SSOR) pour lissage spectral rapide
  integer, parameter, public :: PRECOND_SSOR   = 2
  ! Factorisation incomplète Incomplete LU sans remplissage (ILU(0))
  integer, parameter, public :: PRECOND_ILU0   = 3

  ! --- 7. Schémas de Discrétisation Spatiale des Termes Convectifs ---
  ! Schéma décentré amont d'ordre 1 (Upwind) inconditionnellement stable mais diffusif
  integer, parameter, public :: CONV_UPWIND       = 1
  ! Schéma centré d'ordre 2 (Centered) sans dissipation artificielle
  integer, parameter, public :: CONV_CENTERED     = 2
  ! Schéma TVD Superbee à haute résolution avec limiteur de flux non linéaire (Roe 1986)
  integer, parameter, public :: CONV_TVD_SUPERBEE = 3

  ! --- 8. Modèles de Viscosité Dynamique du Fluide ---
  ! Fluide newtonien homogène à viscosité dynamique constante mu = mu0
  integer, parameter, public :: VISC_CONSTANT    = 1
  ! Fluide à viscosité variable spatiale exponentielle mu(y) = mu0 * exp(a * (y-e)^2 / e^2)
  integer, parameter, public :: VISC_EXPONENTIAL = 2

  ! --- 9. Types de Configurations Initiales et Cas de Vérification ---
  ! Écoulement monophasique uniforme avec fluide au repos (u=0, v=0, C=1)
  integer, parameter, public :: IC_UNIFORM          = 1
  ! Profil parabolique d'entrée de Poiseuille établi
  integer, parameter, public :: IC_POISEUILLE       = 2
  ! Tâche circulaire d'interface transportée dans un vortex cisaillé réversible
  integer, parameter, public :: IC_VORTEX_SHEAR     = 3
  ! Instabilité diphasique gravitationnelle de Rayleigh-Taylor avec interface perturbée
  integer, parameter, public :: IC_RAYLEIGH_TAYLOR  = 4
  ! Gouttelette liquide statique immergée pour la vérification de la loi de Laplace
  integer, parameter, public :: IC_DROPLET          = 5
  ! Bulle de gaz ascendante déformable dans un liquide visqueux (Benchmark de Hysing)
  integer, parameter, public :: IC_BUBBLE           = 6
  ! Effondrement d'une colonne d'eau sous gravité (Rupture de barrage de Martin & Moyce)
  integer, parameter, public :: IC_DAM_BREAK        = 7
  ! Test analytique d'ordre de convergence de l'opérateur Laplacien discret à 5 points
  integer, parameter, public :: IC_LAPLACIAN_TEST   = 8

  ! --- 10. Formes Géométriques d'Obstacles Immergés (Brinkman) ---
  ! Aucun obstacle immergé dans le domaine fluide
  integer, parameter, public :: OBS_NONE     = 0
  ! Obstacle rectangulaire défini par ses coordonnées minimales et maximales (x0, y0, x1, y1)
  integer, parameter, public :: OBS_BOX      = 1
  ! Obstacle circulaire cylindrique défini par son centre (xc, yc) et son rayon R
  integer, parameter, public :: OBS_CYLINDER = 2

end module module_constants
