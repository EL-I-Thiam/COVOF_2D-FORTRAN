! [Procédure/Module] module module_json
module module_json
  use module_precision, only: dp
  use module_constants, only: OBS_NONE, OBS_BOX, OBS_CYLINDER
  implicit none
  private

  public :: json_read_file
  public :: json_get_string
  public :: json_get_real
  public :: json_get_int
  public :: json_get_bool
  public :: json_get_obstacle

contains

  ! Lit un fichier complet dans une chaîne de caractères allouable en ignorant les commentaires //
  subroutine json_read_file(filename, buffer, ierr)
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: filename
    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable, intent(out) :: buffer
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(out) :: ierr

    ! Déclaration de variable locale / paramètre de calcul
    integer :: u, io
    ! Déclaration de variable locale / paramètre de calcul
    character(len=4096) :: line
    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable :: temp
    ! Déclaration de variable locale / paramètre de calcul
    integer :: comment_pos

    ! Évaluation de la norme du résidu ou de l'erreur numérique
    ierr = 0
    ! Affectation et calcul numérique
    buffer = ""
    ! Affectation et calcul numérique
    open(newunit=u, file=trim(filename), status='old', action='read', iostat=io)
    ! Évaluation de la condition logique
    if (io /= 0) then
       ! Évaluation de la norme du résidu ou de l'erreur numérique
       ierr = io
       ! Retour à la procédure appelante
       return
    end if

    do
       ! Affectation et calcul numérique
       read(u, '(A)', iostat=io) line
       ! Affectation et calcul numérique
       if (io /= 0) exit
       ! Affectation et calcul numérique
       comment_pos = index(line, '//')
       ! Affectation et calcul numérique
       if (comment_pos > 0) line = line(1:comment_pos-1)
       ! Affectation et calcul numérique
       line = adjustl(line)
       ! Évaluation de la condition logique
       if (len_trim(line) > 0) then
          ! Évaluation de la condition logique
          if (len(buffer) == 0) then
             ! Affectation et calcul numérique
             buffer = trim(line)
          ! Branche par défaut si condition non satisfaite
          else
             ! Affectation et calcul numérique
             temp = buffer // " " // trim(line)
             ! Affectation et calcul numérique
             buffer = temp
          end if
       end if
    end do
    close(u)
  end subroutine json_read_file

  ! Recherche d'un chemin séparé par des points dans le buffer JSON
  subroutine json_find_value_str(buffer, path, val_str, found)
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: buffer
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: path
    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable, intent(out) :: val_str
    ! Déclaration de variable locale / paramètre de calcul
    logical, intent(out) :: found

    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable :: cur_buf
    ! Déclaration de variable locale / paramètre de calcul
    integer :: p_dot, pos_key, pos_colon
    ! Déclaration de variable locale / paramètre de calcul
    character(len=64) :: key, sub_path
    ! Déclaration de variable locale / paramètre de calcul
    character(len=70) :: search_key
    ! Déclaration de variable locale / paramètre de calcul
    integer :: brace_count, start_pos, end_pos

    ! Affectation et calcul numérique
    found = .false.
    ! Affectation et calcul numérique
    val_str = ""
    ! Affectation et calcul numérique
    cur_buf = buffer

    ! Affectation et calcul numérique
    sub_path = adjustl(trim(path))

    ! Boucle d'itération numérique
    do while (len_trim(sub_path) > 0)
       ! Affectation et calcul numérique
       p_dot = index(sub_path, ".")
       ! Évaluation de la condition logique
       if (p_dot > 0) then
          ! Affectation et calcul numérique
          key = sub_path(1:p_dot-1)
          ! Affectation et calcul numérique
          sub_path = adjustl(sub_path(p_dot+1:))
       ! Branche par défaut si condition non satisfaite
       else
          ! Affectation et calcul numérique
          key = trim(sub_path)
          ! Affectation et calcul numérique
          sub_path = ""
       end if

       ! Affectation et calcul numérique
       search_key = '"' // trim(key) // '"'
       ! Affectation et calcul numérique
       pos_key = index(cur_buf, trim(search_key))
       ! Évaluation de la condition logique
       if (pos_key <= 0) then
          ! Affectation et calcul numérique
          found = .false.
          ! Retour à la procédure appelante
          return
       end if

       ! Trouver les deux-points ':' après la clé
       pos_colon = pos_key + len_trim(search_key)
       ! Boucle d'itération numérique
       do while (pos_colon <= len(cur_buf) .and. cur_buf(pos_colon:pos_colon) /= ':')
          ! Affectation et calcul numérique
          pos_colon = pos_colon + 1
       end do
       ! Évaluation de la condition logique
       if (pos_colon > len(cur_buf)) then
          ! Affectation et calcul numérique
          found = .false.
          ! Retour à la procédure appelante
          return
       end if

       ! Si ce n'est pas le dernier segment de chemin, chercher l'objet { ... } ou [ ... ] correspondant
       if (len_trim(sub_path) > 0) then
          ! Affectation et calcul numérique
          start_pos = pos_colon + 1
          ! Boucle d'itération numérique
          do while (start_pos <= len(cur_buf) .and. cur_buf(start_pos:start_pos) /= '{' .and. cur_buf(start_pos:start_pos) /= '[')
             ! Affectation et calcul numérique
             start_pos = start_pos + 1
          end do
          ! Évaluation de la condition logique
          if (start_pos > len(cur_buf)) then
             ! Affectation et calcul numérique
             found = .false.
             ! Retour à la procédure appelante
             return
          end if

          ! Affectation et calcul numérique
          brace_count = 1
          ! Affectation et calcul numérique
          end_pos = start_pos + 1
          ! Boucle d'itération numérique
          do while (end_pos <= len(cur_buf) .and. brace_count > 0)
             ! Affectation et calcul numérique
             if (cur_buf(end_pos:end_pos) == '{' .or. cur_buf(end_pos:end_pos) == '[') brace_count = brace_count + 1
             ! Affectation et calcul numérique
             if (cur_buf(end_pos:end_pos) == '}' .or. cur_buf(end_pos:end_pos) == ']') brace_count = brace_count - 1
             ! Affectation et calcul numérique
             end_pos = end_pos + 1
          end do
          ! Affectation et calcul numérique
          cur_buf = cur_buf(start_pos:end_pos-1)
       ! Branche par défaut si condition non satisfaite
       else
          ! Dernier élément du chemin : extraire la valeur
          start_pos = pos_colon + 1
          ! Boucle d'itération numérique
          do while (start_pos <= len(cur_buf) .and. (cur_buf(start_pos:start_pos) == ' ' .or. cur_buf(start_pos:start_pos) == char(9)))
             ! Affectation et calcul numérique
             start_pos = start_pos + 1
          end do
          ! Évaluation de la condition logique
          if (start_pos > len(cur_buf)) then
             ! Affectation et calcul numérique
             found = .false.
             ! Retour à la procédure appelante
             return
          end if

          ! Évaluation de la condition logique
          if (cur_buf(start_pos:start_pos) == '"') then
             ! Chaîne entre guillemets
             end_pos = start_pos + 1
             ! Boucle d'itération numérique
             do while (end_pos <= len(cur_buf) .and. cur_buf(end_pos:end_pos) /= '"')
                ! Affectation et calcul numérique
                end_pos = end_pos + 1
             end do
             ! Affectation et calcul numérique
             val_str = cur_buf(start_pos+1:end_pos-1)
             ! Affectation et calcul numérique
             found = .true.
             ! Retour à la procédure appelante
             return
          ! Branche par défaut si condition non satisfaite
          else
             ! Nombre, booléen, etc.
             end_pos = start_pos
             ! Boucle d'itération numérique
             do while (end_pos <= len(cur_buf) .and. &
                       ! Affectation et calcul numérique
                       cur_buf(end_pos:end_pos) /= ',' .and. &
                       ! Affectation et calcul numérique
                       cur_buf(end_pos:end_pos) /= '}' .and. &
                       ! Affectation et calcul numérique
                       cur_buf(end_pos:end_pos) /= ']' .and. &
                       ! Affectation et calcul numérique
                       cur_buf(end_pos:end_pos) /= ' ' .and. &
                       ! Affectation et calcul numérique
                       cur_buf(end_pos:end_pos) /= char(9) .and. &
                       ! Affectation et calcul numérique
                       cur_buf(end_pos:end_pos) /= char(10) .and. &
                       ! Affectation et calcul numérique
                       cur_buf(end_pos:end_pos) /= char(13))
                ! Affectation et calcul numérique
                end_pos = end_pos + 1
             end do
             ! Affectation et calcul numérique
             val_str = trim(adjustl(cur_buf(start_pos:end_pos-1)))
             ! Affectation et calcul numérique
             found = .true.
             ! Retour à la procédure appelante
             return
          end if
       end if
    end do
  end subroutine json_find_value_str

  ! [Procédure/Module] subroutine json_get_string
  subroutine json_get_string(buffer, path, val, default_val)
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: buffer, path, default_val
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(out) :: val
    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable :: s
    ! Déclaration de variable locale / paramètre de calcul
    logical :: found

    ! Appel de sous-programme modulaire
    call json_find_value_str(buffer, path, s, found)
    ! Évaluation de la condition logique
    if (found .and. allocated(s)) then
       ! Affectation et calcul numérique
       val = trim(s)
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       val = trim(default_val)
    end if
  end subroutine json_get_string

  ! [Procédure/Module] subroutine json_get_real
  subroutine json_get_real(buffer, path, val, default_val)
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: buffer, path
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(out) :: val
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(in) :: default_val
    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable :: s
    ! Déclaration de variable locale / paramètre de calcul
    logical :: found
    ! Déclaration de variable locale / paramètre de calcul
    integer :: io

    ! Appel de sous-programme modulaire
    call json_find_value_str(buffer, path, s, found)
    ! Évaluation de la condition logique
    if (found .and. allocated(s)) then
       ! Affectation et calcul numérique
       read(s, *, iostat=io) val
       ! Affectation et calcul numérique
       if (io /= 0) val = default_val
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       val = default_val
    end if
  end subroutine json_get_real

  ! [Procédure/Module] subroutine json_get_int
  subroutine json_get_int(buffer, path, val, default_val)
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: buffer, path
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(out) :: val
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(in) :: default_val
    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable :: s
    ! Déclaration de variable locale / paramètre de calcul
    logical :: found
    ! Déclaration de variable locale / paramètre de calcul
    integer :: io

    ! Appel de sous-programme modulaire
    call json_find_value_str(buffer, path, s, found)
    ! Évaluation de la condition logique
    if (found .and. allocated(s)) then
       ! Affectation et calcul numérique
       read(s, *, iostat=io) val
       ! Affectation et calcul numérique
       if (io /= 0) val = default_val
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       val = default_val
    end if
  end subroutine json_get_int

  ! [Procédure/Module] subroutine json_get_bool
  subroutine json_get_bool(buffer, path, val, default_val)
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: buffer, path
    ! Déclaration de variable locale / paramètre de calcul
    logical, intent(out) :: val
    ! Déclaration de variable locale / paramètre de calcul
    logical, intent(in) :: default_val
    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable :: s
    ! Déclaration de variable locale / paramètre de calcul
    logical :: found

    ! Appel de sous-programme modulaire
    call json_find_value_str(buffer, path, s, found)
    ! Évaluation de la condition logique
    if (found .and. allocated(s)) then
       ! Évaluation de la condition logique
       if (index(s, "true") > 0 .or. index(s, "TRUE") > 0 .or. index(s, "1") > 0) then
          ! Affectation et calcul numérique
          val = .true.
       ! Branche par défaut si condition non satisfaite
       else
          ! Affectation et calcul numérique
          val = .false.
       end if
    ! Branche par défaut si condition non satisfaite
    else
       ! Affectation et calcul numérique
       val = default_val
    end if
  end subroutine json_get_bool

  ! Détection et extraction des obstacles (BOX ou CYLINDER)
  subroutine json_get_obstacle(buffer, obs_type, x0, y0, x1, y1, xc, yc, radius, perm, found)
    ! Déclaration de variable locale / paramètre de calcul
    character(len=*), intent(in) :: buffer
    ! Déclaration de variable locale / paramètre de calcul
    integer, intent(out) :: obs_type
    ! Déclaration de variable locale / paramètre de calcul
    real(dp), intent(out) :: x0, y0, x1, y1, xc, yc, radius, perm
    ! Déclaration de variable locale / paramètre de calcul
    logical, intent(out) :: found

    ! Déclaration de variable locale / paramètre de calcul
    integer :: pos_obs, pos_end
    ! Déclaration de variable locale / paramètre de calcul
    character(len=:), allocatable :: obs_str
    ! Déclaration de variable locale / paramètre de calcul
    character(len=32) :: type_str

    ! Affectation et calcul numérique
    obs_type = OBS_NONE
    ! Affectation et calcul numérique
    x0 = 0.0_dp; y0 = 0.0_dp; x1 = 0.0_dp; y1 = 0.0_dp
    ! Affectation et calcul numérique
    xc = 0.0_dp; yc = 0.0_dp; radius = 0.0_dp; perm = 1.0e-40_dp
    ! Affectation et calcul numérique
    found = .false.

    ! Affectation et calcul numérique
    pos_obs = index(buffer, '"obstacles"')
    ! Affectation et calcul numérique
    if (pos_obs <= 0) return

    ! Affectation et calcul numérique
    pos_end = pos_obs + 11
    ! Boucle d'itération numérique
    do while (pos_end <= len(buffer) .and. buffer(pos_end:pos_end) /= ']')
       ! Affectation et calcul numérique
       pos_end = pos_end + 1
    end do
    if (pos_end > len(buffer)) return

    ! Affectation et calcul numérique
    obs_str = buffer(pos_obs:pos_end)
    ! Appel de sous-programme modulaire
    call json_get_string(obs_str, "type", type_str, "")

    ! Évaluation de la condition logique
    if (trim(type_str) == "BOX") then
       ! Affectation et calcul numérique
       obs_type = OBS_BOX
       ! Appel de sous-programme modulaire
       call json_get_real(obs_str, "x0", x0, 0.0_dp)
       ! Appel de sous-programme modulaire
       call json_get_real(obs_str, "y0", y0, 0.0_dp)
       ! Appel de sous-programme modulaire
       call json_get_real(obs_str, "x1", x1, 0.0_dp)
       ! Appel de sous-programme modulaire
       call json_get_real(obs_str, "y1", y1, 0.0_dp)
       ! Appel de sous-programme modulaire
       call json_get_real(obs_str, "permeability", perm, 1.0e-40_dp)
       ! Affectation et calcul numérique
       found = .true.
    ! Branche conditionnelle alternative
    else if (trim(type_str) == "CYLINDER") then
       ! Affectation et calcul numérique
       obs_type = OBS_CYLINDER
       ! Appel de sous-programme modulaire
       call json_get_real(obs_str, "xc", xc, 0.0_dp)
       ! Appel de sous-programme modulaire
       call json_get_real(obs_str, "yc", yc, 0.0_dp)
       ! Appel de sous-programme modulaire
       call json_get_real(obs_str, "radius", radius, 0.0_dp)
       ! Appel de sous-programme modulaire
       call json_get_real(obs_str, "permeability", perm, 1.0e-40_dp)
       ! Affectation et calcul numérique
       found = .true.
    end if
  end subroutine json_get_obstacle

end module module_json
