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
    character(len=*), intent(in) :: filename
    character(len=:), allocatable, intent(out) :: buffer
    integer, intent(out) :: ierr

    integer :: u, io
    character(len=4096) :: line
    character(len=:), allocatable :: temp
    integer :: comment_pos

    ierr = 0
    buffer = ""
    open(newunit=u, file=trim(filename), status='old', action='read', iostat=io)
    if (io /= 0) then
       ierr = io
       return
    end if

    do
       read(u, '(A)', iostat=io) line
       if (io /= 0) exit
       comment_pos = index(line, '//')
       if (comment_pos > 0) line = line(1:comment_pos-1)
       line = adjustl(line)
       if (len_trim(line) > 0) then
          if (len(buffer) == 0) then
             buffer = trim(line)
          else
             temp = buffer // " " // trim(line)
             buffer = temp
          end if
       end if
    end do
    close(u)
  end subroutine json_read_file

  ! Recherche d'un chemin séparé par des points dans le buffer JSON
  subroutine json_find_value_str(buffer, path, val_str, found)
    character(len=*), intent(in) :: buffer
    character(len=*), intent(in) :: path
    character(len=:), allocatable, intent(out) :: val_str
    logical, intent(out) :: found

    character(len=:), allocatable :: cur_buf
    integer :: p_dot, pos_key, pos_colon
    character(len=64) :: key, sub_path
    character(len=70) :: search_key
    integer :: brace_count, start_pos, end_pos

    found = .false.
    val_str = ""
    cur_buf = buffer

    sub_path = adjustl(trim(path))

    do while (len_trim(sub_path) > 0)
       p_dot = index(sub_path, ".")
       if (p_dot > 0) then
          key = sub_path(1:p_dot-1)
          sub_path = adjustl(sub_path(p_dot+1:))
       else
          key = trim(sub_path)
          sub_path = ""
       end if

       search_key = '"' // trim(key) // '"'
       pos_key = index(cur_buf, trim(search_key))
       if (pos_key <= 0) then
          found = .false.
          return
       end if

       ! Trouver les deux-points ':' après la clé
       pos_colon = pos_key + len_trim(search_key)
       do while (pos_colon <= len(cur_buf) .and. cur_buf(pos_colon:pos_colon) /= ':')
          pos_colon = pos_colon + 1
       end do
       if (pos_colon > len(cur_buf)) then
          found = .false.
          return
       end if

       ! Si ce n'est pas le dernier segment de chemin, chercher l'objet { ... } ou [ ... ] correspondant
       if (len_trim(sub_path) > 0) then
          start_pos = pos_colon + 1
          do while (start_pos <= len(cur_buf) .and. cur_buf(start_pos:start_pos) /= '{' .and. cur_buf(start_pos:start_pos) /= '[')
             start_pos = start_pos + 1
          end do
          if (start_pos > len(cur_buf)) then
             found = .false.
             return
          end if

          brace_count = 1
          end_pos = start_pos + 1
          do while (end_pos <= len(cur_buf) .and. brace_count > 0)
             if (cur_buf(end_pos:end_pos) == '{' .or. cur_buf(end_pos:end_pos) == '[') brace_count = brace_count + 1
             if (cur_buf(end_pos:end_pos) == '}' .or. cur_buf(end_pos:end_pos) == ']') brace_count = brace_count - 1
             end_pos = end_pos + 1
          end do
          cur_buf = cur_buf(start_pos:end_pos-1)
       else
          ! Dernier élément du chemin : extraire la valeur
          start_pos = pos_colon + 1
          do while (start_pos <= len(cur_buf) .and. (cur_buf(start_pos:start_pos) == ' ' .or. cur_buf(start_pos:start_pos) == char(9)))
             start_pos = start_pos + 1
          end do
          if (start_pos > len(cur_buf)) then
             found = .false.
             return
          end if

          if (cur_buf(start_pos:start_pos) == '"') then
             ! Chaîne entre guillemets
             end_pos = start_pos + 1
             do while (end_pos <= len(cur_buf) .and. cur_buf(end_pos:end_pos) /= '"')
                end_pos = end_pos + 1
             end do
             val_str = cur_buf(start_pos+1:end_pos-1)
             found = .true.
             return
          else
             ! Nombre, booléen, etc.
             end_pos = start_pos
             do while (end_pos <= len(cur_buf) .and. &
                       cur_buf(end_pos:end_pos) /= ',' .and. &
                       cur_buf(end_pos:end_pos) /= '}' .and. &
                       cur_buf(end_pos:end_pos) /= ']' .and. &
                       cur_buf(end_pos:end_pos) /= ' ' .and. &
                       cur_buf(end_pos:end_pos) /= char(9) .and. &
                       cur_buf(end_pos:end_pos) /= char(10) .and. &
                       cur_buf(end_pos:end_pos) /= char(13))
                end_pos = end_pos + 1
             end do
             val_str = trim(adjustl(cur_buf(start_pos:end_pos-1)))
             found = .true.
             return
          end if
       end if
    end do
  end subroutine json_find_value_str

  subroutine json_get_string(buffer, path, val, default_val)
    character(len=*), intent(in) :: buffer, path, default_val
    character(len=*), intent(out) :: val
    character(len=:), allocatable :: s
    logical :: found

    call json_find_value_str(buffer, path, s, found)
    if (found .and. allocated(s)) then
       val = trim(s)
    else
       val = trim(default_val)
    end if
  end subroutine json_get_string

  subroutine json_get_real(buffer, path, val, default_val)
    character(len=*), intent(in) :: buffer, path
    real(dp), intent(out) :: val
    real(dp), intent(in) :: default_val
    character(len=:), allocatable :: s
    logical :: found
    integer :: io

    call json_find_value_str(buffer, path, s, found)
    if (found .and. allocated(s)) then
       read(s, *, iostat=io) val
       if (io /= 0) val = default_val
    else
       val = default_val
    end if
  end subroutine json_get_real

  subroutine json_get_int(buffer, path, val, default_val)
    character(len=*), intent(in) :: buffer, path
    integer, intent(out) :: val
    integer, intent(in) :: default_val
    character(len=:), allocatable :: s
    logical :: found
    integer :: io

    call json_find_value_str(buffer, path, s, found)
    if (found .and. allocated(s)) then
       read(s, *, iostat=io) val
       if (io /= 0) val = default_val
    else
       val = default_val
    end if
  end subroutine json_get_int

  subroutine json_get_bool(buffer, path, val, default_val)
    character(len=*), intent(in) :: buffer, path
    logical, intent(out) :: val
    logical, intent(in) :: default_val
    character(len=:), allocatable :: s
    logical :: found

    call json_find_value_str(buffer, path, s, found)
    if (found .and. allocated(s)) then
       if (index(s, "true") > 0 .or. index(s, "TRUE") > 0 .or. index(s, "1") > 0) then
          val = .true.
       else
          val = .false.
       end if
    else
       val = default_val
    end if
  end subroutine json_get_bool

  ! Détection et extraction des obstacles (BOX ou CYLINDER)
  subroutine json_get_obstacle(buffer, obs_type, x0, y0, x1, y1, xc, yc, radius, perm, found)
    character(len=*), intent(in) :: buffer
    integer, intent(out) :: obs_type
    real(dp), intent(out) :: x0, y0, x1, y1, xc, yc, radius, perm
    logical, intent(out) :: found

    integer :: pos_obs, pos_end
    character(len=:), allocatable :: obs_str
    character(len=32) :: type_str

    obs_type = OBS_NONE
    x0 = 0.0_dp; y0 = 0.0_dp; x1 = 0.0_dp; y1 = 0.0_dp
    xc = 0.0_dp; yc = 0.0_dp; radius = 0.0_dp; perm = 1.0e-40_dp
    found = .false.

    pos_obs = index(buffer, '"obstacles"')
    if (pos_obs <= 0) return

    pos_end = pos_obs + 11
    do while (pos_end <= len(buffer) .and. buffer(pos_end:pos_end) /= ']')
       pos_end = pos_end + 1
    end do
    if (pos_end > len(buffer)) return

    obs_str = buffer(pos_obs:pos_end)
    call json_get_string(obs_str, "type", type_str, "")

    if (trim(type_str) == "BOX") then
       obs_type = OBS_BOX
       call json_get_real(obs_str, "x0", x0, 0.0_dp)
       call json_get_real(obs_str, "y0", y0, 0.0_dp)
       call json_get_real(obs_str, "x1", x1, 0.0_dp)
       call json_get_real(obs_str, "y1", y1, 0.0_dp)
       call json_get_real(obs_str, "permeability", perm, 1.0e-40_dp)
       found = .true.
    else if (trim(type_str) == "CYLINDER") then
       obs_type = OBS_CYLINDER
       call json_get_real(obs_str, "xc", xc, 0.0_dp)
       call json_get_real(obs_str, "yc", yc, 0.0_dp)
       call json_get_real(obs_str, "radius", radius, 0.0_dp)
       call json_get_real(obs_str, "permeability", perm, 1.0e-40_dp)
       found = .true.
    end if
  end subroutine json_get_obstacle

end module module_json
