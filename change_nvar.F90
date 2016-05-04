program reduce_nlevelmax_program
    use amr_utils

    character (LEN=:), allocatable :: output_basedir
    character (LEN=:), allocatable :: new_output_basedir
    character (LEN=:), allocatable :: integer_str
    integer                        :: ioutput
    integer                        :: length
    integer                        :: new_nvar
    integer                        :: old_nvar
    logical                        :: set_grav_form

    integer                        :: icpu
    character (LEN=12)             :: output_path_aux
    character (LEN=18)             :: amr_filename_aux
    character (LEN=20)             :: hydro_filename_aux
    character (LEN=19)             :: grav_filename_aux
    character (LEN=:), allocatable :: output_path
    
    real (kind=dp), allocatable    :: u_temp(:,:)

    if (command_argument_count() /= 4) then
        write (6,*) "Wrong number of command line arguments (expecting 3)!"
        write (6,*) "First argument: directory containing RAMSES outputs"
        write (6,*) "Second argument: output number of desired output"
        write (6,*) "Third argument: directory to put new RAMSES output"
        write (6,*) "Fourth argument: new nvar. If negative, new gravity format."
        stop 1
    end if

    call get_command_argument(1, length=length)
    allocate(character(LEN=length) :: output_basedir)
    call get_command_argument(1, output_basedir)
    call get_command_argument(2, length=length)
    allocate(character(LEN=length) :: integer_str)
    call get_command_argument(2, integer_str)
    read(integer_str, *) ioutput
    deallocate(integer_str)
    call get_command_argument(3, length=length)
    allocate(character(LEN=length) :: new_output_basedir)
    call get_command_argument(3, new_output_basedir)
    allocate(character(LEN=length) :: integer_str)
    call get_command_argument(4, integer_str)
    read(integer_str, *) new_nvar
    deallocate(integer_str)
    
    set_grav_form = (new_nvar < 0)
    if (set_grav_form) new_nvar = -new_nvar

    call output_dir(ioutput, output_path_aux)
    if (len(output_basedir) == 0) then
        output_path = output_path_aux // '/'
    else
        output_path = output_basedir // '/' // output_path_aux // '/'
    end if

    call read_info_header(output_path, ioutput)

    do icpu=1,ncpu
        
        if (icpu > 1) nvar = old_nvar
        !ivar_max_use = 0
        
        call amr_filename(ioutput, icpu, amr_filename_aux)
        write (6,*) "Reading file ", output_path//amr_filename_aux
        call read_amr(output_path//amr_filename_aux, icpu)
        
        call hydro_filename(ioutput, icpu, hydro_filename_aux)
        write (6,*) "Reading file ", output_path//hydro_filename_aux
        call read_hydro(output_path//hydro_filename_aux, icpu)
        
        call grav_filename(ioutput, icpu, grav_filename_aux)
        write (6,*) "Reading file ", output_path//grav_filename_aux
        call read_grav(output_path//grav_filename_aux, icpu)
        
        old_nvar = nvar
        nvar = new_nvar
        ivar_max_use = nvar
        if (set_grav_form) grav_form = 1
        
        allocate(u_temp(1:ncell, 1:ivar_max_use))
        u_temp(:,1:min(old_nvar, new_nvar)) = &
            & u_cpu(:,1:min(old_nvar, new_nvar))
        call move_alloc(u_temp, u_cpu)

        write (6,*) "Writing file ", new_output_basedir//'/'//amr_filename_aux
        call write_amr(new_output_basedir//'/'//amr_filename_aux, icpu)
        write (6,*) "Writing file ", new_output_basedir//'/'//hydro_filename_aux
        call write_hydro(new_output_basedir//'/'//hydro_filename_aux, icpu)
        write (6,*) "Writing file ", new_output_basedir//'/'//grav_filename_aux
        call write_grav(new_output_basedir//'/'//grav_filename_aux, icpu)
        call deallocate_amr()
        call deallocate_hydro()
        call deallocate_grav()
        
    end do


end program reduce_nlevelmax_program
