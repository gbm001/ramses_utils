program reduce_nlevelmax_program
    use amr_utils

    character (LEN=:), allocatable :: output_basedir
    character (LEN=:), allocatable :: new_output_basedir
    character (LEN=:), allocatable :: integer_str
    integer                        :: ioutput
    integer                        :: length

    integer                        :: new_nlevelmax

    integer                        :: icpu
    character (LEN=12)             :: output_path_aux
    character (LEN=18)             :: amr_filename_aux
    character (LEN=20)             :: hydro_filename_aux
    character (LEN=:), allocatable :: output_path
    
    integer, allocatable           :: temp_int2D(:, :)

    if (command_argument_count() /= 4) then
        write (6,*) "Wrong number of command line arguments (expecting 4)!"
        write (6,*) "First argument: directory containing RAMSES outputs"
        write (6,*) "Second argument: output number of desired output"
        write (6,*) "Third argument: directory to put new RAMSES output"
        write (6,*) "Fourth argument: new nlevelmax"
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
    read(integer_str, *) new_nlevelmax
    deallocate(integer_str)

    call output_dir(ioutput, output_path_aux)
    if (len(output_basedir) == 0) then
        output_path = output_path_aux // '/'
    else
        output_path = output_basedir // '/' // output_path_aux // '/'
    end if

    call read_info_header(output_path, ioutput)

    do icpu=1,ncpu
        call amr_filename(ioutput, icpu, amr_filename_aux)
        write (6,*) "Reading file ", output_path//amr_filename_aux
        call read_amr(output_path//amr_filename_aux, icpu)
        
        call hydro_filename(ioutput, icpu, hydro_filename_aux)
        write (6,*) "Reading file ", output_path//hydro_filename_aux
        call read_hydro(output_path//hydro_filename_aux, icpu)
        
        allocate(temp_int2D(1:ncpu, 1:new_nlevelmax))
        temp_int2D = headl_cpu(:, 1:new_nlevelmax)
        call move_alloc(temp_int2D, headl_cpu)
        
        allocate(temp_int2D(1:ncpu, 1:new_nlevelmax))
        temp_int2D = taill_cpu(:, 1:new_nlevelmax)
        call move_alloc(temp_int2D, taill_cpu)
        
        allocate(temp_int2D(1:ncpu, 1:new_nlevelmax))
        temp_int2D = numbl_cpu(:, 1:new_nlevelmax)
        call move_alloc(temp_int2D, numbl_cpu)
        
        allocate(temp_int2D(1:10, 1:new_nlevelmax))
        temp_int2D = numbtot_cpu(:, 1:new_nlevelmax)
        call move_alloc(temp_int2D, numbtot_cpu)

        nlevelmax_cpu = new_nlevelmax

        write (6,*) "Writing file ", new_output_basedir//'/'//amr_filename_aux
        call write_amr(new_output_basedir//'/'//amr_filename_aux, icpu)
        write (6,*) "Writing file ", new_output_basedir//'/'//hydro_filename_aux
        call write_hydro(new_output_basedir//'/'//hydro_filename_aux, icpu)
        call deallocate_amr()
        call deallocate_hydro()
    end do


end program reduce_nlevelmax_program
