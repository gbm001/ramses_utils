program dump_level_cells_program
    use amr_utils
    implicit none
    
    character (LEN=:), allocatable :: input_filename
    character (LEN=:), allocatable :: temp_str
    character (LEN=256)            :: out_filename
    integer                        :: cmd_args
    integer                        :: length
    integer                        :: select_level, cells_size
    integer                        :: select_ivar_min, select_ivar_max
    double precision, allocatable  :: dump_cells(:,:,:,:)
    
    ! Validate number of command line arguments
    cmd_args = command_argument_count()
    if ((cmd_args /= 2) .AND. (cmd_args /= 4)) then
        write (6,*) "Wrong number of command line arguments (expecting 2 or 4)!"
        write (6,*) "First argument: filename of combined output file"
        write (6,*) "Second argument: desired level"
        write (6,*) "Third/fourth arguments: min/max hydro column of data"
        stop 1
    end if
    
    ! Read command line arguments
    call get_command_argument(1, length=length)
    allocate(character(LEN=length) :: input_filename)
    call get_command_argument(1, input_filename)
    call get_command_argument(2, length=length)
    allocate(character(LEN=length) :: temp_str)
    call get_command_argument(2, temp_str)
    read(temp_str, *) select_level
    
    if (select_level < 1) then
        stop "Must use a level >= 1!"
    end if
    
    if (cmd_args == 4) then
        call get_command_argument(3, length=length)
        deallocate(temp_str)
        allocate(character(LEN=length) :: temp_str)
        call get_command_argument(3, temp_str)
        read(temp_str, *) select_ivar_min
        call get_command_argument(4, length=length)
        deallocate(temp_str)
        allocate(character(LEN=length) :: temp_str)
        call get_command_argument(4, temp_str)
        read(temp_str, *) select_ivar_max
    end if
    
    ! Read single combined output file
    call load_single(input_filename)
    
    ! Select hydro columns and validate
    if (cmd_args == 2) then
       select_ivar_min = min_ivar
       select_ivar_max = max_ivar
    else
       if (select_ivar_min < min_ivar) then
           stop "Invalid selected hydro columns (lower limit)!"
       end if
       if (select_ivar_max > max_ivar) then
           stop "Invalid selected hydro columns (upper limit)"
       end if
    end if
    if (select_level > nlevelmax) then
        stop "Invalid level selected!"
    end if
    
    ! Allocate space
    cells_size = 2**select_level
    if (ndim==1) then
         allocate(dump_cells(1:1,1:1,1:cells_size,&
                             &select_ivar_min:select_ivar_max))
    else if (ndim==2) then
         allocate(dump_cells(1:1,1:cells_size,1:cells_size,&
                             &select_ivar_min:select_ivar_max))
    else if (ndim==3) then
         allocate(dump_cells(1:cells_size,1:cells_size,1:cells_size,&
                             &select_ivar_min:select_ivar_max))
    else
        stop 'ndim must be 1, 2 or 3?'
    end if
    
    ! Set a default in case there are gaps in this level (i.e. AMR)
    dump_cells = huge(0.0d0)
    
    ! Load data
    call dump_level_cells(select_level, select_ivar_min,&
                         &select_ivar_max, dump_cells, leaf_only=.TRUE.)
    
    write (out_filename,'(A,I0,A)') 'level_', select_level, '.dat'
    
    ! Write to level_N.dat file
    open(1, file=out_filename, access='stream', form='unformatted')
    write(1) ndim, cells_size, select_ivar_min, select_ivar_max
    write(1) dump_cells
    close(1)
    
    stop
    
end program dump_level_cells_program
