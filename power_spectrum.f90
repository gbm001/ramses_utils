program power_spectrum_program
    use amr_utils
    implicit none
    
    double precision, parameter    :: pc=3.0856776d+18      ! pc in cm
    double precision, parameter    :: km_s=100000.d0        ! km/s in cm/s
    double precision               :: v_scale
    
    character (LEN=:), allocatable :: filename
    character (LEN=:), allocatable :: temp_str
    integer                        :: cmd_args, length
    integer                        :: select_level
    
    integer                        :: ilevel, iilevel, ind
    integer                        :: max_level, max_grids, sample_size
    integer                        :: cells_size
    double precision, allocatable  :: dump_cells(:,:,:,:)
   
    lowmem_tables=.TRUE.
    keep_father=.FALSE.
    keep_next=.TRUE.
    keep_prev=.FALSE.
    keep_level=.FALSE.
    keep_flag1=.FALSE.
    keep_cpu_map=.FALSE.
    keep_nbor=.FALSE.
    keep_xg=.FALSE.
    keep_u=.TRUE.
    
    min_ivar=2
    max_ivar=4
    
    ! Validate number of command line arguments
    cmd_args = command_argument_count()
    if (cmd_args /= 2) then
        write (6,*) "Wrong number of command line arguments (expecting 2)!"
        write (6,*) "First argument: filename of combined output file"
        write (6,*) "Second argument: desired level"
        stop 1
    end if
    
    ! Read command line arguments
    call get_command_argument(1, length=length)
    allocate(character(LEN=length) :: filename)
    call get_command_argument(1, filename)
    call get_command_argument(2, length=length)
    allocate(character(LEN=length) :: temp_str)
    call get_command_argument(2, temp_str)
    read(temp_str, *) select_level
    
    if (select_level < 1) then
        stop "Must use a level >= 1!"
    end if
    
    call load_single(filename)
    
    if (select_level > nlevelmax) then
        stop "Invalid level selected! (level > nlevelmax)"
    end if
    
    v_scale = (unit_l / unit_t) / km_s
    
    ! Find out how many levels are actually being used
    do ilevel=nlevelmax,1,-1
        if (numbl(ilevel) /= 0) exit
    end do
    max_level = ilevel
    max_grids = sum(numbl)
    
!     if (select_level > max_level) then
!         stop "Invalid level selected! (level > maximum level)"
!     end if
    
    ! Allocate space
    cells_size = 2**select_level
    if (ndim==1) then
         allocate(dump_cells(1:1, 1:1, 1:cells_size, 2:4))
    else if (ndim==2) then
         allocate(dump_cells(1:1, 1:cells_size, 1:cells_size, 2:4))
    else if (ndim==3) then
         allocate(dump_cells(1:cells_size, 1:cells_size, 1:cells_size, 2:4))
    else
        stop 'ndim must be 1, 2 or 3?'
    end if
    
    ! Set a default in case there are gaps in this level (i.e. AMR)
    dump_cells = huge(0.0d0)
    
    ! Load data
    call dump_level_cells(select_level, 2,&
                         &4, dump_cells, leaf_only=.FALSE.)
    
    
    stop

end program power_spectrum_program