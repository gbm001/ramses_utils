program rms_vel_3d_program
    use amr_utils
    
    character (LEN=:), allocatable :: output_basedir
    character (LEN=:), allocatable :: ioutput_str
    integer                        :: ioutput
    integer                        :: length
    
    integer                        :: ilevel, ind, i
    integer                        :: igrid, ngrid
    
    double precision               :: vsqd_sum
    
    ! We only need 'next' and 'u'
    keep_father=.FALSE.
    keep_next=.TRUE.
    keep_prev=.FALSE.
    keep_level=.FALSE.
    keep_flag1=.FALSE.
    keep_cpu_map=.FALSE.
    keep_nbor=.FALSE.
    keep_xg=.FALSE.
    keep_u=.TRUE.
    ! Only store u,v,w
    min_ivar = 2
    max_ivar = 4
    
    if (command_argument_count() /= 2) then
        write (6,*) "Wrong number of command line arguments (expecting 2)!"
        write (6,*) "First argument: directory containing RAMSES outputs"
        write (6,*) "Second argument: output number of desired output"
        stop 1
    end if
    
    call get_command_argument(1, length=length)
    allocate(character(LEN=length) :: output_basedir)
    call get_command_argument(1, output_basedir)
    call get_command_argument(2, length=length)
    allocate(character(LEN=length) :: ioutput_str)
    call get_command_argument(2, ioutput_str)
    read(ioutput_str, *) ioutput
    
    call read_all(output_basedir, ioutput)
    
    if (.NOT. (keep_next .AND. keep_u)) then
        stop "Need 'next' and 'u' data!"
    end if
    
    vsqd_sum = 0.d0
    
    do ilevel=nlevelmax, 1, -1
        ngrid = numbl(ilevel)
        if (ngrid==0) cycle
        ! Otherwise we have found the lowest level with cells
        igrid = headl(ilevel)
        do i=1,ncache
            do ind=1,twotondim
                vsqd_sum = vsqd_sum + sum(u(igrid, ind, 2:4)**2)
            end do
            igrid = next(igrid)
        end do
        exit ! only do the lowest occupied level (assume fixed grid)
    end do
    
    vsqd_sum = vsqd_sum / dble(twotondim)
    vsqd_sum = vsqd_sum / dble(ngrid)
    
    write(6,*) "rms velocity: ", sqrt(vsqd_sum)
    
!     integer, intent(in)               :: n_cells
!     integer, intent(in)               :: ivar
!     double precision, intent(out)     :: cells(1:n_cells)
!     
!     integer                           :: ilevel, ind, i
!     integer                           :: igrid, ncache, slot
!     
!     if (.NOT. (keep_next .AND. keep_u)) then
!         stop "Cannot get_hydro without 'next' and 'u' data!"
!     end if
!     
!     slot = 1
!     do ilevel=1, nlevelmax
!         igrid = headl(ilevel)
!         ncache = numbl(ilevel)
!         do i=1,ncache
!             do ind=1,twotondim
!                 cells(slot) = u(igrid, ind, ivar)
!                 slot = slot + 1
!             end do
!             igrid = next(igrid)
!         end do
!     end do
    
end program rms_vel_3d_program