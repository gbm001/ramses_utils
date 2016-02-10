program power_spectrum_program
    use amr_utils
    use binning
    use fftw_module
    implicit none
    
!     double precision, parameter    :: pc=3.0856776d+18      ! pc in cm
    double precision, parameter    :: km_s=100000.d0        ! km/s in cm/s
    double precision               :: v_scale
    
    character (LEN=:), allocatable :: filename, helmholtz_filename
    character (LEN=:), allocatable :: temp_str
    character (LEN=5)              :: suffix
    character (LEN=2)              :: extra
    character (LEN=:), allocatable :: output_filename
    integer                        :: cmd_args, length
    integer                        :: select_level
    
    integer                        :: i, j, k, ilevel, ioutput, outputs
    integer                        :: N, Nbins
    integer                        :: max_level, max_grids
    integer                        :: gmin, gmax
    double precision, allocatable :: dump_cells(:,:,:,:)
    double precision, allocatable  :: vmag(:,:,:,:)
    complex(kind=cdp), allocatable :: vfft(:,:,:)
    double precision, allocatable  :: power1D(:,:)
    double precision, allocatable  :: power2D(:,:,:)
    double precision, allocatable  :: power3D(:,:,:,:)
    double precision, allocatable  :: bin_avg(:,:), bin_centres(:,:)
    double precision, allocatable  :: tr_bin_avg(:,:)
    integer, allocatable           :: counts(:,:)
    
    real (kind=dp), allocatable, target :: v(:,:,:)
    real (kind=dp), allocatable, target :: vcomp(:,:,:)
    real (kind=dp), allocatable, target :: vsol(:,:,:)
    
    lowmem_tables=.TRUE.
    keep_father=.FALSE.
    keep_next=.TRUE.
    keep_prev=.FALSE.
    keep_level=.FALSE.
    keep_flag1=.FALSE.
    keep_cpu_map=.FALSE.
    keep_nbor=.FALSE.
    keep_xg=.TRUE.
    keep_u=.TRUE.
    
    min_ivar=2
    max_ivar=4
    
    ! Validate number of command line arguments
    cmd_args = command_argument_count()
    if (cmd_args < 2 .OR. cmd_args > 3) then
        write (6,*) "Wrong number of command line arguments (expecting 2 or 3)!"
        write (6,*) "First argument: filename of combined output file"
        write (6,*) "Optional second argument: filename of Helmholtz decomposition file"
        write (6,*) "Final argument: desired level"
        write (6,*) "If using Helmholtz file, all velocities in centre-of-velocity frame"
        stop 1
    end if
    
    if (cmd_args == 2) outputs = 1
    if (cmd_args == 3) outputs = 3
    
    ! Read command line arguments
    call get_command_argument(1, length=length)
    allocate(character(LEN=length) :: filename)
    call get_command_argument(1, filename)
    if (cmd_args == 2) then
        call get_command_argument(2, length=length)
        allocate(character(LEN=length) :: temp_str)
        call get_command_argument(2, temp_str)
        read(temp_str, *) select_level
    else
        call get_command_argument(2, length=length)
        allocate(character(LEN=length) :: helmholtz_filename)
        call get_command_argument(2, helmholtz_filename)
        call get_command_argument(3, length=length)
        allocate(character(LEN=length) :: temp_str)
        call get_command_argument(3, temp_str)
        read(temp_str, *) select_level
    end if
    
    if (select_level < 1) then
        stop "Must use a level >= 1!"
    end if
    
    call load_single(filename)
    
    if (select_level > nlevelmax) then
        stop "Invalid level selected! (level > nlevelmax)"
    end if
    
    v_scale = (unit_l / unit_t) / km_s
    
    if (cmd_args == 3) then
        allocate(v(1:ngridmax, 1:8, 1:3))
        allocate(vcomp(1:ngridmax, 1:8, 1:3))
        allocate(vsol(1:ngridmax, 1:8, 1:3))
        open(10, file=helmholtz_filename, access='stream')
            read(10) v
            read(10) vcomp
            read(10) vsol
        close(10)
    end if

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
    N = 2**select_level
    gmin = -N / 2
    gmax = (N-1) / 2
    if (ndim==1) then
         allocate(dump_cells(1:1, 1:1, 1:N, 2:4))
         allocate(vmag(1:1, 1:1, 1:N, 1:outputs))
         allocate(vfft(1:1, 1:1, 1:N))
         allocate(power1D(gmin:gmax, 1:outputs))
    else if (ndim==2) then
         allocate(dump_cells(1:1, 1:N, 1:N, 2:4))
         allocate(vmag(1:1, 1:N, 1:N, 1:outputs))
         allocate(vfft(1:1, 1:N, 1:N))
         allocate(power2D(gmin:gmax, gmin:gmax, 1:outputs))
    else if (ndim==3) then
         allocate(dump_cells(1:N, 1:N, 1:N, 2:4))
         allocate(vmag(1:N, 1:N, 1:N, 1:outputs))
         allocate(vfft(1:N, 1:N, 1:N))
         allocate(power3D(gmin:gmax, gmin:gmax, gmin:gmax, 1:outputs))
    else
        stop 'ndim must be 1, 2 or 3?'
    end if
    
    ! Set a default in case there are gaps in this level (i.e. AMR)
    dump_cells = huge(0.0d0)
    
    ! Load data
    if (cmd_args == 2) then
        call dump_level_cells(select_level, u(:,:,2:4),&
                             &dump_cells, leaf_only=.TRUE.)
        vmag(:,:,:,1) = sqrt(&
           &sum(dump_cells(:, :, :, 2:4)**2, dim=ndim+1)) * v_scale
    else
        call dump_level_cells(select_level, v,&
                             &dump_cells, leaf_only=.TRUE.)
        vmag(:,:,:,1) = sqrt(sum(dump_cells(:, :, :, 2:4)**2, dim=ndim+1))
        
        call dump_level_cells(select_level, vcomp,&
                             &dump_cells, leaf_only=.TRUE.)
        vmag(:,:,:,2) = sqrt(sum(dump_cells(:, :, :, 2:4)**2, dim=ndim+1))
        
        call dump_level_cells(select_level, vsol,&
                             &dump_cells, leaf_only=.TRUE.)
        vmag(:,:,:,3) = sqrt(sum(dump_cells(:, :, :, 2:4)**2, dim=ndim+1))
    end if
    
    Nbins = ceiling(sqrt(dble(ndim * N**2)))
    allocate(counts(0:Nbins, 1:outputs))
    allocate(bin_centres(0:Nbins, 1:outputs))
    allocate(bin_avg(0:Nbins, 1:outputs))
    allocate(tr_bin_avg(0:Nbins, 1:outputs))
    do ioutput=1,outputs
        ! Note magic sign-flipping change in image domain
        ! so we don't need to shift later... (only works for even domains)
        ! http://dsp.stackexchange.com/questions/9039/centering-zero-frequency-for-discrete-fourier-transform
        if (ndim==1) then
            do k=1,N
                vmag(1,1,k,ioutput) = vmag(1,1,k,ioutput) * (-1)**(k-1)
            end do
            call DFT_1D(vmag(1,1,:,ioutput), vfft(1,1,:), N)
            power1D(:,ioutput) = reshape(dble(vfft * conjg(vfft)), (/N/))
        else if (ndim==2) then
            do k=1,N
                do j=1,N
                    vmag(1,j,k,ioutput) = vmag(1,j,k,ioutput) * (-1)**(j+k)
                end do
            end do
            call DFT_2D(vmag(1,:,:,ioutput), vfft(1,:,:), N)
            power2D(:,:,ioutput) = reshape(dble(vfft * conjg(vfft)), (/N,N/))
        else if (ndim==3) then
            do k=1,N
                do j=1,N
                    do i=1,N
                        vmag(i,j,k,ioutput) = vmag(i,j,k,ioutput) * (-1)**(i+j+k-1)
                    end do
                end do
            end do
            call DFT_3D(vmag(:,:,:,ioutput), vfft, N)
            power3D(:,:,:,ioutput) = dble(vfft * conjg(vfft))
            
        else
            stop 'ndim /= 1,2 or 3!'
        end if
        
        if (ndim==1) then
            call square_binning(power1D(:,ioutput), Nbins, counts(:,ioutput), bin_avg(:,ioutput), bin_centres(:,ioutput))
            call triangle_binning(power1D(:,ioutput), Nbins, tr_bin_avg(:,ioutput))
        else if (ndim==2) then
            call square_binning(power2D(:,:,ioutput), Nbins, counts(:,ioutput), bin_avg(:,ioutput), bin_centres(:,ioutput))
            call triangle_binning(power2D(:,:,ioutput), Nbins, tr_bin_avg(:,ioutput))
        else if (ndim==3) then
            call square_binning(power3D(:,:,:,ioutput), Nbins, counts(:,ioutput), bin_avg(:,ioutput), bin_centres(:,ioutput))
            call triangle_binning(power3D(:,:,:,ioutput), Nbins, tr_bin_avg(:,ioutput))
        else
            stop 'ndim /= 1,2 or 3!'
        end if
    end do
    
    suffix = get_single_suffix(filename)
    
    extra = ""
    if (outputs == 3) extra = '_h'
    if (suffix == "") then
        output_filename = "power_spec"//trim(extra)//".dat"
    else
        output_filename = "power_spec"//trim(extra)//"_"//suffix//".dat"
    end if
    
    open (10, file=output_filename)
    
    write (10, '(A)') '# (square binning) bin centre; data bin centre; &
                     &data bin average; (triangle binning) bin centre; &
                     &data bin average'
    
    do ioutput=1,outputs
        if (ioutput /= 1) write (10,*)
        do i=1,N/2
            write (10,'(5(G0,:,1X))') dble(i)+0.5d0, bin_centres(i,ioutput), bin_avg(i,ioutput), i, tr_bin_avg(i,ioutput)
        end do
    end do
    
    close(10)
    
    stop

end program power_spectrum_program