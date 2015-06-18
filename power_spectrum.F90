module fft
   implicit none
#include "fftw3.f"

end module fft

program power_spectrum_program
    use amr_utils
    use binning
    implicit none
    
!     double precision, parameter    :: pc=3.0856776d+18      ! pc in cm
    double precision, parameter    :: km_s=100000.d0        ! km/s in cm/s
    double precision               :: v_scale
    
    character (LEN=:), allocatable :: filename
    character (LEN=:), allocatable :: temp_str
    character (LEN=5)              :: suffix
    character (LEN=:), allocatable :: output_filename
    integer                        :: cmd_args, length
    integer                        :: select_level
    
    integer                        :: i, j, k, ilevel
    integer                        :: N, Nbins, shift
    integer                        :: max_level, max_grids
    integer                        :: gmin, gmax
    double precision, allocatable  :: dump_cells(:,:,:,:)
    double precision, allocatable  :: vmag(:,:,:)
    complex(kind=cdp), allocatable :: vfft(:,:,:)
    double precision, allocatable  :: power1D(:), power2D(:,:), power3D(:,:,:)
    double precision, allocatable  :: bin_avg(:), bin_centres(:), tr_bin_avg(:)
    integer, allocatable           :: counts(:)
    
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
    N = 2**select_level
    gmin = -N / 2
    gmax = (N-1) / 2
    if (ndim==1) then
         allocate(dump_cells(1:1, 1:1, 1:N, 2:4))
         allocate(vmag(1:1, 1:1, 1:N))
         allocate(vfft(1:1, 1:1, 1:N))
         allocate(power1D(gmin:gmax))
    else if (ndim==2) then
         allocate(dump_cells(1:1, 1:N, 1:N, 2:4))
         allocate(vmag(1:1, 1:N, 1:N))
         allocate(vfft(1:1, 1:N, 1:N))
         allocate(power2D(gmin:gmax, gmin:gmax))
    else if (ndim==3) then
         allocate(dump_cells(1:N, 1:N, 1:N, 2:4))
         allocate(vmag(1:N, 1:N, 1:N))
         allocate(vfft(1:N, 1:N, 1:N))
         allocate(power3D(gmin:gmax, gmin:gmax, gmin:gmax))
    else
        stop 'ndim must be 1, 2 or 3?'
    end if
    
    ! Set a default in case there are gaps in this level (i.e. AMR)
    dump_cells = huge(0.0d0)
    
    ! Load data
    call dump_level_cells(select_level, 2,&
                         &4, dump_cells, leaf_only=.TRUE.)
    
    vmag = sqrt(sum(dump_cells(:, :, :, 2:4)**2, dim=ndim+1))
    
    shift = N / 2
    
    ! Note magic input-array change so we don't need to shift later...
    if (ndim==1) then
        do k=1,N
            vmag(1,1,k) = vmag(1,1,k) * (-1)**(k-1)
        end do
        call FFT_1D(vmag(1,1,:), vfft(1,1,:), N)
        power1D = reshape(dble(vfft * conjg(vfft)), (/N/))
!         power1D = cshift(power1D, -shift)
    else if (ndim==2) then
        do k=1,N
            do j=1,N
                vmag(1,j,k) = vmag(1,j,k) * (-1)**(j+k)
            end do
        end do
        call FFT_2D(vmag(1,:,:), vfft(1,:,:), N)
        power2D = reshape(dble(vfft * conjg(vfft)), (/N,N/))
!         power2D = cshift(power2D, -shift, 1)
!         power2D = cshift(power2D, -shift, 2)
    else if (ndim==3) then
        do k=1,N
            do j=1,N
                do i=1,N
                    vmag(i,j,k) = vmag(i,j,k) * (-1)**(i+j+k-1)
                end do
            end do
        end do
        call FFT_3D(vmag, vfft, N)
        power3D = dble(vfft * conjg(vfft))
!         power3D = cshift(power3D, -shift, 1)
!         power3D = cshift(power3D, -shift, 2)
!         power3D = cshift(power3D, -shift, 3)
        
    else
        stop 'ndim /= 1,2 or 3!'
    end if
    
    Nbins = ceiling(sqrt(dble(ndim * N**2)))
    allocate(counts(0:Nbins))
    allocate(bin_centres(0:Nbins))
    allocate(bin_avg(0:Nbins))
    allocate(tr_bin_avg(0:Nbins))
    
    if (ndim==1) then
        call square_binning(power1D, Nbins, counts, bin_avg, bin_centres)
        call triangle_binning(power1D, Nbins, tr_bin_avg)
    else if (ndim==2) then
        call square_binning(power2D, Nbins, counts, bin_avg, bin_centres)
        call triangle_binning(power2D, Nbins, tr_bin_avg)
    else if (ndim==3) then
        call square_binning(power3D, Nbins, counts, bin_avg, bin_centres)
        call triangle_binning(power3D, Nbins, tr_bin_avg)
    else
        stop 'ndim /= 1,2 or 3!'
    end if
    
    suffix = get_single_suffix(filename)
    
    if (suffix == "") then
        output_filename = "power_spec.dat"
    else
        output_filename = "power_spec_"//suffix//".dat"
    end if
    
    open (10, file=output_filename)
    
    do i=1,N/2
        write (10,'(5(G0,:,1X))') dble(i)+0.d5, bin_centres(i), bin_avg(i), i, tr_bin_avg(i)
    end do
    
    close(10)
    
    stop
    
    contains
    
    subroutine FFT_1D(real_field, complex_field, N)
       use fft
       implicit none
       ! Transform purely real field into complex field

       integer, intent(in)            :: N     ! Size of grid
       real(kind=dp), intent(in)      :: real_field(0:N-1)
                                               ! Scalar field
       complex(kind=cdp), intent(out) :: complex_field(0:N-1)
                                               ! Complex field output
       
       integer (kind=i8b)             :: plan  ! FFTW plan pointer

       call dfftw_plan_dft_1d(plan, N, complex_field, &
          & complex_field, FFTW_FORWARD, FFTW_ESTIMATE)

       complex_field = cmplx(real_field, kind=cdp)

       !write (6,*) "Performing FFT..."

       call dfftw_execute_dft(plan, complex_field, complex_field)
       complex_field = complex_field / real(N, dp)
       
       call dfftw_destroy_plan(plan)
       
       return
    end subroutine FFT_1D

    subroutine FFT_2D(real_field, complex_field, N)
       use fft
       implicit none
       ! Transform purely real field into complex field

       integer, intent(in)            :: N     ! Size of grid
       real(kind=dp), intent(in)      :: real_field(0:N-1, 0:N-1)
                                               ! Scalar field
       complex(kind=cdp), intent(out) :: complex_field(0:N-1, 0:N-1)
                                               ! Complex field output
       
       integer (kind=i8b)             :: plan  ! FFTW plan pointer

       call dfftw_plan_dft_2d(plan, N, N, complex_field, &
          & complex_field, FFTW_FORWARD, FFTW_ESTIMATE)

       complex_field = cmplx(real_field, kind=cdp)

       !write (6,*) "Performing FFT..."
       call dfftw_execute_dft(plan, complex_field, complex_field)
       complex_field = complex_field / real(N**2, dp)
   
       call dfftw_destroy_plan(plan)

       return
    end subroutine FFT_2D

    subroutine FFT_3D(real_field, complex_field, N)
       use fft
       implicit none
       ! Transform purely real field into complex field

       integer, intent(in)            :: N     ! Size of grid
       real(kind=dp), intent(in)      :: real_field(0:N-1, 0:N-1, 0:N-1)
                                               ! Scalar field
       complex(kind=cdp), intent(out) :: complex_field(0:N-1, 0:N-1, 0:N-1)
                                               ! Complex field output
       
       integer (kind=i8b)             :: plan  ! FFTW plan pointer
       
       call dfftw_plan_dft_3d(plan, N, N, N, complex_field, &
          & complex_field, FFTW_FORWARD, FFTW_ESTIMATE)

       complex_field = cmplx(real_field, kind=cdp)

       !write (6,*) "Performing FFT..."
       call dfftw_execute_dft(plan, complex_field, complex_field)
       complex_field = complex_field / real(N**3, dp)
   
       call dfftw_destroy_plan(plan)

       return
    end subroutine FFT_3D

end program power_spectrum_program