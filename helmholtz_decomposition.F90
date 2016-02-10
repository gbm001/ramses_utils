program helmholtz_cells
    use, intrinsic :: iso_c_binding
    use amr_utils
    use fftw_module
    
    logical, parameter               :: mass_weighted_averaging = .TRUE.
    integer, parameter               :: ckind=C_DOUBLE_COMPLEX
    double precision, parameter      :: km_s=100000.d0        ! km/s in cm/s
    double precision                 :: v_scale

    character (LEN=:), allocatable   :: filename
    character (LEN=:), allocatable   :: temp_str
    character (LEN=5)                :: suffix
    character (LEN=:), allocatable   :: output_filename
    integer                          :: cmd_args, length
    integer                          :: select_level
    
    integer                          :: i, j, k, ilevel
    integer                          :: ki, kj, kk
    double precision                 :: k_vec(1:3), k_mag
    complex(kind=ckind)              :: k_unit(1:3)

    double precision, allocatable    :: v_grid(:,:,:,:)
    integer, allocatable             :: grid_cell(:,:,:,:)

    complex(kind=ckind), allocatable :: v_fft(:,:,:,:)
    complex(kind=ckind), allocatable :: vhelm_fft(:,:,:,:)
    complex(kind=ckind), allocatable :: vcomp_cmplx(:,:,:,:)
    complex(kind=ckind), allocatable :: vsol_cmplx(:,:,:,:)
    
    real (kind=dp), allocatable      :: v(:,:,:)
    real (kind=dp), allocatable      :: vcomp(:,:,:)
    real (kind=dp), allocatable      :: vsol(:,:,:)
    
    real (kind=dp)                   :: v_sum(1:3)
    real (kind=dp)                   :: vcomp_sum(1:3)
    real (kind=dp)                   :: vsol_sum(1:3)
    real (kind=dp)                   :: rho_sum
    
    integer                          :: N, grid_size
    integer                          :: gx, gy, gz
    integer                          :: igrid, ind, igrid_son, ind_son
    double precision                 :: xg_loc(1:3)
    
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

    min_ivar=1
    max_ivar=4

    ! Validate number of command line arguments

    cmd_args = command_argument_count()
    if (cmd_args /= 2) then
        write (6,*) "Wrong number of command line arguments (expecting 2)!"
        write (6,*) "First argument: filename of combined output file"
        write (6,*) "Second argument: maximum desired level to decompose"
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
    
    if (ndim /= 3) then
        stop 'ndim /= 3!'
    end if

    if (select_level > nlevelmax) then
        stop "Invalid level selected! (level > nlevelmax)"
    end if

    v_scale = (unit_l / unit_t) / km_s
    
    allocate(v(1:ngridmax, 1:twotondim, 1:3))
    allocate(vcomp(1:ngridmax, 1:twotondim, 1:3))
    allocate(vsol(1:ngridmax, 1:twotondim, 1:3))
    
    grid_size = 2**(select_level-1)
    N = 2**select_level

    allocate(v_grid(1:N, 1:N, 1:N, 1:3))
    allocate(grid_cell(1:N, 1:N, 1:N, 1:2))
    v_grid = huge(v_grid)
    grid_cell = -1
    
    allocate(v_fft(1:N, 1:N, 1:N, 1:3))
    allocate(vhelm_fft(1:N, 1:N, 1:N, 1:3))
    allocate(vcomp_cmplx(1:N, 1:N, 1:N, 1:3))
    allocate(vsol_cmplx(1:N, 1:N, 1:N, 1:3))

    ! Put the velocities in a grid, keep track of where they came from
    igrid = headl(select_level)
    do
        xg_loc(1:ndim) = xg(igrid, 1:ndim)
        gx = nint(xg_loc(1) * grid_size + 0.5)
        gy = nint(xg_loc(2) * grid_size + 0.5)
        gz = nint(xg_loc(3) * grid_size + 0.5)

        v_grid(2*gx-1, 2*gy-1, 2*gz-1, :) = u(igrid, 1, 2:4)
        v_grid(2*gx  , 2*gy-1, 2*gz-1, :) = u(igrid, 2, 2:4)
        v_grid(2*gx-1, 2*gy  , 2*gz-1, :) = u(igrid, 3, 2:4)
        v_grid(2*gx  , 2*gy  , 2*gz-1, :) = u(igrid, 4, 2:4)
        v_grid(2*gx-1, 2*gy-1, 2*gz  , :) = u(igrid, 5, 2:4)
        v_grid(2*gx  , 2*gy-1, 2*gz  , :) = u(igrid, 6, 2:4)
        v_grid(2*gx-1, 2*gy  , 2*gz  , :) = u(igrid, 7, 2:4)
        v_grid(2*gx  , 2*gy  , 2*gz  , :) = u(igrid, 8, 2:4)

        grid_cell(2*gx-1, 2*gy-1, 2*gz-1, :) = (/igrid, 1/)
        grid_cell(2*gx  , 2*gy-1, 2*gz-1, :) = (/igrid, 2/)
        grid_cell(2*gx-1, 2*gy  , 2*gz-1, :) = (/igrid, 3/)
        grid_cell(2*gx  , 2*gy  , 2*gz-1, :) = (/igrid, 4/)
        grid_cell(2*gx-1, 2*gy-1, 2*gz  , :) = (/igrid, 5/)
        grid_cell(2*gx  , 2*gy-1, 2*gz  , :) = (/igrid, 6/)
        grid_cell(2*gx-1, 2*gy  , 2*gz  , :) = (/igrid, 7/)
        grid_cell(2*gx  , 2*gy  , 2*gz  , :) = (/igrid, 8/)

        if (igrid == taill(select_level)) exit
        igrid = next(igrid)
    end do
    
    if (any(grid_cell(:,:,:,1) == -1)) then
        stop 'Selected level is not fully populated; this will fail'
    end if
    
    ! Subtract average velocity and scale
    v_grid(:,:,:,1) = v_grid(:,:,:,1) - (sum(v_grid(:,:,:,1)) / N**3)
    v_grid(:,:,:,2) = v_grid(:,:,:,2) - (sum(v_grid(:,:,:,2)) / N**3)
    v_grid(:,:,:,3) = v_grid(:,:,:,3) - (sum(v_grid(:,:,:,3)) / N**3)
    v_grid = v_grid * v_scale

    ! Helmholtz decomposition
    call DFT_3D(v_grid(:,:,:,1), v_fft(:,:,:,1), N)
    call DFT_3D(v_grid(:,:,:,2), v_fft(:,:,:,2), N)
    call DFT_3D(v_grid(:,:,:,3), v_fft(:,:,:,3), N)

    ! Find the purely compressive field (find the longitudinal field)
    do k=1,N
        kk = k-1
        if (kk > N/2) kk = kk-N
        do j=1,N
            kj = j-1
            if (kj > N/2) kj = kj-N
            do i=1,N
                ki = i-1
                if (ki > N/2) ki = ki-N
                k_vec = dble((/ki, kj, kk/))
                k_mag = max(1.d0, sqrt(sum(k_vec**2)))
                k_unit = dcmplx(k_vec / k_mag)
                vhelm_fft(i,j,k,:) = &
                    & dot_product(k_unit, v_fft(i,j,k,:)) * k_unit
            end do
        end do
    end do
    
    vhelm_fft(1,1,1,:) = 0.d0  ! We don't care about the constant part

    call IFT_3D(vhelm_fft(:,:,:,1), vcomp_cmplx(:,:,:,1), N)
    call IFT_3D(vhelm_fft(:,:,:,2), vcomp_cmplx(:,:,:,2), N)
    call IFT_3D(vhelm_fft(:,:,:,3), vcomp_cmplx(:,:,:,3), N)

    ! Find the purely solenoidal field (find the transverse field)
    vhelm_fft = v_fft - vhelm_fft
    vhelm_fft(1,1,1,:) = 0.d0  ! We don't care about the constant part
    call IFT_3D(vhelm_fft(:,:,:,1), vsol_cmplx(:,:,:,1), N)
    call IFT_3D(vhelm_fft(:,:,:,2), vsol_cmplx(:,:,:,2), N)
    call IFT_3D(vhelm_fft(:,:,:,3), vsol_cmplx(:,:,:,3), N)

    ! Store v, vcomp, vsol in cell grids
    do k=1,N
        do j=1,N
            do i=1,N
                igrid = grid_cell(i,j,k,1)
                ind = grid_cell(i,j,k,2)
                v(igrid, ind, 1:3) = v_grid(i,j,k,1:3)
                vcomp(igrid, ind, 1:3) = dble(vcomp_cmplx(i,j,k,1:3))
                vsol(igrid, ind, 1:3) = dble(vsol_cmplx(i,j,k,1:3))
            end do
        end do
    end do
    
    deallocate(v_grid)
    deallocate(grid_cell)
    deallocate(v_fft, vhelm_fft, vcomp_cmplx, vsol_cmplx)

    ! Now loop over levels, averaging upwards
    do ilevel=select_level-1, 1, -1
        igrid = headl(ilevel)
        do
            do ind=1,8
                igrid_son = son(igrid, ind)
                v_sum = 0.d0; vcomp_sum = 0.d0; vsol_sum = 0.d0
                if (mass_weighted_averaging) rho_sum = 0.d0
                do ind_son=1,8
                    if (mass_weighted_averaging) then
                        rho_sum = rho_sum + u(igrid_son, ind_son, 1)
                        v_sum = v_sum + u(igrid_son, ind_son, 1) * v(igrid_son, ind_son, 1:3)
                        vcomp_sum = vcomp_sum + u(igrid_son, ind_son, 1) * vcomp(igrid_son, ind_son, 1:3)
                        vsol_sum = vsol_sum + u(igrid_son, ind_son, 1) * vsol(igrid_son, ind_son, 1:3)
                    else
                        v_sum = v_sum + v(igrid_son, ind_son, 1:3)
                        vcomp_sum = vcomp_sum + vcomp(igrid_son, ind_son, 1:3)
                        vsol_sum = vsol_sum + vsol(igrid_son, ind_son, 1:3)
                    end if
                end do
                if (mass_weighted_averaging) then
                    v(igrid, ind, 1:3) = v_sum / rho_sum
                    vcomp(igrid, ind, 1:3) = vcomp_sum / rho_sum
                    vsol(igrid, ind, 1:3) = vsol_sum / rho_sum
                else
                    v(igrid, ind, 1:3) = v_sum / 8.d0
                    vcomp(igrid, ind, 1:3) = vcomp_sum / 8.d0
                    vsol(igrid, ind, 1:3) = vsol_sum / 8.d0
                end if
            end do

            if (igrid == taill(ilevel)) exit
            igrid = next(igrid)
        end do
    end do
    
    ! Now write to a file
    
    suffix = get_single_suffix(filename)

    if (suffix == "") then
        output_filename = "helmholtz_velocities.dat"
    else
        output_filename = "helmholtz_velocities_"//suffix//".dat"
    end if

    open (10, file=output_filename, access='STREAM')
    
    write(10) v
    write(10) vcomp
    write(10) vsol
    
    close(10)

end program helmholtz_cells