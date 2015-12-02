module amr_utils
    use, intrinsic :: ISO_FORTRAN_ENV, only: FILE_STORAGE_SIZE, &
                                             & INT64, REAL128
    implicit none
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Description of combined tree format
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! This module allows normal (periodic) RAMSES outputs to be combined into
    ! a single complete tree. This tree is stored slightly differently to
    ! that in RAMSES. Variables relating to RAMSES files are given their
    ! usual names but with the suffix _cpu. Variables relating the the combined
    ! tree are not given a suffix. It is possible to store only some of the
    ! available data in the combined trees.
    !
    ! Like the RAMSES trees, the combined tree has the usual son, father, next,
    ! prev, flag1, cpu_map, nbor, xg and (for hydro) u variables. An extra
    ! variable level is available, storing the level of that oct. The tree, as
    ! in RAMSES, consists of a collection of octs (in 3D). However, the RAMSES
    ! mechanism of addressing individual cells by ncoarse + ind*gridmax is not
    ! used. First there is no single coarse cell at the start (i.e. ncoarse=0);
    ! oct number 1 is the root oct of the tree. Secondly, the 8 octs are
    ! addressed by indexing the last index of the relevant variables.
    ! Variables:
    !
    !   note - [1] refers to the nth index of each variable e.g. nbor([1], [2])
    ! twondim = 2.0 * ndim;   twotondim = 2.0**ndim
    ! nlevelmax - at least as large as the maximum level in the tree
    ! ngridmax - at least as large as the maximum number of octs
    ! headl(1:nlevelmax) - start of the linked list of octs on level [1]
    ! taill(1:nlevelmax) - end of the linked list of octs on level [1]
    ! numbl(1:nlevelmax) - number of octs on level [1]
    ! son(1:ngridmax, 1:twotondim) - [2]nd child of oct [1] (0 for no child)
    ! father(1:ngridmax) - father oct of oct [1]
    ! next(1:ngridmax) - next oct in headl/taill linked list for oct [1]
    ! prev(1:ngridmax) - prev oct in headl/taill linked list for oct [1]
    ! level(1:ngridmax) - level of oct [1]
    ! flag1(1:ngridmax) - flag1 value of oct [1]
    ! cpu_map(1:ngridmax) - original cpu of oct [1]
    ! nbor(1:ngridmax, 1:twondim) - the [2]nd neighbour oct of oct [1]
    ! xg(1:ngridmax, 1:ndim) - the 3D grid centre of oct [1]
    ! u(1:ngridmax, 1:twotondim, ivar_min_use:ivar_max_use) - the hydro
    !              variables for oct [1], [2]nd child, columns [3]
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Settings to be modified
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Store translation tables in scratch files to save memory?
    logical :: lowmem_tables=.FALSE.
    ! In extremely brief testing on small files this didn't make a massive
    ! difference (< 4%)
    
    ! Which data to keep in the combined tree ('son' is already kept)
    logical :: keep_father=.TRUE.
    logical :: keep_next=.TRUE.
    logical :: keep_prev=.TRUE.
    logical :: keep_level=.TRUE.
    logical :: keep_flag1=.TRUE.
    logical :: keep_cpu_map=.TRUE.
    logical :: keep_nbor=.TRUE.
    logical :: keep_xg=.TRUE.
    logical :: keep_u=.TRUE.
    ! Further extremely brief testing showed that making these parameters
    ! made no discernable difference; YMMV. Left as variables for easier
    ! python integration.
    
    integer            :: min_ivar=0  ! Minimum ivar, change externally
    integer            :: max_ivar=0  ! Maximum ivar, change externally
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Data storage
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Assume double precision
    integer, parameter :: sp=kind(1.0E0)
    integer, parameter :: dp=kind(1.0D0)   ! real*8
    integer, parameter :: qdp=REAL128      ! real*16 ! QUADHILBERT ONLY
    integer, parameter :: i8b=INT64
    integer, parameter :: csp=kind((1.0_sp, 1.0_sp))
    integer, parameter :: cdp=kind((1.0_dp, 1.0_dp))
    integer, parameter :: maxout=1000
    integer, parameter :: maxlevel=100
    integer, parameter :: overload=1  ! assume this isn't changed
    
    ! File implementation details
    integer, parameter :: int_iosize = STORAGE_SIZE(1) / FILE_STORAGE_SIZE
    integer, parameter :: sp_iosize = STORAGE_SIZE(1.0_sp) / FILE_STORAGE_SIZE
    integer, parameter :: dp_iosize = STORAGE_SIZE(1.0_dp) / FILE_STORAGE_SIZE
    
    ! RAMSES header parameters
    integer :: ncpu, ndim, nx, ny, nz, nlevelmax_cpu, ngridmax_cpu
    integer :: nboundary, ngrid_current_cpu
    integer :: noutput, iout, ifout, nstep, nstep_coarse
    
    ! Values derived from RAMSES header
    integer :: ndomain, ncoarse, ncell, tot_cells_cpu
    integer :: twondim, twotondim, ncache
    
    ! Hydro variables
    integer :: nvar, ivar_min_use, ivar_max_use
    
    ! Basic RAMSES floating point parameters
    ! From amr files:
    real (kind=dp) :: boxlen, tout(1:maxout), aout(1:maxout), t
    real (kind=dp) :: dtold(1:maxlevel), dtnew(1:maxlevel)
    real (kind=dp) :: const, mass_tot_0, rho_tot
    real (kind=dp) :: omega_m, omega_l, omega_k, omega_b, h0, aexp_ini
    real (kind=dp) :: boxlen_ini, aexp, hexp, aexp_old
    real (kind=dp) :: epot_tot_int, epot_tot_old, mass_sph
    ! From hydro files:
    real (kind=dp) :: gamma
    ! From info header:
    real (kind=dp) :: unit_l=1.d0, unit_d=1.d0, unit_t=1.d0
    
    ! Level and boundary linked list info
    integer, allocatable :: headl_cpu(:,:), taill_cpu(:,:), numbl_cpu(:,:), numbtot_cpu(:,:)
    integer, allocatable :: headb(:,:), tailb(:,:), numbb(:,:)
    integer :: headf, tailf, numbf, used_mem, used_mem_tot
    
    ! Other RAMSES header information
    character (LEN=128) :: ordering
    logical :: use_qdp
    real (kind=dp), allocatable :: bound_key_dp(:)
    real (kind=qdp), allocatable :: bound_key_qdp(:)
    
    ! RAMSES AMR data
    integer, allocatable :: son_cpu(:), father_cpu(:), next_cpu(:), prev_cpu(:)
    integer, allocatable :: flag1_cpu(:), cpu_map_cpu(:), nbor_cpu(:,:)
    
    real (kind=dp), allocatable :: xg_cpu(:,:), u_cpu(:,:)
    
    ! Combined AMR tree data
    integer :: ngridmax, nlevelmax, tot_cells, free
    integer, allocatable :: headl(:), taill(:), numbl(:)
    integer, allocatable :: son(:,:), father(:)
    integer, allocatable :: next(:), prev(:)
    integer, allocatable :: level(:), flag1(:), cpu_map(:)
    integer, allocatable :: nbor(:,:)
    integer, allocatable :: translation_table(:)
    integer, allocatable :: translation_tables(:,:)
    
    real (kind=dp), allocatable :: xg(:,:)
    real (kind=dp), allocatable :: u(:,:,:)
    real (kind=dp), allocatable :: unscaled_u(:,:,:)
    
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! File helper routines
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine amr_filename(ioutput, icpu, filename)
        ! Construct RAMSES-style AMR filename
        integer, intent(in)                   :: ioutput, icpu
        character (LEN=18), intent(out)       :: filename
        
        write(filename, '(A,I5.5,A,I5.5)') 'amr_', ioutput, '.out', icpu
        
        return
    end subroutine amr_filename
    
    subroutine hydro_filename(ioutput, icpu, filename)
        ! Construct RAMSES-style AMR hydro filename
        integer, intent(in)                   :: ioutput, icpu
        character (LEN=20), intent(out)       :: filename
        
        write(filename, '(A,I5.5,A,I5.5)') 'hydro_', ioutput, '.out', icpu
        
        return
    end subroutine hydro_filename
    
    subroutine info_filename(ioutput, filename)
        ! Construct RAMSES-style info header filename
        integer, intent(in)                   :: ioutput
        character (LEN=14), intent(out)       :: filename
        
        write(filename, '(A,I5.5,A)') 'info_', ioutput, '.txt'
        
        return
    end subroutine info_filename
    
    subroutine output_dir(ioutput, output_dir_out)
        ! Construct RAMSES output directory name
        integer, intent(in)                :: ioutput
        character (LEN=12), intent(out)    :: output_dir_out
        
        write(output_dir_out, '(A,I5.5)') 'output_', ioutput
        
        return
    end subroutine output_dir
    
    subroutine read_info_header(output_path, ioutput)
        character (LEN=*), intent(in)      :: output_path
        integer, intent(in)                :: ioutput
        
        character (LEN=14)                 :: header_name_aux
        character (LEN=:), allocatable     :: header_file_path
        
        character (LEN=2000)               :: line
        character (LEN=:), allocatable     :: item, quantity_str
        
        call info_filename(ioutput, header_name_aux)
        header_file_path = output_path//'/'//header_name_aux
        
        ! Primitive parsing of some of the info file
        ncpu = 0
        write (6,'(A,A,A)') "'", header_file_path, "'"
        open(10, file=header_file_path, status='old', form='formatted')
        do
            line = ""
            read(10, fmt='(A)', end=200) line
            ! Check if this is a data entry
            if (line(13:13) /= '=') cycle
            ! What item?
            item = adjustl(trim(line(1:12)))
            quantity_str = adjustl(trim(line(14:)))
            select case (item)
                case ('ncpu')
                    read(quantity_str, *) ncpu
                case ('unit_l')
                    read(quantity_str, *) unit_l
                case ('unit_d')
                    read(quantity_str, *) unit_d
                case ('unit_t')
                    read(quantity_str, *) unit_t
            end select
        end do
        200 continue
        
        close(10)
        
        if (ncpu == 0) stop 'Error reading info file, ncpu not found!'
    
        return
    end subroutine read_info_header
    
    subroutine read_all(output_basedir, ioutput)
        character (LEN=*), intent(in)      :: output_basedir
        integer, intent(in)                :: ioutput
        
        integer                            :: icpu
        integer                            :: nlevelmax, ngrid_current
        character (LEN=12)                 :: output_path_aux
        character (LEN=18)                 :: amr_filename_aux
        character (LEN=20)                 :: hydro_filename_aux
        character (LEN=:), allocatable     :: output_path
        
        call output_dir(ioutput, output_path_aux)
        if (len(output_basedir) == 0) then
            output_path = output_path_aux // '/'
        else
            output_path = output_basedir // '/' // output_path_aux // '/'
        end if
        
        call read_info_header(output_path, ioutput)
        
        tot_cells = 0
        nlevelmax = 0
        ngrid_current = 0
        
        do icpu=1,ncpu
            call amr_filename(ioutput, icpu, amr_filename_aux)
            write (6,*) "Reading header from file ", &
                &output_path//amr_filename_aux
            call read_amr_header(output_path//amr_filename_aux, icpu)
            
            tot_cells = tot_cells + tot_cells_cpu
            ngrid_current = ngrid_current + ngrid_current_cpu
            nlevelmax = max(nlevelmax, nlevelmax_cpu)
            
            call hydro_filename(ioutput, icpu, hydro_filename_aux)
            call read_hydro_header(output_path//hydro_filename_aux, icpu)
            deallocate(headl_cpu, taill_cpu, numbl_cpu, numbtot_cpu)
            if (use_qdp) then
                deallocate(bound_key_qdp)
            else
                deallocate(bound_key_dp)
            end if
        end do
        
        call allocate_combined(nlevelmax, ngrid_current, tot_cells)
        
        do icpu=1,ncpu
            call amr_filename(ioutput, icpu, amr_filename_aux)
            write (6,*) "Reading file ", output_path//amr_filename_aux
            call read_amr(output_path//amr_filename_aux, icpu)
            call add_tree(icpu)
            call deallocate_amr()
        end do
        
        do icpu=1,ncpu
            call amr_filename(ioutput, icpu, amr_filename_aux)
            call hydro_filename(ioutput, icpu, hydro_filename_aux)
            write (6,*) "Reading file ", output_path//amr_filename_aux
            call read_amr(output_path//amr_filename_aux, icpu)
            write (6,*) "Reading file ", output_path//hydro_filename_aux
            if (keep_u) then
                call read_hydro(output_path//hydro_filename_aux, icpu)
            end if
            
            call add_neighbours_hydro(icpu)
            call deallocate_amr()
            call deallocate_hydro()
        end do
        
        call deallocate_translations()
        
    end subroutine read_all
    
    subroutine convert_to_single(output_basedir, ioutput)
        character (LEN=*), intent(in)      :: output_basedir
        integer, intent(in)                :: ioutput
        
        character (LEN=12)                 :: output_path_aux
        character (LEN=16)                 :: single_filename
        character (LEN=:), allocatable     :: output_path
        
        call read_all(output_basedir, ioutput)
        
        call output_dir(ioutput, output_path_aux)
        output_path = output_basedir // output_path_aux // '/'
        if (len(output_basedir) == 0) then
            output_path = output_path_aux // '/'
        else
            output_path = output_basedir // '/' // output_path_aux // '/'
        end if
        
        write (single_filename, '(A,I5.5,A)') 'output_', ioutput, '.dat'
        
        output_path = output_path // single_filename
        
        call save_single(output_path)
        write (6,*) "Created output file ", output_path
    end subroutine convert_to_single
    
    function get_single_suffix(filename)
        ! Take a conmbined tree single file filename and separate out the
        ! RAMSES-style suffix - e.g. for output_00029.dat, return '00029'
        character (LEN=*), intent(in)       :: filename
        character (LEN=5)                   :: get_single_suffix
        
        integer                             :: len_filename
        character (LEN=10)                  :: string_end
        
        get_single_suffix = ""
        len_filename = len(filename)
        
        if (len_filename < 11) return
        
        string_end = filename(len_filename-9:len_filename)
    
        if (string_end(1:1) /= '_') return
        if ((string_end(7:10) /= '.dat') .AND.&
            &(string_end(7:10) /= '.DAT')) return
        
        if (verify(string_end(2:6), '0123456789') /= 0) return
        
        get_single_suffix = string_end(2:6)
        
    end function get_single_suffix
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! File input/output routines (CPU trees)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine read_amr_header(filename, icpu)
        ! Read header from RAMSES AMR cpu file (calls read_amr_header_internal)
        character (LEN=*), intent(in)           :: filename
        integer, intent(in)                     :: icpu
        
        integer                                 :: iunit
        
        iunit = 10 + icpu
        
        open(unit=iunit, file=trim(filename), status='old', form='unformatted')
        call read_amr_header_internal(iunit)
        close(unit=iunit)
        
        tot_cells_cpu = sum(numbl_cpu(icpu, :))
        
        !deallocate(headl_cpu, taill_cpu, numbl_cpu, numbtot_cpu)
        
        return
    end subroutine read_amr_header
    
    subroutine read_amr_header_internal(iunit)
        ! Read header from RAMSES AMR cpu file
        integer, intent(in)                     :: iunit
        integer                                 :: ierr
    
        read(unit=iunit) ncpu
        read(unit=iunit) ndim
        read(unit=iunit) nx, ny, nz
        read(unit=iunit) nlevelmax_cpu
        read(unit=iunit) ngridmax_cpu
        read(unit=iunit) nboundary
        read(unit=iunit) ngrid_current_cpu
        read(unit=iunit) boxlen
        read(unit=iunit) noutput, iout, ifout
        read(unit=iunit) tout(1:noutput)
        read(unit=iunit) aout(1:noutput)
        read(unit=iunit) t
        read(unit=iunit) dtold(1:nlevelmax_cpu)
        read(unit=iunit) dtnew(1:nlevelmax_cpu)
        read(unit=iunit) nstep, nstep_coarse
        read(unit=iunit) const, mass_tot_0, rho_tot
        read(unit=iunit) omega_m, omega_l, omega_k, omega_b, h0,&
                      & aexp_ini, boxlen_ini
        read(unit=iunit) aexp, hexp, aexp_old, epot_tot_int, epot_tot_old
        read(unit=iunit) mass_sph
        
        twondim = 2 * ndim
        twotondim = 2**ndim
        ncoarse = nx * ny * nz
        ncell = ncoarse + twotondim*ngridmax_cpu
        ndomain = ncpu * overload
        
        allocate(headl_cpu(1:ncpu, 1:nlevelmax_cpu), taill_cpu(1:ncpu, 1:nlevelmax_cpu))
        allocate(numbl_cpu(1:ncpu, 1:nlevelmax_cpu), numbtot_cpu(1:10, 1:nlevelmax_cpu))
        read(unit=iunit) headl_cpu(1:ncpu, 1:nlevelmax_cpu)
        read(unit=iunit) taill_cpu(1:ncpu, 1:nlevelmax_cpu)
        read(unit=iunit) numbl_cpu(1:ncpu, 1:nlevelmax_cpu)
        read(unit=iunit) numbtot_cpu(1:10, 1:nlevelmax_cpu)
!         allocate(headb(1:nboundary, 1:nlevelmax_cpu))
!         allocate(tailb(1:nboundary, 1:nlevelmax_cpu))
!         allocate(numbb(1:nboundary, 1:nlevelmax_cpu))
!         read(unit=iunit) headb(1:nboundary, 1:nlevelmax_cpu)
!         read(unit=iunit) tailb(1:nboundary, 1:nlevelmax_cpu)
!         read(unit=iunit) numbb(1:nboundary, 1:nlevelmax_cpu)
        read(unit=iunit) headf, tailf, numbf, used_mem, used_mem_tot
        read(unit=iunit) ordering
        if (trim(ordering) /= 'hilbert') stop 'Need hilbert ordering!'
        allocate(bound_key_dp(0:ndomain))
        use_qdp = .FALSE.
        read(unit=iunit, iostat=ierr) bound_key_dp(0:ndomain)
        if (ierr /= 0) then
            deallocate(bound_key_dp)
            allocate(bound_key_qdp(0:ndomain))
            backspace(iunit)
            read(unit=iunit) bound_key_qdp(0:ndomain)
            use_qdp = .TRUE.
        end if
        !deallocate(bound_key)
    
        return
    end subroutine read_amr_header_internal
    
    subroutine write_amr_header_internal(iunit)
        ! Write header to RAMSES AMR cpu file
        integer, intent(in)                     :: iunit
    
        write(unit=iunit) ncpu
        write(unit=iunit) ndim
        write(unit=iunit) nx, ny, nz
        write(unit=iunit) nlevelmax_cpu
        write(unit=iunit) ngridmax_cpu
        write(unit=iunit) nboundary
        write(unit=iunit) ngrid_current_cpu
        write(unit=iunit) boxlen
        write(unit=iunit) noutput, iout, ifout
        write(unit=iunit) tout(1:noutput)
        write(unit=iunit) aout(1:noutput)
        write(unit=iunit) t
        write(unit=iunit) dtold(1:nlevelmax_cpu)
        write(unit=iunit) dtnew(1:nlevelmax_cpu)
        write(unit=iunit) nstep, nstep_coarse
        write(unit=iunit) const, mass_tot_0, rho_tot
        write(unit=iunit) omega_m, omega_l, omega_k, omega_b, h0,&
                      & aexp_ini, boxlen_ini
        write(unit=iunit) aexp, hexp, aexp_old, epot_tot_int, epot_tot_old
        write(unit=iunit) mass_sph
        
        write(unit=iunit) headl_cpu(1:ncpu, 1:nlevelmax_cpu)
        write(unit=iunit) taill_cpu(1:ncpu, 1:nlevelmax_cpu)
        write(unit=iunit) numbl_cpu(1:ncpu, 1:nlevelmax_cpu)
        write(unit=iunit) numbtot_cpu(1:10, 1:nlevelmax_cpu)
        write(unit=iunit) headf, tailf, numbf, used_mem, used_mem_tot
        write(unit=iunit) ordering
        if (use_qdp) then
            write(unit=iunit) bound_key_qdp
        else
            write(unit=iunit) bound_key_dp
        end if
    
        return
    end subroutine write_amr_header_internal

    subroutine read_amr(filename, icpu)
        ! Read data from RAMSES AMR cpu file
        character(LEN=*), intent(in)            :: filename
        integer, intent(in)                     :: icpu
        
        integer                                 :: iunit
        real (kind=dp), allocatable             :: xdp(:)
        integer, allocatable                    :: iig(:), ind_grid(:)
        integer                                 :: i, ind
        integer                                 :: ilevel, ibound, iskip
        
        iunit = 10 + icpu
        
        open(unit=iunit, file=trim(filename), status='old', form='unformatted')
        
        call read_amr_header_internal(iunit)
        
        allocate(son_cpu(1:ncell))
        if (keep_flag1) allocate(flag1_cpu(0:ncell))
        if (keep_cpu_map) allocate(cpu_map_cpu(1:ncell))
        if (keep_prev) allocate(prev_cpu(1:ngridmax_cpu))
        if (keep_father) allocate(father_cpu(1:ngridmax_cpu))
        if (keep_nbor) allocate(nbor_cpu(1:ngridmax_cpu, 1:twondim))
        if (keep_next) allocate(next_cpu(1:ngridmax_cpu))
        if (keep_xg) allocate(xg_cpu(1:ngridmax_cpu, 1:ndim))
        
        read(unit=iunit) son_cpu(1:ncoarse)
        if (keep_flag1) then
            read(unit=iunit) flag1_cpu(1:ncoarse)
        else
            read(unit=iunit)
        end if
        if (keep_cpu_map) then
            read(unit=iunit) cpu_map_cpu(1:ncoarse)
        else
            read(unit=iunit)
        end if
        do ilevel=1,nlevelmax_cpu
            do ibound=1,nboundary + ncpu
                if (ibound <= ncpu) then
                    ncache = numbl_cpu(ibound, ilevel)
!                     istart = headl_cpu(ibound, ilevel)
                else
                    ncache = numbb(ibound - ncpu, ilevel)
!                     istart = headb(ibound - ncpu, ilevel)
                end if
                if (ncache <= 0) cycle
                allocate(ind_grid(1:ncache))
                allocate(xdp(1:ncache), iig(1:ncache))
                read(unit=iunit) ind_grid
                read(unit=iunit) iig
                if (keep_next) then
                    do i=1,ncache
                        next_cpu(ind_grid(i)) = iig(i)
                    end do
                end if
                read(unit=iunit) iig
                if (keep_prev) then
                    do i=1,ncache
                        prev_cpu(ind_grid(i)) = iig(i)
                    end do
                end if
                do ind=1,ndim
                    read(unit=iunit) xdp
                    if (keep_xg) then
                        do i=1,ncache
                            xg_cpu(ind_grid(i), ind) = xdp(i)
                        end do
                    end if
                end do
                read(unit=iunit) iig
                if (keep_father) then
                    do i=1,ncache
                        father_cpu(ind_grid(i)) = iig(i)
                    end do
                end if
                do ind=1,twondim
                    read(unit=iunit) iig(1:ncache)
                    if (keep_nbor) then
                        do i=1,ncache
                            nbor_cpu(ind_grid(i), ind) = iig(i)
                        end do
                    end if
                end do
                do ind=1,twotondim
                    iskip=ncoarse + (ind-1)*ngridmax_cpu
                    read(unit=iunit) iig
                    do i=1,ncache
                       son_cpu(ind_grid(i) + iskip) = iig(i)
                    end do
                end do
                do ind=1,twotondim
                    iskip=ncoarse + (ind-1)*ngridmax_cpu
                    read(unit=iunit) iig
                    if (keep_cpu_map) then
                        do i=1,ncache
                            cpu_map_cpu(ind_grid(i) + iskip) = iig(i)
                        end do
                    end if
                end do
                do ind=1,twotondim
                    iskip=ncoarse + (ind-1)*ngridmax_cpu
                    read(unit=iunit) iig
                    if (keep_flag1) then
                        do i=1,ncache
                            flag1_cpu(ind_grid(i) + iskip) = iig(i)
                        end do
                    end if
                end do
                deallocate(ind_grid, iig, xdp)
            end do
        end do
        
        close(unit=iunit)
        
        return
    end subroutine read_amr

    subroutine write_amr(filename, icpu)
        ! Read data to RAMSES AMR cpu file
        character(LEN=*), intent(in)            :: filename
        integer, intent(in)                     :: icpu
        
        integer                                 :: iunit
        real (kind=dp), allocatable             :: xdp(:)
        integer, allocatable                    :: iig(:), ind_grid(:)
        integer                                 :: i, ind, igrid, istart
        integer                                 :: ncache
        integer                                 :: ilevel, ibound, iskip
        
        iunit = 10 + icpu
        
        if (.NOT. (keep_father .AND. keep_next .AND. keep_prev .AND. &
                  &keep_flag1 .AND. keep_cpu_map .AND. keep_nbor .AND. &
                  &keep_xg)) then
            stop 'Need all AMR data (except keep_level) to make AMR cpu file!'
        end if
        
        open(unit=iunit, file=trim(filename), status='unknown', form='unformatted')
        
        call write_amr_header_internal(iunit)
        
        write(unit=iunit) son_cpu(1:ncoarse)
        write(unit=iunit) flag1_cpu(1:ncoarse)
        write(unit=iunit) cpu_map_cpu(1:ncoarse)
        
        do ilevel=1,nlevelmax_cpu
            do ibound=1,nboundary + ncpu
                if (ibound <= ncpu) then
                    ncache = numbl_cpu(ibound, ilevel)
                    istart = headl_cpu(ibound, ilevel)
                else
                    ncache = numbb(ibound - ncpu, ilevel)
                    istart = headb(ibound - ncpu, ilevel)
                end if
                if (ncache <= 0) cycle
                allocate(ind_grid(1:ncache))
                allocate(xdp(1:ncache), iig(1:ncache))
                igrid=istart
                do i=1,ncache
                    ind_grid(i)=igrid
                    igrid=next_cpu(igrid)
                end do
                write(unit=iunit) ind_grid

                do i=1,ncache
                    iig(i) = next_cpu(ind_grid(i))
                end do
                write(unit=iunit) iig

                do i=1,ncache
                    iig(i) = prev_cpu(ind_grid(i))
                end do
                write(unit=iunit) iig

                do ind=1,ndim
                    do i=1,ncache
                        xdp(i) = xg_cpu(ind_grid(i), ind)
                    end do
                    write(unit=iunit) xdp
                end do

                do i=1,ncache
                    iig(i) = father_cpu(ind_grid(i))
                end do
                write(unit=iunit) iig
                
                do ind=1,twondim
                    do i=1,ncache
                        iig(i) = nbor_cpu(ind_grid(i), ind)
                    end do
                    write(unit=iunit) iig
                end do
                
                do ind=1,twotondim
                    iskip=ncoarse+(ind-1)*ngridmax
                    do i=1,ncache
                        iig(i) = son_cpu(ind_grid(i) + iskip)
                    end do
                    write(unit=iunit) iig
                end do
                
                do ind=1,twotondim
                    iskip=ncoarse+(ind-1)*ngridmax
                    do i=1,ncache
                        iig(i) = cpu_map_cpu(ind_grid(i) + iskip)
                    end do
                    write(unit=iunit) iig
                end do
                
                do ind=1,twotondim
                    iskip=ncoarse+(ind-1)*ngridmax
                    do i=1,ncache
                        iig(i) = flag1_cpu(ind_grid(i) + iskip)
                    end do
                    write(unit=iunit) iig
                end do
                
                deallocate(ind_grid, iig, xdp)
            end do
        end do
        
        close(unit=iunit)
        
        return
    end subroutine write_amr
    
    subroutine read_hydro_header(filename, icpu)
        ! Read header from RAMSES AMR hydro cpu file
        ! (calls read_hydro_header_internal)
        
        character (LEN=*), intent(in)           :: filename
        integer, intent(in)                     :: icpu
        
        integer                                 :: iunit
        
        iunit = 100000 + icpu
        
        open(unit=iunit, file=trim(filename), status='old', form='unformatted')
        call read_hydro_header_internal(iunit)
        close(unit=iunit)
        
        !deallocate(headl_cpu, taill_cpu, numbl_cpu, numbtot_cpu)
        
        return
    end subroutine read_hydro_header
    
    subroutine read_hydro_header_internal(iunit)
        ! Read header from RAMSES AMR hydro cpu file
        integer, intent(in)              :: iunit
        
        integer                          :: ncpu_temp, ndim_temp
        integer                          :: nlevelmax_cpu_temp, nboundary_temp
        
        read(unit=iunit) ncpu_temp
        if (ncpu /= ncpu_temp) stop "ncpu does not match!"
        read(unit=iunit) nvar
        read(unit=iunit) ndim_temp
        if (ndim /= ndim_temp) stop "ndim does not match!"
        read(unit=iunit) nlevelmax_cpu_temp
        if (nlevelmax_cpu /= nlevelmax_cpu_temp) stop "nlevelmax_cpu does not match!"
        read(unit=iunit) nboundary_temp
        if (nboundary /= nboundary_temp) stop "nlevelmax_cpu does not match!"
        read(unit=iunit) gamma
        
        return
    end subroutine read_hydro_header_internal
    
    subroutine write_hydro_header_internal(iunit)
        ! Write header to RAMSES AMR hydro cpu file
        integer, intent(in)              :: iunit
        
        write(unit=iunit) ncpu
        write(unit=iunit) nvar
        write(unit=iunit) ndim
        write(unit=iunit) nlevelmax_cpu
        write(unit=iunit) nboundary
        write(unit=iunit) gamma
        
        return
    end subroutine write_hydro_header_internal
    
    subroutine read_hydro(filename, icpu)
        ! Read data from RAMSES AMR hydro cpu file
        character(LEN=*), intent(in)     :: filename
        integer, intent(in)              :: icpu
        
        integer                          :: iunit
        real (kind=dp), allocatable      :: xdp(:)
        integer, allocatable             :: ind_grid(:)
        integer                          :: i, ind
        integer                          :: ilevel, ibound, ivar
        integer                          :: istart, iskip, igrid
        integer                          :: ilevel_temp, ncache_temp
        
        iunit = 100000 + icpu
        
        if (.NOT. (keep_u .AND. keep_next)) then
            stop "Cannot read hydro data without keeping 'u' and 'next' data!"
        end if
        
        open(unit=iunit, file=trim(filename), status='old', form='unformatted')
        
        call read_hydro_header_internal(iunit)
        
        if (min_ivar == 0) then
            ivar_min_use = 1
        else
            ivar_min_use = min_ivar
        end if
        if (max_ivar == 0) then
            ivar_max_use = nvar
        else
            ivar_max_use = max_ivar
        end if
        
        allocate(u_cpu(1:ncell, ivar_min_use:ivar_max_use))
        
        do ilevel=1,nlevelmax_cpu
            do ibound=1,nboundary + ncpu
                read(unit=iunit) ilevel_temp
                if (ilevel_temp /= ilevel) stop "Inconsistent levels!"
                read(unit=iunit) ncache_temp
                if (ibound <= ncpu) then
                    ncache = numbl_cpu(ibound, ilevel)
                    istart = headl_cpu(ibound, ilevel)
                else
                    ncache = numbb(ibound - ncpu, ilevel)
                    istart = headb(ibound - ncpu, ilevel)
                end if
                if (ncache_temp /= ncache) stop "Inconsistent ncache!"
                if (ncache <= 0) cycle
                allocate(ind_grid(1:ncache), xdp(1:ncache))
                igrid = istart
                do i=1,ncache
                    ind_grid(i)=igrid
                    igrid=next_cpu(igrid)
                end do
                do ind=1,twotondim
                    iskip = ncoarse + (ind-1)*ngridmax_cpu
                    do ivar=1, ivar_min_use-1
                        read(unit=iunit)
                    end do
                    do ivar=ivar_min_use,ivar_max_use
                        read(unit=iunit) xdp
                        do i=1,ncache
                           u_cpu(ind_grid(i)+iskip, ivar) = xdp(i)
                        end do
                    end do
                    do ivar=ivar_max_use+1, nvar
                        read(unit=iunit)
                    end do
                end do
                deallocate(xdp, ind_grid)
            end do
        end do
        
        close(unit=iunit)
    
    end subroutine read_hydro
    
    subroutine write_hydro(filename, icpu)
        ! Write data to RAMSES AMR hydro cpu file
        character(LEN=*), intent(in)     :: filename
        integer, intent(in)              :: icpu
        
        integer                          :: iunit
        real (kind=dp), allocatable      :: xdp(:)
        integer, allocatable             :: ind_grid(:)
        integer                          :: i, ind
        integer                          :: ilevel, ibound, ivar
        integer                          :: istart, iskip, igrid
        
        iunit = 100000 + icpu
        
        if (.NOT. (keep_u .AND. keep_next)) then
            stop "Cannot write hydro data without 'u' and 'next' data!"
        end if
        
        if (ivar_min_use /= 1 .OR. ivar_max_use /= nvar) then
            stop "Using a custom min_ivar or max_ivar &
                 &is probably a bad idea..."
        end if
        
        open(unit=iunit, file=trim(filename), status='unknown', form='unformatted')
        
        call write_hydro_header_internal(iunit)
        
        do ilevel=1,nlevelmax_cpu
            do ibound=1,nboundary + ncpu
                if (ibound <= ncpu) then
                    ncache = numbl_cpu(ibound, ilevel)
                    istart = headl_cpu(ibound, ilevel)
                else
                    ncache = numbb(ibound - ncpu, ilevel)
                    istart = headb(ibound - ncpu, ilevel)
                end if
                write(unit=iunit) ilevel
                write(unit=iunit) ncache
                if (ncache <= 0) cycle
                allocate(ind_grid(1:ncache), xdp(1:ncache))
                igrid = istart
                do i=1,ncache
                    ind_grid(i)=igrid
                    igrid=next_cpu(igrid)
                end do
                do ind=1,twotondim
                    iskip = ncoarse + (ind-1)*ngridmax_cpu
                    do ivar=1,nvar
                        do i=1,ncache
                           xdp(i) = u_cpu(ind_grid(i)+iskip, ivar)
                        end do
                        write(unit=iunit) xdp
                    end do
                end do
                deallocate(xdp, ind_grid)
            end do
        end do
        
        close(unit=iunit)
    
    end subroutine write_hydro
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! File input/output routines (combined tree)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine load_single(filename)
        ! Load a AMR combined tree
        ! Data format:
        ! 1) Logical variables detailing information included
        ! 2) Header of scalar variables
        ! 3) headl, taill and numbl (for use with next/prev)
        ! 4) some of: son, father, next, prev, level, flag1, cpu_map,
        !    nbor, xg, u
        character (LEN=*), intent(in)        :: filename
        integer                              :: maxout_in, maxlevel_in
        integer                              :: sp_in, dp_in, qdp_in
        integer                              :: min_ivar_in, max_ivar_in, ivar
        
        integer                              :: cur_pos, data_iosize
        integer                              :: logicals(10)
        logical :: file_son, file_father, file_next, file_prev, file_level
        logical :: file_flag1, file_cpu_map, file_nbor, file_xg, file_u
        
        open(unit=50, file=trim(filename),&
            &status='old', access='stream')
        
        ! Read logicals which are stored as integers as different compilers
        ! do logicals in different ways...
        read (50) logicals
        file_son = (logicals(1) == 1)
        file_father = (logicals(2) == 1)
        file_next = (logicals(3) == 1)
        file_prev = (logicals(4) == 1)
        file_level = (logicals(5) == 1)
        file_flag1 = (logicals(6) == 1)
        file_cpu_map = (logicals(7) == 1)
        file_nbor = (logicals(8) == 1)
        file_xg = (logicals(9) == 1)
        file_u = (logicals(10) == 1)
        
        read (50) min_ivar_in, max_ivar_in
        ! Verify we have the data we want
        ! (but we could add it later if not, so only a warning)
        if (.NOT. file_son) stop "ERROR - file has no 'son'!"
        if ((.NOT. file_father) .AND. keep_father) &
            &write(6,*) "WARNING - file has no 'father'!"
        if ((.NOT. file_next) .AND. keep_next) &
            &write(6,*) "WARNING - file has no 'next'!"
        if ((.NOT. file_prev) .AND. keep_prev) &
            &write(6,*) "WARNING - file has no 'prev'!"
        if ((.NOT. file_level) .AND. keep_level) &
            &write(6,*) "WARNING - file has no 'level'!"
        if ((.NOT. file_flag1) .AND. keep_flag1) &
            &write(6,*) "WARNING - file has no 'flag1'!"
        if ((.NOT. file_cpu_map) .AND. keep_cpu_map) &
            &write(6,*) "WARNING - file has no 'cpu_map'!"
        if ((.NOT. file_nbor) .AND. keep_nbor) &
            &write(6,*) "WARNING - file has no 'nbor'!"
        if ((.NOT. file_xg) .AND. keep_xg) &
            &write(6,*) "WARNING - file has no 'xg'!"
        if (keep_u) then
            if (file_u) then
                ! If we have not set min/max_ivar, use file settings,
                ! otherwise use what we asked for
                if (min_ivar == 0) min_ivar = min_ivar_in
                if (max_ivar == 0) max_ivar = max_ivar_in
                ! Check what we asked for is actually in the file
                if (min_ivar < min_ivar_in .OR.&
                   &max_ivar > max_ivar_in) then
                    write(6,*) "WARNING - file does not have full&
                                &range of hydro values!"
                    write (6,'(A,I0,A,I0)') &
                        &"File: ", min_ivar_in, " : ", max_ivar_in
                    write (6,'(A,I0,A,I0)') &
                        &"Requested: ", min_ivar, " : ", max_ivar
                end if
            else
                write(6,*) "WARNING - file has no 'u'!"
            end if
        end if
        
        read (50) ndim, maxout_in, maxlevel_in, sp_in, dp_in, qdp_in, nvar
        if (maxout_in /= maxout) stop "Wrong maxout!"
        if (maxlevel_in /= maxlevel) stop "Wrong maxlevel!"
        if (sp_in /= sp) stop "Wrong sp!"
        if (dp_in /= dp) stop "Wrong dp!"
!         if (qdp_in /= qdp) then                  ! we don't use qdp (yet?)
!             write (6,*) "Warning - wrong qdp!"
!         end if
        twondim = 2 * ndim
        twotondim = 2**ndim
        read (50) nlevelmax, ngridmax, free
        call allocate_combined(0, 0, 0)
        read (50) boxlen, tout, aout, t, dtold, dtnew, const, mass_tot_0, &
                 &rho_tot, omega_m, omega_l, omega_k, omega_b, h0, aexp_ini, &
                 &boxlen_ini, aexp, hexp, aexp_old, epot_tot_int, &
                 &epot_tot_old, mass_sph, gamma, unit_l, unit_d, unit_t
        read (50) headl, taill, numbl
        if (file_son) then
            read (50) son
        end if
        if (file_father) then
            if (keep_father) then
                read (50) father
            else
                inquire(unit=50, pos=cur_pos)
                data_iosize = int_iosize * ngridmax
                read (50, pos=cur_pos+data_iosize)
            end if
        end if
        if (file_next) then
            if (keep_next) then
                read (50) next
            else
                inquire(unit=50, pos=cur_pos)
                data_iosize = int_iosize * ngridmax
                read (50, pos=cur_pos+data_iosize)
            end if
        end if
        if (file_prev) then
            if (keep_prev) then
                read (50) prev
            else
                inquire(unit=50, pos=cur_pos)
                data_iosize = int_iosize * ngridmax
                read (50, pos=cur_pos+data_iosize)
            end if
        end if
        if (file_level) then
            if (keep_level) then
                read (50) level
            else
                inquire(unit=50, pos=cur_pos)
                data_iosize = int_iosize * ngridmax
                read (50, pos=cur_pos+data_iosize)
            end if
        end if
        if (file_flag1) then
            if (keep_flag1) then
                read (50) flag1
            else
                inquire(unit=50, pos=cur_pos)
                data_iosize = int_iosize * ngridmax
                read (50, pos=cur_pos+data_iosize)
            end if
        end if
        if (file_cpu_map) then
            if (keep_cpu_map) then
                read (50) cpu_map
            else
                inquire(unit=50, pos=cur_pos)
                data_iosize = int_iosize * ngridmax
                read (50, pos=cur_pos+data_iosize)
            end if
        end if
        if (file_nbor) then
            if (keep_nbor) then
                read (50) nbor
            else
                inquire(unit=50, pos=cur_pos)
                data_iosize = int_iosize * ngridmax * twondim
                read (50, pos=cur_pos+data_iosize)
            end if
        end if
        if (file_xg) then
            if (keep_xg) then
                read (50) xg
            else
                inquire(unit=50, pos=cur_pos)
                data_iosize = dp_iosize * ngridmax * ndim
                read (50, pos=cur_pos+data_iosize)
            end if
        end if
        if (file_u) then
            if (keep_u) then
                do ivar=min_ivar_in, max_ivar
                    if (ivar < min_ivar) then
                        inquire(unit=50, pos=cur_pos)
                        data_iosize = dp_iosize * ngridmax * twotondim
                        read (50, pos=cur_pos+data_iosize)
                    else
                        read (50) u(:, :, ivar)
                    end if
                end do
            end if
        end if
        
        close(50)
        
        return
    end subroutine load_single
    
    subroutine save_single(filename)
        ! Load a AMR combined tree
        ! Data format:
        ! 1) Logical variables detailing information included
        ! 2) Header of scalar variables
        ! 3) headl, taill and numbl (for use with next/prev)
        ! 4) some of: son, father, next, prev, level, flag1, cpu_map,
        !    nbor, xg, u
        integer                              :: ivar
        character (LEN=*), intent(in)        :: filename
        integer                              :: logicals(10)
        
        open(unit=50, file=trim(filename),&
            &status='replace', access='stream')
        
        ! Write logicals, storing as integers as different compilers
        ! do logicals in different ways...
        logicals = 0
        logicals(1) = 1
        if (keep_father) logicals(2) = 1
        if (keep_next) logicals(3) = 1
        if (keep_prev) logicals(4) = 1
        if (keep_level) logicals(5) = 1
        if (keep_flag1) logicals(6) = 1
        if (keep_cpu_map) logicals(7) = 1
        if (keep_nbor) logicals(8) = 1
        if (keep_xg) logicals(9) = 1
        if (keep_u) logicals(10) = 1
        write (50) logicals
        
        write (50) ivar_min_use, ivar_max_use
        
        write (50) ndim, maxout, maxlevel, sp, dp, qdp, nvar
        write (50) nlevelmax, ngridmax, free
        write (50) boxlen, tout, aout, t, dtold, dtnew, const, mass_tot_0, &
                  &rho_tot, omega_m, omega_l, omega_k, omega_b, h0, aexp_ini, &
                  &boxlen_ini, aexp, hexp, aexp_old, epot_tot_int, &
                  &epot_tot_old, mass_sph, gamma, unit_l, unit_d, unit_t
        write (50) headl, taill, numbl
        write (50) son
        if (keep_father) write (50) father
        if (keep_next) write (50) next
        if (keep_prev) write (50) prev
        if (keep_level) write (50) level
        if (keep_flag1) write (50) flag1
        if (keep_cpu_map) write (50) cpu_map
        if (keep_nbor) write (50) nbor
        if (keep_xg) write (50) xg
        if (keep_u) then
            do ivar=ivar_min_use, ivar_max_use
                write (50) u(:, :, ivar)
            end do
        end if
        
        close(50)
        
        return
    end subroutine save_single
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Translation tree storage routines
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine store_translation_table(icpu)
        integer, intent(in)                 :: icpu
        
        integer                             :: iunit
        
        if (lowmem_tables) then
            iunit = 200000 + icpu
            open(unit=iunit, status='scratch', form='unformatted')
            write (unit=iunit) translation_table
            rewind(unit=iunit)
        else
            translation_tables(:, icpu) = translation_table
        end if
    
        return
    end subroutine store_translation_table
    
    subroutine load_translation_table(icpu)
        integer, intent(in)                 :: icpu
        
        integer                             :: iunit
        
        if (lowmem_tables) then
            iunit = 200000 + icpu
            read (unit=iunit) translation_table
        else
            translation_table = translation_tables(:, icpu)
        end if
    
        return
    end subroutine load_translation_table
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate/deallocate routines
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine deallocate_amr()
        ! Deallocate data for AMR cpu file
        deallocate(headl_cpu, taill_cpu, numbl_cpu, numbtot_cpu)
!         deallocate(headb, tailb, numbb)
        if (use_qdp) then
            deallocate(bound_key_qdp)
        else
            deallocate(bound_key_dp)
        end if
        deallocate(son_cpu)
        if (allocated(father_cpu)) deallocate(father_cpu)
        if (allocated(prev_cpu)) deallocate(prev_cpu)
        if (allocated(flag1_cpu)) deallocate(flag1_cpu)
        if (allocated(cpu_map_cpu)) deallocate(cpu_map_cpu)
        if (allocated(nbor_cpu)) deallocate(nbor_cpu)
        if (allocated(next_cpu)) deallocate(next_cpu)
        if (allocated(xg_cpu)) deallocate(xg_cpu)
    
        return
    end subroutine deallocate_amr
    
    subroutine deallocate_hydro()
        ! Deallocate hydro data for AMR hydro cpu file
        if (allocated(u_cpu)) deallocate(u_cpu)
    
        return
    end subroutine deallocate_hydro
    
    subroutine allocate_combined(nlevelmax_in,&
                                &ngrid_current_in, tot_cells_in)
        ! Allocate combined AMR tree
        integer, intent(in) :: nlevelmax_in
        integer, intent(in) :: ngrid_current_in
        integer, intent(in) :: tot_cells_in
        
        integer             :: aux
        
        ! If you pass in nlevelmax_in == 0 and ngrid_current_in == 0, does
        ! not allocate translation tables or prepare arrays (assume you are
        ! going to load data from a file)
    
        if (nlevelmax_in > 0) nlevelmax = nlevelmax_in
        if (ngrid_current_in > 0) ngridmax = tot_cells_in + 1
        
        allocate(headl(1:nlevelmax))
        allocate(taill(1:nlevelmax))
        allocate(numbl(1:nlevelmax))
        allocate(son(1:ngridmax, 1:twotondim))
        if (keep_father) allocate(father(1:ngridmax))
        if (keep_next) allocate(next(1:ngridmax))
        if (keep_prev) allocate(prev(1:ngridmax))
        if (keep_level) allocate(level(1:ngridmax))
        if (keep_flag1) allocate(flag1(1:ngridmax))
        if (keep_cpu_map) allocate(cpu_map(1:ngridmax))
        if (keep_nbor) allocate(nbor(1:ngridmax, 1:twondim))
        if (keep_xg) allocate(xg(1:ngridmax, 1:ndim))
        
        if (min_ivar == 0) then
            ivar_min_use = 1
        else
            ivar_min_use = min_ivar
        end if
        if (max_ivar == 0) then
            ivar_max_use = nvar
        else
            ivar_max_use = max_ivar
        end if
        if (keep_u) allocate(u(1:ngridmax, 1:twotondim,&
                            &ivar_min_use:ivar_max_use))
        
        if (nlevelmax_in == 0 .AND. ngrid_current_in == 0) return
        
        allocate(translation_table(1:ngrid_current_in))
        if (.NOT. lowmem_tables) then
            allocate(translation_tables(1:ngrid_current_in, 1:ncpu))
        end if
        
        numbl = 0
        free = 1
        aux = new_grid(1) ! aux should be 1; we just need to call new_grid once
        son(1:ncoarse, :) = 0
        if (keep_father) father(1) = 0
        if (keep_nbor) nbor(1, :) = 1
        if (keep_xg) xg(1, :) = 0.5
        if (keep_u) u = -1.0
        
        return
    end subroutine allocate_combined
    
    subroutine deallocate_combined()
        ! Deallocate combined AMR tree
        if (allocated(headl)) deallocate(headl)
        if (allocated(taill)) deallocate(taill)
        if (allocated(numbl)) deallocate(numbl)
        if (allocated(son)) deallocate(son)
        if (allocated(father)) deallocate(father)
        if (allocated(next)) deallocate(next)
        if (allocated(prev)) deallocate(prev)
        if (allocated(level)) deallocate(level)
        if (allocated(flag1)) deallocate(flag1)
        if (allocated(cpu_map)) deallocate(cpu_map)
        if (allocated(nbor)) deallocate(nbor)
        if (allocated(xg)) deallocate(xg)
        if (allocated(u)) deallocate(u)
        
        if (allocated(translation_table)) deallocate(translation_table)
        if (allocated(translation_tables)) deallocate(translation_tables)
        
        return
    end subroutine deallocate_combined
    
    subroutine deallocate_translations()
        ! Deallocate the translation tables
        ! (used only in construction of combined AMR tree)
        integer                               :: icpu, iunit
        
        deallocate(translation_table)
        
        if (lowmem_tables) then
            ! Delete temporary translation table files
            do icpu=1,ncpu
                iunit = 200000 + icpu
                close(unit=iunit)
            end do
        else
            deallocate(translation_tables)
        end if
        
        return
    end subroutine deallocate_translations
    
    subroutine deallocate_level_lists()
        ! Deallocate level linked lists for combined AMR tree
        deallocate(headl)
        deallocate(taill)
        deallocate(numbl)
        if (allocated(next)) deallocate(next)
        if (allocated(prev)) deallocate(prev)
    
        return
    end subroutine deallocate_level_lists
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Combined tree construction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine add_tree(icpu)
        ! Add a RAMSES AMR cpu tree to the combined tree by walking
        ! Build up a translation table to allow later addition of neighbours
        integer, intent(in)        :: icpu
        
        integer                    :: ind, iskip, ilevel
        integer                    :: cur_grid_cpu, cur_grid
        integer                    :: file_son_stack(1:twotondim, 1:maxlevel)
        integer                    :: new_son_stack(1:twotondim, 1:maxlevel)
        integer                    :: son_no_stack(1:maxlevel)
        integer                    :: file_son, new_son
!         logical, parameter         :: debug=.FALSE.
        
        if (ncoarse > 1) stop 'We are assuming there is only one&
                              &coarse cell in the simulation'
        
!         if (debug) write (6,*) "**************** NEW TREE ADDITION ****************"
!         if (debug) write (6,*) "cpu ", icpu
        
        ! Loop over son_cpu cells, recursing into cells
        ! If a cell is deeper in the icpu/file tree than the
        ! combined tree, load into the combined tree
        ! Otherwise skip over entries
        ! ALL REFERENCES TO CURGRID ARE GRID NUMBER NOT MEMORY INDEX
        ! (i.e. ADD NCOARSE, and possibly (ind-1)*ngridmax_cpu)
        ilevel = 1
        cur_grid_cpu = 1
        cur_grid = 1
        son_no_stack(1) = 0
        new_son_stack(:, 1) = 0
        translation_table(1) = 1
        
        mainloop: do
!             if (debug) write (6,*) "doing main loop, ilevel=", ilevel
!             if (debug) write (6,*) "cur_grid_cpu = ", cur_grid_cpu
!             if (debug) write (6,*) "cur_grid = ", cur_grid
            if (son_no_stack(ilevel) == 0) then
!                 if (debug) write (6,*) "son_no_stack(ilevel) = 0"
                son_no_stack(ilevel) = 1
                ! Descent
                do ind=1,twotondim
                    iskip = ncoarse + (ind-1)*ngridmax_cpu
                    file_son_stack(ind, ilevel) = son_cpu(cur_grid_cpu+iskip)
!                     if (debug) write (6,*) "explicit: son_cpu(", cur_grid_cpu+iskip, ") = ", son_cpu(cur_grid_cpu+iskip)
                end do
!                 if (debug) write (6,*) "file_son_stack: ", file_son_stack(:, ilevel)
                new_son_stack(1:twotondim, ilevel) = &
                    &son(cur_grid, 1:twotondim)
!                 if (debug) write (6,*) "new_son_stack: ", new_son_stack(:, ilevel)
            end if
!             if (debug) write (6,*) "Begin iteration through son_cpus"
            do ind=son_no_stack(ilevel),twotondim
!                 if (debug) write (6,*) "ind = ", ind
                file_son = file_son_stack(ind, ilevel)
                new_son = new_son_stack(ind, ilevel)
!                 if (debug) write (6,*) "file_son = ", file_son
!                 if (debug) write (6,*) "new_son = ", new_son
                if ((file_son == 0) .AND. (new_son == 0)) then
                    ! Not present in either tree
                    cycle
                else if (file_son == 0) then
                    ! No new information in file tree
                    cycle
                else if (new_son == 0) then
                    ! This child cell is not represented in the current tree
                    ! Add it
                    new_son = new_grid(ilevel+1)
!                     if (debug) write (6,*) "Adding new son_cpu to slot ", new_son
                    son(cur_grid, ind) = new_son
                    new_son_stack(ind, ilevel) = new_son
!                     if (debug) write (6,*) "new_son_stack (",ind,",",ilevel,&
!                         &") becomes: ", new_son_stack(:, ilevel)
                    if (keep_father) father(new_son) = cur_grid
                    if (keep_level) level(new_son) = ilevel + 1
                    if (keep_flag1) flag1(new_son) = flag1_cpu(file_son)
                    if (keep_cpu_map) cpu_map(new_son) = cpu_map_cpu(file_son)
                    if (keep_xg) xg(new_son, :) = xg_cpu(file_son, :)
                end if
!                 if (debug) write (6,*) &
!                     &"translation_table(", file_son, ") = ", new_son
                translation_table(file_son) = new_son
                ! Descend
                cur_grid_cpu = file_son
                cur_grid = new_son
!                 if (debug) write (6,*) "cur_grid_cpu, cur_grid (descending) &
!                                        &become ", cur_grid_cpu, cur_grid
                son_no_stack(ilevel) = ind
!                 if (debug) write (6,*) "son_no_stack(ilevel) = ", son_no_stack(ilevel)
                ilevel = ilevel + 1
                son_no_stack(ilevel) = 0
!                 if (debug) write (6,*) "We descend to level ", ilevel ,"!"
!                 if (debug) write (6,*) "son_no_stack(ilevel) = ", son_no_stack(ilevel)
                cycle mainloop
                
            end do
            ! Ascent
!             if (debug) write (6,*) "We have finished, ascending to level ", ilevel-1
            ilevel = ilevel - 1
            if (ilevel == 0) then
                exit ! All done
            else if (ilevel == 1) then
                ! Special case for first level
                cur_grid_cpu = 1
                cur_grid = 1
            else
!                 if (debug) write (6,*) "son_no_stack(ilevel-1) = ", son_no_stack(ilevel-1)
!                 if (debug) write (6,*) "file_son_stack(:, ilevel-1) = ", file_son_stack(:, ilevel-1)
!                 if (debug) write (6,*) "new_son_stack(:, ilevel-1) = ", new_son_stack(:, ilevel-1)
                cur_grid_cpu = file_son_stack(son_no_stack(ilevel-1), ilevel-1)
                cur_grid = new_son_stack(son_no_stack(ilevel-1), ilevel-1)
!                 if (debug) write (6,*) "cur_grid_cpu, cur_grid (ascending) become ",&
!                             &cur_grid_cpu, cur_grid
            end if
            son_no_stack(ilevel) = son_no_stack(ilevel) + 1
        end do mainloop
        
!         ! TEST WE HAVE ALL CELLS
!         do ilevel=1, nlevelmax_cpu
!             ncache = numbl_cpu(icpu, ilevel)
!             cur_grid_cpu = headl_cpu(icpu, ilevel)
!             do file_son=1,ncache
!                 cur_grid = translation_table(cur_grid_cpu)
!                 if (cur_grid_cpu == 0) then
!                     write (6,*) "========================================"
!                     write (6,*) "We seem to have missed cell ", cur_grid, "!"
!                     write (6,*) "Level ", ilevel
!                     write (6,*) "cpu: ", icpu
!                     write (6,*) "========================================"
!                 end if
!                 do ind=1,twondim
!                     if (nbor_cpu(cur_grid_cpu, ind) == 0) then
!                         write (6,*) "**************************************"
!                         write (6,*) "Grid ", cur_grid_cpu, " has no neighbour (ind = ", ind, ")!"
!                         write (6,*) "**************************************"
!                     end if
!                 end do
!                 cur_grid_cpu = next_cpu(cur_grid_cpu)
!             end do
!         end do

        call store_translation_table(icpu)
        
        return
    end subroutine add_tree
    
    subroutine add_neighbours_hydro(icpu)
        ! Add neighbour information to the combined tree already built
        ! using the translation tables stored earlier
        ! Also add hydro information using translation tables
        integer, intent(in)        :: icpu
        
        integer                    :: ilevel, ivar, i, dir, ind, iskip
        integer                    :: igrid, ncache, newid
        integer                    :: nbor_aux, nbor_grid, nbor_ind
        
        call load_translation_table(icpu)
        
        if (.NOT. (keep_nbor .OR. keep_u)) then
            stop "Cannot add_neighbours_hydro without 'nbor' and/or 'u' data!"
        end if
        
        if (min_ivar == 0) then
            ivar_min_use = 1
        else
            ivar_min_use = min_ivar
        end if
        if (max_ivar == 0) then
            ivar_max_use = nvar
        else
            ivar_max_use = max_ivar
        end if
        
        if (keep_nbor) then
            do dir=1,twondim
                nbor_ind = nbor_cpu(1, dir) / ngridmax_cpu   ! this is actually ind - 1
                nbor_grid = nbor_cpu(1, dir) - (nbor_ind * ngridmax_cpu)
                if (nbor_grid /= 1) stop "Non-periodic!"
            end do
            nbor(1, :) = 0 ! we can't process this anyway
        end if
        
        if (keep_u) then
            do ivar=ivar_min_use,ivar_max_use
                do ind=1,twotondim
                    iskip = ncoarse + (ind-1)*ngridmax_cpu
                    u(1, ind, ivar) = u_cpu(iskip + 1, ivar)
                end do
            end do
        end if
        
        do ilevel=2, nlevelmax_cpu
            ncache = numbl_cpu(icpu, ilevel)
            ! First do the hydro part
            if (keep_u) then
                igrid = headl_cpu(icpu, ilevel)
                do i=1,ncache
                    newid = translation_table(igrid)
                    do ivar=ivar_min_use,ivar_max_use
                        do ind=1,twotondim
                            iskip = ncoarse + (ind-1)*ngridmax_cpu
                            u(newid, ind, ivar) = u_cpu(iskip + igrid, ivar)
                        end do
                    end do
                    igrid = next_cpu(igrid)
                end do
            end if
            ! Then the neighbours
            if (.NOT. keep_nbor) cycle
            igrid = headl_cpu(icpu, ilevel)
            do i=1,ncache
                newid = translation_table(igrid)
                do dir=1,twondim
                    if (newid <= 0) then
                        write (6,*) "Error in add_neighbours_hydro"
                        write (6,*) "igrid = ", igrid
                        write (6,*) "icpu = ", icpu
                        write (6,*) "dir = ", dir
                        write (6,*) "nbor_cpu(igrid, :) = ", nbor_cpu(igrid, :)
                        write (6,*) "ilevel = ", ilevel
                        write (6,*) "translation_table(igrid) = ", &
                            & translation_table(igrid)
                        write (6,*) "ngrid_current_cpu = ", ngrid_current_cpu
                    end if
                    nbor_aux = nbor_cpu(igrid, dir) - ncoarse
                    nbor_ind = nbor_aux / ngridmax_cpu   ! actually ind - 1
                    nbor_grid = nbor_aux - (nbor_ind * ngridmax_cpu)
                    nbor(newid, dir) = &
                        & translation_table(nbor_grid) + &
                        &(nbor_ind * ngridmax)
                end do
                igrid = next_cpu(igrid)
            end do
        end do
        
        return
    end subroutine add_neighbours_hydro
    
    function new_grid(ilevel)
        ! Add a new cell to the combined tree on level ilevel
        ! Update the appropriate linked lists
        integer, intent(in)        :: ilevel
        integer                    :: new_grid
        integer                    :: old_tail
        
        if (numbl(ilevel) == 0) then
            numbl(ilevel) = 1
            headl(ilevel) = free
            taill(ilevel) = free
            if (keep_next) next(free) = 0
            if (keep_prev) prev(free) = 0
            new_grid = free
            free = free + 1
        else
            numbl(ilevel) = numbl(ilevel) + 1
            old_tail = taill(ilevel)
            taill(ilevel) = free
            if (keep_next) next(old_tail) = free
            if (keep_prev) prev(free) = old_tail
            new_grid = free
            free = free + 1
        end if
        
        return
    end function new_grid
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Data retrieval routines
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine dump_level_cells(level_in, iv_min, iv_max, &
                               &dump_cells, leaf_only)
        ! Scan through all the grids on level level_in. Scan all the hydro
        ! quantities for their cells between columns ivar_min_in and
        ! ivar_max_in. Store these in the flat grid dump_grid. Leave gaps from
        ! AMR unchanged. If leaf_only, only change the value if the cell
        ! is a leaf cell.
        integer, intent(in)               :: level_in
        integer, intent(in)               :: iv_min, iv_max
        double precision, intent(inout)   :: dump_cells(:,:,:,:)
        logical, intent(in), optional     :: leaf_only
        
        integer                           :: grid_size, cell_size
        integer                           :: gx, gy, gz
        integer                           :: igrid
        double precision                  :: xg_loc(1:3)
        logical                           :: use_leaf_only
        logical                           :: lt(1:8) ! leaf table
        
        ! just to stop gfortran warnings...
        xg_loc = 0.0
        gy = 0
        gz = 0
        
        use_leaf_only = .FALSE.
        if (present(leaf_only)) then
            if (leaf_only) then
                use_leaf_only = .TRUE.
            end if
        end if
        
        if (.NOT. (keep_next .AND. keep_xg .AND. keep_u)) then
            stop "Cannot dump_level_cells without 'next', 'xg' and 'u' data!"
        end if
        
        if ((iv_min < min_ivar) .OR. (iv_max > max_ivar)) then
            stop "Requested hydro columns outside of those in file!"
        end if
        
        ! We will identify the position of the grids (each with 2/4/8 cells)
        ! and then place the children into dump_cells
        grid_size = 2**(level_in-1)
        cell_size = 2**level_in
        
        igrid = headl(level_in)
        do
            lt = .TRUE.
            xg_loc(1:ndim) = xg(igrid, 1:ndim)
            gx = nint(xg_loc(1) * grid_size + 0.5)
            if (ndim>1) gy = nint(xg_loc(2) * grid_size + 0.5)
            if (ndim>2) gz = nint(xg_loc(3) * grid_size + 0.5)
            
            if (use_leaf_only) then
                lt(1:twotondim) = (son(igrid, 1:twotondim) == 0)
                if (ndim==3) then
                    if (lt(1)) dump_cells(2*gx-1, 2*gy-1, 2*gz-1, :) = u(igrid, 1, iv_min:iv_max)
                    if (lt(2)) dump_cells(2*gx  , 2*gy-1, 2*gz-1, :) = u(igrid, 2, iv_min:iv_max)
                    if (lt(3)) dump_cells(2*gx-1, 2*gy  , 2*gz-1, :) = u(igrid, 3, iv_min:iv_max)
                    if (lt(4)) dump_cells(2*gx  , 2*gy  , 2*gz-1, :) = u(igrid, 4, iv_min:iv_max)
                    if (lt(5)) dump_cells(2*gx-1, 2*gy-1, 2*gz  , :) = u(igrid, 5, iv_min:iv_max)
                    if (lt(6)) dump_cells(2*gx  , 2*gy-1, 2*gz  , :) = u(igrid, 6, iv_min:iv_max)
                    if (lt(7)) dump_cells(2*gx-1, 2*gy  , 2*gz  , :) = u(igrid, 7, iv_min:iv_max)
                    if (lt(8)) dump_cells(2*gx  , 2*gy  , 2*gz  , :) = u(igrid, 8, iv_min:iv_max)
                
                else if (ndim==2) then
                    if (lt(1)) dump_cells(1, 2*gx-1, 2*gy-1, :) = u(igrid, 1, iv_min:iv_max)
                    if (lt(2)) dump_cells(1, 2*gx  , 2*gy-1, :) = u(igrid, 2, iv_min:iv_max)
                    if (lt(3)) dump_cells(1, 2*gx-1, 2*gy  , :) = u(igrid, 3, iv_min:iv_max)
                    if (lt(4)) dump_cells(1, 2*gx  , 2*gy  , :) = u(igrid, 4, iv_min:iv_max)
                
                else
                    if (lt(1)) dump_cells(1, 1, 2*gx-1, :) = u(igrid, 1, iv_min:iv_max)
                    if (lt(2)) dump_cells(1, 1, 2*gx  , :) = u(igrid, 2, iv_min:iv_max)
                end if
            else
                if (ndim==3) then
                    dump_cells(2*gx-1, 2*gy-1, 2*gz-1, :) = u(igrid, 1, iv_min:iv_max)
                    dump_cells(2*gx  , 2*gy-1, 2*gz-1, :) = u(igrid, 2, iv_min:iv_max)
                    dump_cells(2*gx-1, 2*gy  , 2*gz-1, :) = u(igrid, 3, iv_min:iv_max)
                    dump_cells(2*gx  , 2*gy  , 2*gz-1, :) = u(igrid, 4, iv_min:iv_max)
                    dump_cells(2*gx-1, 2*gy-1, 2*gz  , :) = u(igrid, 5, iv_min:iv_max)
                    dump_cells(2*gx  , 2*gy-1, 2*gz  , :) = u(igrid, 6, iv_min:iv_max)
                    dump_cells(2*gx-1, 2*gy  , 2*gz  , :) = u(igrid, 7, iv_min:iv_max)
                    dump_cells(2*gx  , 2*gy  , 2*gz  , :) = u(igrid, 8, iv_min:iv_max)
                
                else if (ndim==2) then
                    dump_cells(1, 2*gx-1, 2*gy-1, :) = u(igrid, 1, iv_min:iv_max)
                    dump_cells(1, 2*gx  , 2*gy-1, :) = u(igrid, 2, iv_min:iv_max)
                    dump_cells(1, 2*gx-1, 2*gy  , :) = u(igrid, 3, iv_min:iv_max)
                    dump_cells(1, 2*gx  , 2*gy  , :) = u(igrid, 4, iv_min:iv_max)
                
                else
                    dump_cells(1, 1, 2*gx-1, :) = u(igrid, 1, iv_min:iv_max)
                    dump_cells(1, 1, 2*gx  , :) = u(igrid, 2, iv_min:iv_max)
                end if
            end if
            
            if (igrid == taill(level_in)) exit
            igrid = next(igrid)
        end do
        
        return
    end subroutine dump_level_cells
    
    subroutine get_hydro(n_cells, ivar, cells)
        ! Return all hydro cells for desired hydro variable/column
        ! Must be called with correct number of cells in tree
        integer, intent(in)               :: n_cells
        integer, intent(in)               :: ivar
        double precision, intent(out)     :: cells(1:n_cells)
        
        integer                           :: ilevel, ind, i
        integer                           :: igrid, ncache, slot
        
        if (.NOT. (keep_next .AND. keep_u)) then
            stop "Cannot get_hydro without 'next' and 'u' data!"
        end if
        
        slot = 1
        do ilevel=1, nlevelmax
            igrid = headl(ilevel)
            ncache = numbl(ilevel)
            do i=1,ncache
                do ind=1,twotondim
                    cells(slot) = u(igrid, ind, ivar)
                    slot = slot + 1
                end do
                igrid = next(igrid)
            end do
        end do
        
        return
    end subroutine get_hydro
    
    subroutine grid_min_max(xg_loc, ilevel, min_r, max_r)
        ! Returns to limits of a hypothetical cell
        ! centred at xg_in on level ilevel
        double precision, intent(in)      :: xg_loc(1:ndim)
        integer, intent(in)               :: ilevel
        double precision, intent(out)     :: min_r(1:ndim, 0:twotondim)
        double precision, intent(out)     :: max_r(1:ndim, 0:twotondim)
        
        integer                           :: ind, ix, iy, iz
        double precision                  :: dx, xgc(1:3), xc(1:3, 1:twotondim)
        
        dx = 0.5**ilevel ! Spacing at CELL level (ilevel+1), not GRID level
        
        do ind=1,twotondim
            iz = (ind-1)/4
            iy = (ind-1-4*iz)/2
            ix = (ind-1-2*iy-4*iz)
            if (ndim>0) xc(1, ind) = (dble(ix) - 0.5D0) * dx
            if (ndim>1) xc(2, ind) = (dble(iy) - 0.5D0) * dx
            if (ndim>2) xc(3, ind) = (dble(iz) - 0.5D0) * dx
        end do
        
        min_r(1:ndim, 0) = xg_loc - dx
        max_r(1:ndim, 0) = xg_loc + dx
        
        do ind=1,twotondim
            xgc = xg_loc + xc(:, ind)
            min_r(:, ind) = xgc - 0.5*dx
            max_r(:, ind) = xgc + 0.5*dx
        end do
        
        return
    end subroutine grid_min_max
    
    function cell_position(igrid, ind, ilevel)
        ! Return the centre position of cell (igrid, ind)
        ! on level ilevel
        double precision                  :: cell_position(1:ndim)
        
        integer, intent(in)               :: igrid, ind
        integer, intent(in)               :: ilevel ! GRID level?
        
        integer                           :: ix, iy, iz
        double precision                  :: xg_loc(1:3)
        double precision                  :: dx, xc(1:3, 1:twotondim)
        
        if (.NOT. keep_xg) then
            stop "Cannot calculate cell_position without 'xg' data!"
        end if
        
        xg_loc(1:ndim) = xg(igrid, :)

        dx = 0.5**ilevel ! Spacing at CELL level (ilevel+1), not GRID level
        
        iz = (ind-1)/4
        iy = (ind-1-4*iz)/2
        ix = (ind-1-2*iy-4*iz)
        if (ndim>0) xc(1, ind) = (dble(ix) - 0.5D0) * dx
        if (ndim>1) xc(2, ind) = (dble(iy) - 0.5D0) * dx
        if (ndim>2) xc(3, ind) = (dble(iz) - 0.5D0) * dx
        
        cell_position = xg_loc(1:ndim) + xc(1:ndim, ind)
        
        return
    end function cell_position
    
    subroutine column_density(ivar, resolution, min_z, max_z, map)
        ! Make a z-axis column density plot of hydro variable ivar
        integer, intent(in)              :: ivar
        integer, intent(in)              :: resolution
        double precision, intent(in)     :: min_z, max_z
        double precision, intent(out)    :: map(1:resolution, 1:resolution)
        
        logical                          :: use_limits
        integer                          :: ilevel, ind, icache, i, j
        integer                          :: igrid, ncache
        double precision                 :: xg_loc(1:ndim)
        double precision                 :: dens, line_length
        double precision                 :: contrib, frac(1:2)
        double precision                 :: min_r(1:ndim, 0:twotondim)
        double precision                 :: max_r(1:ndim, 0:twotondim)
        integer                          :: min_grid(1:2), max_grid(1:2)
        
        if (.NOT. keep_xg) then
            stop "Cannot calculate column_density without 'xg' data!"
        end if
        
        map = 0.0
        use_limits = (min_z > 0.0 .OR. (max_z > 0.0 .AND. max_z < 1.0))
        
        do ilevel=1, nlevelmax
!             write (6,*) "ilevel = ", ilevel
            igrid = headl(ilevel)
            ncache = numbl(ilevel)
            do icache=1,ncache
                xg_loc = xg(igrid, :)
                call grid_min_max(xg_loc, ilevel, min_r, max_r)
!                 write (6,*) "igrid = ", igrid
!                 write (6,*) "xg_loc = ", xg_loc
                min_r(1:2,:) = min_r(1:2,:) * resolution
                max_r(1:2,:) = max_r(1:2,:) * resolution
                do ind=1,twotondim
                    if (son(igrid, ind) == 0) then
                        ! Leaf cell
                        if (use_limits) then
                            if (min_r(3,ind) > max_z .OR.&
                                &max_r(3,ind) < min_z) cycle
                        end if
                        dens = u(igrid, ind, ivar)
                        line_length = max_r(3, ind) - min_r(3,ind)
                        contrib = dens * line_length
                        min_grid(1:2) = floor(min_r(1:2, ind))
!                         min_grid(1:2) = max(min_grid(1:2), (/1, 1/))
!                         min_grid(1:2) = min(min_grid(1:2),&
!                                            &(/resolution, resolution/))
                        max_grid(1:2) = ceiling(max_r(1:2, ind))
!                         max_grid(1:2) = max(max_grid(1:2), (/1, 1/))
!                         max_grid(1:2) = min(max_grid(1:2),&
!                                            &(/resolution, resolution/))
                        do i=min_grid(1)+1, max_grid(1)
                            if (i==min_grid(1)+1 .AND. i==max_grid(1)) then
                                frac(1) = max_r(1, ind) - min_r(1, ind)
                            else if (i==min_grid(1)+1) then
                                frac(1) = dble(min_grid(1)+1) - min_r(1, ind)
                            else if (i==max_grid(1)) then
                                frac(1) = max_r(1, ind) - dble(max_grid(1)-1)
                            else
                                frac(1) = 1.0d0
                            end if
                            do j=min_grid(2)+1, max_grid(2)
                                if (j==min_grid(2)+1 .AND. j==max_grid(2)) then
                                    frac(2) = max_r(2, ind) - min_r(2, ind)
                                else if (j==min_grid(2)+1) then
                                    frac(2) = dble(min_grid(2)+1) - min_r(2, ind)
                                else if (j==max_grid(2)) then
                                    frac(2) = max_r(2, ind) - dble(max_grid(2)-1)
                                else
                                    frac(2) = 1.0d0
                                end if
                                map(i, j) = map(i, j) + contrib * product(frac)
                            end do
                        end do
                    end if
                end do
                igrid = next(igrid)
            end do
        end do
        
        return
    end subroutine column_density
    
    function get_left_neighbours(igrid, ind)
        ! Return the -x, -y, -z neighbouring cells (igrid, ind format)
        ! of cell (igrid, ind), accounting for refinement
        integer                        :: get_left_neighbours(1:12)
        
        integer, intent(in)            :: igrid, ind
        
        integer                        :: neighbours(1:2, 1:12)
        integer                        :: parent_nbors(1:twondim)
        integer                        :: parent_nbors_ind(1:twondim)
        integer                        :: i, n_neib, nbor_grid, nbor_ind
        integer                        :: ix, iy, iz
        
        ! We want the 3-12 neighbours of cell (igrid, ind) in the -x, -y and -z
        ! directions (number of neighbours depends on refinement)
        
        if (.NOT. keep_nbor) then
            stop "Cannot get_neighbours without 'nbor' data!"
        end if
        
        get_left_neighbours = 0
        neighbours = 0
        n_neib = 0
        
        if (igrid == 1) stop "no grid?" ! no neighbours for initial cell
        
        parent_nbors = nbor(igrid, 1:twondim)
        do i=1,twondim
            parent_nbors_ind(i) = (parent_nbors(i) / ngridmax) + 1
            parent_nbors(i) = parent_nbors(i) - &
                & (parent_nbors_ind(i)-1) * ngridmax
        end do
        
        iz = (ind-1)/4
        iy = (ind-1-4*iz)/2
        ix = (ind-1-2*iy-4*iz)
        
        ! -x neighbour
        if (ix == 0) then
            ! Look to left (-x) grid
            nbor_grid = parent_nbors(1)
            nbor_ind = parent_nbors_ind(1)
            ! WHAT THIS MEANS - nbor_grid is the neighbouring grid of our PARENT
            ! grid; nbor_grid cell nbor_ind is the neighbouring cell of the SAME
            ! refinement level as our GRID. We check if this cell is refined,
            ! in which case we take the appropriate cell or cells of this grid,
            ! otherwise we just use the (larger) grid/ind as a neighbour
            if (son(nbor_grid, nbor_ind) == 0) then
                ! not refined, use over-sized cell as neighbour
                n_neib = n_neib + 1
                neighbours(1:2, n_neib) = (/nbor_grid, nbor_ind/)
            else
                ! refined; see if singly or doubly refined
                nbor_grid = son(nbor_grid, nbor_ind)
                nbor_ind = ind + 1 ! ix = ix + 1 for opposite cell to ours
                if (son(nbor_grid, nbor_ind) == 0) then
                    ! same refinement as us
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, nbor_ind/)
                else
                    ! more refined than us; four boxes to take
                    nbor_grid = son(nbor_grid, nbor_ind)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 2/)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 4/)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 6/)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 8/)
                end if
            end if
            
        else
            ! Neighbour is within own oct
            n_neib = n_neib + 1
            neighbours(1:2, n_neib) = (/igrid, ind-1/)
        end if
        
        ! -y neighbour
        if (iy == 0) then
            ! Look to left (-y) grid
            nbor_grid = parent_nbors(3)
            nbor_ind = parent_nbors_ind(3)
            ! WHAT THIS MEANS - nbor_grid is the neighbouring grid of our PARENT
            ! grid; nbor_grid cell nbor_ind is the neighbouring cell of the SAME
            ! refinement level as our GRID. We check if this cell is refined,
            ! in which case we take the appropriate cell or cells of this grid,
            ! otherwise we just use the (larger) grid/ind as a neighbour
            if (son(nbor_grid, nbor_ind) == 0) then
                ! not refined, use over-sized cell as neighbour
                n_neib = n_neib + 1
                neighbours(1:2, n_neib) = (/nbor_grid, nbor_ind/)
            else
                ! refined; see if singly or doubly refined
                nbor_grid = son(nbor_grid, nbor_ind)
                nbor_ind = ind + 2 ! iy = iy + 1 for opposite cell to ours
                if (son(nbor_grid, nbor_ind) == 0) then
                    ! same refinement as us
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, nbor_ind/)
                else
                    ! more refined than us; four boxes to take
                    nbor_grid = son(nbor_grid, nbor_ind)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 3/)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 4/)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 7/)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 8/)
                end if
            end if
        else
            ! Neighbour is within own oct
            n_neib = n_neib + 1
            neighbours(1:2, n_neib) = (/igrid, ind-2/)
        end if
        
        ! -z neighbour
        if (iz == 0) then
            ! Look to left (-z) grid
            nbor_grid = parent_nbors(5)
            nbor_ind = parent_nbors_ind(5)
            ! WHAT THIS MEANS - nbor_grid is the neighbouring grid of our PARENT
            ! grid; nbor_grid cell nbor_ind is the neighbouring cell of the SAME
            ! refinement level as our GRID. We check if this cell is refined,
            ! in which case we take the appropriate cell or cells of this grid,
            ! otherwise we just use the (larger) grid/ind as a neighbour
            if (son(nbor_grid, nbor_ind) == 0) then
                ! not refined, use over-sized cell as neighbour
                n_neib = n_neib + 1
                neighbours(1:2, n_neib) = (/nbor_grid, nbor_ind/)
            else
                ! refined; see if singly or doubly refined
                nbor_grid = son(nbor_grid, nbor_ind)
                nbor_ind = ind + 4 ! iz = iz + 1 for opposite cell to ours
                if (son(nbor_grid, nbor_ind) == 0) then
                    ! same refinement as us
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, nbor_ind/)
                else
                    ! more refined than us; four boxes to take
                    nbor_grid = son(nbor_grid, nbor_ind)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 5/)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 6/)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 7/)
                    n_neib = n_neib + 1
                    neighbours(1:2, n_neib) = (/nbor_grid, 8/)
                end if
            end if
        else
            ! Neighbour is within own oct
            n_neib = n_neib + 1
            neighbours(1:2, n_neib) = (/igrid, ind-4/)
        end if
        
        do i=1,n_neib
            get_left_neighbours(i) = neighbours(1,i) + ngridmax*neighbours(2,i)
        end do
        
        return
    end function get_left_neighbours
    
end module amr_utils
