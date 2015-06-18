module stats_scratch_data
    use quickselect_mod, calc_median => median
    implicit none

    type stats_type
        double precision              :: mean
        double precision              :: median
        double precision              :: std_dev
        double precision              :: mad
    end type stats_type
    
    type sample_data_type
        integer                       :: n
        double precision, allocatable :: values(:)
        double precision, allocatable :: weights(:)
        double precision, allocatable :: scratch(:)
    
        contains
        procedure :: add => add_sample
        procedure :: alloc => allocate_sample
        procedure :: reset => reset_sample
        procedure :: stats => stats_sample
    end type sample_data_type
    
    type results_type
        type(sample_data_type) :: means
        type(sample_data_type) :: medians
        type(sample_data_type) :: std_devs
        type(sample_data_type) :: mads
        
        contains
        procedure :: add => add_results
        procedure :: alloc => allocate_results
    end type results_type
    
    contains
    
    ! sample_data methods
    subroutine add_sample(self, value, weight)
        class(sample_data_type), intent(inout) :: self
        double precision, intent(in)           :: value
        double precision, intent(in)           :: weight
        self%n = self%n + 1
        self%values(self%n) = value
        self%weights(self%n) = weight
    end subroutine add_sample
    
    subroutine allocate_sample(self, sample_size)
        class(sample_data_type), intent(inout) :: self
        integer, intent(in)                    :: sample_size
        if (allocated(self%values)) deallocate(self%values)
        if (allocated(self%weights)) deallocate(self%weights)
        if (allocated(self%scratch)) deallocate(self%scratch)
        allocate(self%values(1:sample_size))
        allocate(self%weights(1:sample_size))
        allocate(self%scratch(1:sample_size))
        self%n = 0
    end subroutine allocate_sample
    
    subroutine reset_sample(self)
        class(sample_data_type), intent(inout) :: self
        self%n = 0
    end subroutine reset_sample
    
    function stats_sample(self)
        class(sample_data_type), intent(inout) :: self
        type(stats_type)                       :: stats_sample
        
        double precision                       :: mean, median, std_dev, mad
        
        integer                                :: M, N
        double precision                       :: bias, sum_weights, aux

        N = self%n
        sum_weights = sum(self%weights(1:N))

        ! mean
        mean = sum(self%values(1:N)*self%weights(1:N)) / sum_weights
        
        ! standard deviation
        M = count(self%weights(1:N) > tiny(1.d0))
        if (M > 1) then
            aux = sum(self%weights(1:N) * (self%values(1:N) - mean)**2)
            bias = (dble(M)-1.d0) / dble(M)
            std_dev = sqrt(aux / (bias * sum_weights))
        else
            std_dev = 0.0d0
        end if
        
        ! median
        self%scratch(1:N) = self%values(1:N)
        median = calc_median(self%scratch(1:N))
        
        ! median absolute deviation from the median
        self%scratch(1:N) = abs(self%values - median)
        mad = calc_median(self%scratch(1:N))
        
        stats_sample%mean = mean
        stats_sample%median = median
        stats_sample%std_dev = std_dev
        stats_sample%mad = mad
    end function stats_sample
    
    ! results methods
    subroutine add_results(self, stats)
        class(results_type), intent(inout)  :: self
        type(stats_type), intent(in)        :: stats
        call self%means%add(stats%mean, 1.0d0)
        call self%std_devs%add(stats%std_dev, 1.0d0)
        call self%medians%add(stats%median, 1.0d0)
        call self%mads%add(stats%mad, 1.0d0)
        
        return
    end subroutine add_results
    
    subroutine allocate_results(self, sample_size)
        class(results_type), intent(inout)  :: self
        integer, intent(in)                 :: sample_size
        call self%means%alloc(sample_size)
        call self%medians%alloc(sample_size)
        call self%std_devs%alloc(sample_size)
        call self%mads%alloc(sample_size)
        return
    end subroutine allocate_results
    
end module stats_scratch_data

program velocity_dispersion
    use stats_scratch_data
    use amr_utils
    implicit none
    
    double precision, parameter    :: pc=3.0856776d+18      ! pc in cm
    double precision, parameter    :: km_s=100000.d0        ! km/s in cm/s
    double precision               :: v_scale
    
    character (LEN=:), allocatable :: filename
    character (LEN=5)              :: suffix
    character (LEN=:), allocatable :: output_filename
    integer                        :: length
    
    integer                        :: ilevel, iilevel, ind
    integer                        :: max_level, max_grids, sample_size
    integer, allocatable           :: son_stack(:)
    integer, allocatable           :: ind_stack(:)
    
    type(sample_data_type), allocatable :: sample_stack(:)
    type(results_type), allocatable     :: results_stack(:)
    type(stats_type)                    :: stats
    type(stats_type)                    :: mean_stats, std_dev_stats
    type(stats_type)                    :: median_stats, mad_stats
    
    double precision                    :: vmag, weight
    
    double precision, allocatable       :: dx_stack(:)
    
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
    
    if (command_argument_count() /= 1) then
        stop "Wrong number of command line arguments (expecting 1)!"
    end if
    
    call get_command_argument(1, length=length)
    allocate(character(LEN=length) :: filename)
    call get_command_argument(1, filename)
    
    call load_single(filename)
    
    v_scale = (unit_l / unit_t) / km_s
    
    allocate(son_stack(1:nlevelmax))
    allocate(ind_stack(1:nlevelmax))
    
    ! Find out how many levels are actually being used
    do ilevel=nlevelmax,1,-1
        if (numbl(ilevel) /= 0) exit
    end do
    max_level = ilevel
    max_grids = sum(numbl)
    
    ! Allocate space on each level for the sampling data
    allocate(sample_stack(1:max_level))
    do ilevel=1,max_level
        sample_size = twotondim**(max_level - ilevel + 1)
        sample_size = min(max_grids*twotondim, sample_size)
        call sample_stack(ilevel)%alloc(sample_size)
    end do
    
    ! Allocate space on each level for the statistics results
    allocate(results_stack(1:max_level))
    do ilevel=1,max_level
        sample_size = twotondim**(ilevel-1)
        sample_size = min(max_grids, sample_size)
        call results_stack(ilevel)%alloc(sample_size)
    end do
    
    ! Loop through cells, collecting data and calculating statistics
    ilevel = 1
    son_stack(1) = 1
    ind_stack(1) = 0
    do
        ! Main loop over indices 1->8 (3D)
        ind_stack(ilevel) = ind_stack(ilevel) + 1
        ind = ind_stack(ilevel)
        if (ind>twotondim) then
            ! We have finished this grid on this level
            ! Ascend - first calculate statistics and reset sample data
            stats = sample_stack(ilevel)%stats()
            call results_stack(ilevel)%add(stats)
            call sample_stack(ilevel)%reset()
            ilevel = ilevel - 1
            if (ilevel == 0) then
                ! We are done
                exit
            end if
    !v_scale = 1.d0
            cycle
        end if
        if (son(son_stack(ilevel), ind) == 0) then
            ! We have a leaf cell, add data
            weight = 2.0**(nlevelmax - ilevel)
            vmag = v_scale * sqrt(dot_product(u(son_stack(ilevel), ind, 2:4), &
                                 &u(son_stack(ilevel), ind, 2:4)))
            ! Add this data to the sample data on every level less than
            ! or equal to the level of this cell
            do iilevel=1,ilevel
                call sample_stack(iilevel)%add(vmag, weight)
            end do
        else
            ! We need to descend a level
            son_stack(ilevel + 1) = son(son_stack(ilevel), ind_stack(ilevel))
            ilevel = ilevel + 1
!             write (6,*) "Descending to ilevel = ", ilevel
            ind_stack(ilevel) = 0
        end if
    end do
    
    ! Final averaging/write results
    allocate(dx_stack(1:max_level))
    
    do ilevel=1,max_level
        dx_stack(ilevel) = dx_level(ilevel) * boxlen * unit_l / pc
    end do
    
    suffix = get_single_suffix(filename)
    
    if (suffix == "") then
        output_filename = "vel_stats.dat"
    else
        output_filename = "vel_stats_"//suffix//".dat"
    end if

    open(1, file=output_filename)
    
    ! NOTE THAT THE MEDIAN BASED STATISTICS ARE NOT WEIGHTED BY GRID SIZE
    ! (this should make no difference for fixed grids)
    ! 'MAD' = Median Absolute Difference (from the median)
    
    ! Results file format:
    ! Each line represents a level
    ! Columns are: dx, mean of means, std_dev of means, mean of std_devs,
    ! std_dev of std_devs, median of medians, mad of medians, median of mads,
    ! mad of mads
    write (1, '(A)') '# dx; mean of means; std dev of means; mean of std devs;&
                    & std dev of std devs; median of medians; median of mads;&
                    & mad of mads'
    
    do ilevel=1,max_level
        mean_stats = results_stack(ilevel)%means%stats()
        std_dev_stats = results_stack(ilevel)%std_devs%stats()
        median_stats = results_stack(ilevel)%medians%stats()
        mad_stats = results_stack(ilevel)%mads%stats()
        write (1, '(9E17.10)') &
                    &dx_stack(ilevel), &         ! 0 this is important
                    &mean_stats%mean, &          ! 1
                    &mean_stats%std_dev, &       ! 2
                    &std_dev_stats%mean, &       ! 3 either this
                    &std_dev_stats%std_dev, &    ! 4
                    &median_stats%median, &      ! 5
                    &median_stats%mad, &         ! 6
                    &mad_stats%median, &         ! 7 or this is also important
                    &mad_stats%mad               ! 8
    end do
    
    close(1)
    
    contains
    
    function dx_level(ilevel)
        integer, intent(in)      :: ilevel
        double precision         :: dx_level
        
        dx_level = 2.0**(1-ilevel)
    end function dx_level

end program velocity_dispersion