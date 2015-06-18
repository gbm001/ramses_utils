program test_binning
    use binning
    implicit none
    
    double precision, allocatable        :: x3D(:,:,:)
    double precision, allocatable        :: y3D(:,:,:)
    double precision                     :: P, kmag, R
    
    integer                              :: i, j, k, l
    integer                              :: grid_size, Nbins
    integer                              :: gmin, gmax
    
    double precision, allocatable        :: x1D(:)
    double precision, allocatable        :: y1D(:)
    
    integer, allocatable                 :: counts(:)
    double precision, allocatable        :: bin_avg(:)
    double precision, allocatable        :: bin_centres(:)
    
    double precision, allocatable        :: tr_bin_avg(:)
    
    grid_size = 64
    gmin = -(grid_size-1) / 2
    gmax = grid_size / 2
    Nbins = ceiling(sqrt(3.0*(64**2)))
    allocate(x3D(gmin:gmax, gmin:gmax, gmin:gmax))
    allocate(y3D(gmin:gmax, gmin:gmax, gmin:gmax))
    
    do k=gmin,gmax
        do j=gmin,gmax
            do i=gmin,gmax
                kmag = sqrt(dble(i**2 + j**2 + k**2))
                P = max(1.d0, kmag)**(-4)
                call random_number(R)
                R = (R * 2.0) - 1.0
                
                x3D(i,j,k) = kmag
                y3D(i,j,k) = P * (1.d0 + R)
            end do
        end do
    end do
    
    x1D = reshape(x3D, (/grid_size**3/))
    y1D = reshape(y3D, (/grid_size**3/))
    
    allocate(counts(0:Nbins))
    allocate(bin_avg(0:Nbins))
    allocate(bin_centres(0:Nbins))
    
    allocate(tr_bin_avg(0:Nbins))
    
    !call square_binning(x1D, y1D, Nbins, counts, bin_avg, bin_centres)
    call square_binning(y3D, Nbins, counts, bin_avg, bin_centres)
    
!     call d_triangle_binning(x1D, y1D, Nbins, tr_bin_avg)
    call triangle_binning(y3D, Nbins, tr_bin_avg)
    
    open(10, file='test_binning.dat', form='formatted')
    
    write (10, *) '# bin number : {bin middle : counts : bin average : bin centre actual } '&
                 &'(square binning) : {bin average} (triangular binning)'
    
    do l=0,Nbins
        write (10, *) l, dble(l+0.5), counts(l), bin_avg(l), bin_centres(l), tr_bin_avg(l)
    end do

    close(10)
    
    stop
end program test_binning
