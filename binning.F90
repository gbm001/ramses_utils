module binning
    implicit none

    interface square_binning
        module procedure d_square_binning
        module procedure d_square_binning1D
        module procedure d_square_binning2D
        module procedure d_square_binning3D
    end interface square_binning

    interface triangle_binning
        module procedure d_triangle_binning
        module procedure d_triangle_binning1D
        module procedure d_triangle_binning2D
        module procedure d_triangle_binning3D
    end interface triangle_binning

    contains
    
    subroutine d_square_binning(x, y, Nbins, counts, bin_avg, bin_centres)
        double precision, intent(in)  :: x(:)
        double precision, intent(in)  :: y(:)
        integer, intent(in)           :: Nbins
        
        integer, intent(out)          :: counts(0:Nbins)
        double precision, intent(out) :: bin_avg(0:Nbins)
        double precision, intent(out) :: bin_centres(0:Nbins)
        
        integer                       :: i, bin_slot
        integer                       :: Ndata
        
        Ndata = size(x)
        
        counts = 0
        bin_avg = 0.d0
        bin_centres = 0.d0
        
        do i=1,Ndata
            bin_slot = floor(x(i))
            counts(bin_slot) = counts(bin_slot) + 1
            bin_avg(bin_slot) = bin_avg(bin_slot) + y(i)
            bin_centres(bin_slot) = bin_centres(bin_slot) + x(i)
        end do
        
        where (counts /= 0)
            bin_avg = bin_avg / dble(counts)
            bin_centres = bin_centres / dble(counts)
        end where
        
        return
    end subroutine d_square_binning
    
    subroutine d_square_binning1D(y, Nbins, counts, bin_avg, bin_centres)
        double precision, intent(in) :: y(:)
        integer, intent(in)           :: Nbins
        
        integer, intent(out)          :: counts(0:Nbins)
        double precision, intent(out) :: bin_avg(0:Nbins)
        double precision, intent(out) :: bin_centres(0:Nbins)
        
        integer                       :: i, bin_slot
        integer                       :: ylbound(1), yubound(1)
        
        ylbound = lbound(y)
        yubound = ubound(y)
        
        counts = 0
        bin_avg = 0.d0
        bin_centres = 0.d0
        
        do i=ylbound(1),yubound(1)
            bin_slot = abs(i)
            counts(bin_slot) = counts(bin_slot) + 1
            bin_avg(bin_slot) = bin_avg(bin_slot) + y(i)
            bin_centres(bin_slot) = bin_centres(bin_slot) + dble(bin_slot)
        end do
        
        where (counts /= 0)
            bin_avg = bin_avg / dble(counts)
            bin_centres = bin_centres / dble(counts)
        end where
        
        return
    end subroutine d_square_binning1D
    
    subroutine d_square_binning2D(y, Nbins, counts, bin_avg, bin_centres)
        double precision, intent(in)  :: y(:,:)
        integer, intent(in)           :: Nbins
        
        integer, intent(out)          :: counts(0:Nbins)
        double precision, intent(out) :: bin_avg(0:Nbins)
        double precision, intent(out) :: bin_centres(0:Nbins)
        
        integer                       :: i, j, bin_slot
        integer                       :: ylbound(2), yubound(2)
        double precision              :: x, jsqd
        
        ylbound = lbound(y)
        yubound = ubound(y)
        
        counts = 0
        bin_avg = 0.d0
        bin_centres = 0.d0
        
        do j=ylbound(2),yubound(2)
            jsqd = dble(j**2)
            do i=ylbound(1),yubound(1)
                x = sqrt(dble(i**2) + jsqd)
                bin_slot = floor(x)
                counts(bin_slot) = counts(bin_slot) + 1
                bin_avg(bin_slot) = bin_avg(bin_slot) + y(i,j)
                bin_centres(bin_slot) = bin_centres(bin_slot) + x
            end do
        end do
        
        where (counts /= 0)
            bin_avg = bin_avg / dble(counts)
            bin_centres = bin_centres / dble(counts)
        end where
        
        return
    end subroutine d_square_binning2D
    
    subroutine d_square_binning3D(y, Nbins, counts, bin_avg, bin_centres)
        double precision, intent(in)  :: y(:,:,:)
        integer, intent(in)           :: Nbins
        
        integer, intent(out)          :: counts(0:Nbins)
        double precision, intent(out) :: bin_avg(0:Nbins)
        double precision, intent(out) :: bin_centres(0:Nbins)
        
        integer                       :: i, j, k, bin_slot
        integer                       :: ylbound(3), yubound(3)
        double precision              :: x, jsqd, ksqd
        
        ylbound = lbound(y)
        yubound = ubound(y)
        
        counts = 0
        bin_avg = 0.d0
        bin_centres = 0.d0
        
        do k=ylbound(3),yubound(3)
            ksqd = dble(k**2)
            do j=ylbound(2),yubound(2)
                jsqd = dble(j**2)
                do i=ylbound(1),yubound(1)
                    x = sqrt(dble(i**2) + jsqd + ksqd)
                    bin_slot = floor(x)
                    counts(bin_slot) = counts(bin_slot) + 1
                    bin_avg(bin_slot) = bin_avg(bin_slot) + y(i,j,k)
                    bin_centres(bin_slot) = bin_centres(bin_slot) + x
                end do
            end do
        end do
        
        where (counts /= 0)
            bin_avg = bin_avg / dble(counts)
            bin_centres = bin_centres / dble(counts)
        end where
        
        return
    end subroutine d_square_binning3D
    
    subroutine d_triangle_binning(x, y, Nbins, bin_avg)
        double precision, intent(in)  :: x(:)
        double precision, intent(in)  :: y(:)
        integer, intent(in)           :: Nbins
        
        double precision, intent(out) :: bin_avg(0:Nbins)
        
        integer                       :: i, lower_bin, upper_bin
        integer                       :: Ndata
        double precision              :: count_frac(0:Nbins)
        double precision              :: lower_bin_frac, upper_bin_frac
        
        Ndata = size(x)
        
        count_frac = 0.d0
        bin_avg = 0.d0
        
        do i=1,Ndata
            lower_bin = floor(x(i))
            upper_bin = lower_bin + 1
            
            lower_bin_frac = upper_bin - x(i)
            upper_bin_frac = 1.d0 - lower_bin_frac
            
            count_frac(lower_bin) = count_frac(lower_bin) + lower_bin_frac
            count_frac(upper_bin) = count_frac(upper_bin) + upper_bin_frac
            
            bin_avg(lower_bin) = bin_avg(lower_bin) + y(i) * lower_bin_frac
            bin_avg(upper_bin) = bin_avg(upper_bin) + y(i) * upper_bin_frac
        end do
        
        where (count_frac > tiny(count_frac))
            bin_avg = bin_avg / count_frac
        end where
        
        return
    end subroutine d_triangle_binning
    
    subroutine d_triangle_binning1D(y, Nbins, bin_avg)
        double precision, intent(in)  :: y(:)
        integer, intent(in)           :: Nbins
        
        double precision, intent(out) :: bin_avg(0:Nbins)
        
        integer                       :: i, lower_bin, upper_bin
        integer                       :: ylbound(1), yubound(1)
        double precision              :: count_frac(0:Nbins)
        double precision              :: lower_bin_frac, upper_bin_frac
        
        ylbound = lbound(y)
        yubound = ubound(y)
        
        count_frac = 0.d0
        bin_avg = 0.d0
        
        do i=ylbound(1),yubound(1)
            lower_bin = abs(i)
            upper_bin = lower_bin + 1
            
            lower_bin_frac = 1.d0
            upper_bin_frac = 0.d0
            
            count_frac(lower_bin) = count_frac(lower_bin) + lower_bin_frac
            count_frac(upper_bin) = count_frac(upper_bin) + upper_bin_frac
            
            bin_avg(lower_bin) = bin_avg(lower_bin) + y(i) * lower_bin_frac
            bin_avg(upper_bin) = bin_avg(upper_bin) + y(i) * upper_bin_frac
        end do
        
        where (count_frac > tiny(count_frac))
            bin_avg = bin_avg / count_frac
        end where
        
        return
    end subroutine d_triangle_binning1D
    
    subroutine d_triangle_binning2D(y, Nbins, bin_avg)
        double precision, intent(in)  :: y(:,:)
        integer, intent(in)           :: Nbins
        
        double precision, intent(out) :: bin_avg(0:Nbins)
        
        integer                       :: i, j, lower_bin, upper_bin
        integer                       :: ylbound(2), yubound(2)
        double precision              :: count_frac(0:Nbins)
        double precision              :: x, jsqd
        double precision              :: lower_bin_frac, upper_bin_frac
        
        ylbound = lbound(y)
        yubound = ubound(y)
        
        count_frac = 0.d0
        bin_avg = 0.d0
        
        do j=ylbound(2),yubound(2)
            jsqd = dble(j**2)
            do i=ylbound(1),yubound(1)
                x = sqrt(dble(i**2) + jsqd)
                lower_bin = floor(x)
                upper_bin = lower_bin + 1
                
                lower_bin_frac = upper_bin - x
                upper_bin_frac = 1.d0 - lower_bin_frac
                
                count_frac(lower_bin) = count_frac(lower_bin) + lower_bin_frac
                count_frac(upper_bin) = count_frac(upper_bin) + upper_bin_frac
                
                bin_avg(lower_bin) = bin_avg(lower_bin) + &
                    y(i,j) * lower_bin_frac
                bin_avg(upper_bin) = bin_avg(upper_bin) + &
                    y(i,j) * upper_bin_frac
            end do
        end do
        
        where (count_frac > tiny(count_frac))
            bin_avg = bin_avg / count_frac
        end where
        
        return
    end subroutine d_triangle_binning2D
    
    subroutine d_triangle_binning3D(y, Nbins, bin_avg)
        double precision, intent(in)  :: y(:,:,:)
        integer, intent(in)           :: Nbins
        
        double precision, intent(out) :: bin_avg(0:Nbins)
        
        integer                       :: i, j, k, lower_bin, upper_bin
        integer                       :: ylbound(3), yubound(3)
        double precision              :: count_frac(0:Nbins)
        double precision              :: x, jsqd, ksqd
        double precision              :: lower_bin_frac, upper_bin_frac
        
        ylbound = lbound(y)
        yubound = ubound(y)
        
        count_frac = 0.d0
        bin_avg = 0.d0
        
        do k=ylbound(3),yubound(3)
            ksqd = dble(k**2)
            do j=ylbound(2),yubound(2)
                jsqd = dble(j**2)
                do i=ylbound(1),yubound(1)
                    x = sqrt(dble(i**2) + jsqd + ksqd)
                    lower_bin = floor(x)
                    upper_bin = lower_bin + 1
                    
                    lower_bin_frac = upper_bin - x
                    upper_bin_frac = 1.d0 - lower_bin_frac
                    
                    count_frac(lower_bin) = count_frac(lower_bin) + &
                        lower_bin_frac
                    count_frac(upper_bin) = count_frac(upper_bin) + &
                        upper_bin_frac
                    
                    bin_avg(lower_bin) = bin_avg(lower_bin) + &
                        y(i,j,k) * lower_bin_frac
                    bin_avg(upper_bin) = bin_avg(upper_bin) + &
                        y(i,j,k) * upper_bin_frac
                end do
            end do
        end do
        
        where (count_frac > tiny(count_frac))
            bin_avg = bin_avg / count_frac
        end where
        
        return
    end subroutine d_triangle_binning3D

end module binning
