module quickselect_mod
    implicit none
    
    contains
    
    function partition(N, array, pivot)
        ! MODIFIES THE ARRAY array
        ! Take an array and separate its elements into those
        ! lower or higher than the pivot. Return the index of the pivot.
        integer, intent(in)             :: N
        double precision, intent(inout) :: array(1:N)
        integer, intent(in)             :: pivot
        integer                         :: partition
        
        integer                         :: i
        double precision                :: temp, pivot_value
        
        ! Save pivot value, copy last value into pivot position
        pivot_value = array(pivot)
        array(pivot) = array(N)
        
        ! Iterate over each element
        partition = 1
        do i=1,N-1
            if (array(i) < pivot_value) then
                temp = array(i)
                array(i) = array(partition)
                array(partition) = temp
                partition = partition + 1
            end if
        end do
        
        ! Finally reinstate pivot
        array(N) = array(partition)
        array(partition) = pivot_value
        
        ! The pivot is now in its correct place at index 'partition'
        ! All elements are correctly divided above or below it
        
    end function partition
    
    function quickselect(array, k)
        ! MODIFIES THE ARRAY array
        ! Take an array and returns the k'th largest element
        double precision, intent(inout) :: array(:)
        integer, intent(in)             :: k
        double precision                :: quickselect
        
        integer                         :: left, right, N, pivot
        
        ! Basic validation
        N = size(array)
        if (k > N) stop "k > number of elements in array!"
        if (N == 1) then
            quickselect = array(1)
            return
        end if
        
        ! Starting positions: whole array
        left = 1
        right = N
        
        do
            ! Set a pivot
            pivot = (left + right) / 2       ! Simplest option
            ! Partition the array around that pivot
            pivot = partition(right-left+1, array(left:right), pivot-left+1)
            pivot = pivot + left-1
            ! Compare our (now correctly placed) pivot with k
            if (pivot==k) then
                ! We have our desired element
                quickselect = array(k)
                exit
            else if (k < pivot) then
                ! Our desired element is in the left-hand partition
                right = pivot - 1
            else
                ! Our desired element is in the right-hand partition
                left = pivot + 1
            end if
        end do
        
        continue
    end function quickselect
    
    function median(array)
        ! MODIFIES THE ARRAY array
        ! Finds the median of the array
        ! Calls quickselect twice for even-sized arrays rather
        ! than a more intelligent strategy
        double precision, intent(inout) :: array(:)
        double precision                :: median
        
        integer                         :: N
        double precision                :: aux1, aux2
        
        N = size(array)
        
        if (N==0) stop 'Array must be longer than 0!'
        
        if (mod(N,2) == 0) then
            ! Even number, run quickselect twice
            aux1 = quickselect(array, (N/2))
            aux2 = quickselect(array, (N/2)+1)
            median = (aux1 + aux2) / 2.d0
        else
            ! Odd number, run quickselect once
            median = quickselect(array, (N/2)+1)
        end if
        
    end function median

end module quickselect_mod