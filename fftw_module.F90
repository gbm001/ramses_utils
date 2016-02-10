module fftw_module
    use amr_utils
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    contains
    
    subroutine DFT_1D(real_field, complex_field, N)
       ! Transform purely real field into complex field

       integer, intent(in)            :: N        ! Size of grid
       integer(kind=C_INT)            :: N_c      ! As C integer
       real(kind=dp), intent(in)      :: real_field(0:N-1)
                                               ! Scalar field
       complex(kind=C_DOUBLE_COMPLEX), intent(out) :: complex_field(0:N-1)
                                               ! Complex field output
       
       type(C_PTR)                    :: plan  ! FFTW plan pointer

       N_c = N
       
       complex_field = cmplx(real_field, kind=C_DOUBLE_COMPLEX)

       plan = fftw_plan_dft_1d(N_c, complex_field, complex_field,&
                              &FFTW_FORWARD, FFTW_ESTIMATE)

       !write (6,*) "Performing FFT..."

       call fftw_execute_dft(plan, complex_field, complex_field)
       complex_field = complex_field / sqrt(real(N, dp))
       
       call fftw_destroy_plan(plan)
       
       return
    end subroutine DFT_1D

    subroutine DFT_2D(real_field, complex_field, N)
       ! Transform purely real field into complex field

       integer, intent(in)             :: N        ! Size of grid
       integer(kind=C_INT)             :: N_c      ! As C integer
       real(kind=dp), intent(in)       :: real_field(0:N-1, 0:N-1)
                                                   ! Scalar field
       complex(kind=C_DOUBLE_COMPLEX), intent(out) :: complex_field(0:N-1, 0:N-1)
                                                   ! Complex field output
       
       type(C_PTR)                    :: plan  ! FFTW plan pointer

       N_c = N
       
       complex_field = cmplx(real_field, kind=C_DOUBLE_COMPLEX)

       plan = fftw_plan_dft_2d(N_c, N_c, complex_field, complex_field,&
                              &FFTW_FORWARD, FFTW_ESTIMATE)

       !write (6,*) "Performing FFT..."
       call fftw_execute_dft(plan, complex_field, complex_field)
       complex_field = complex_field / sqrt(real(N**2, dp))
   
       call fftw_destroy_plan(plan)

       return
    end subroutine DFT_2D

    subroutine DFT_3D(real_field, complex_field, N)
       ! Transform purely real field into complex field

       integer, intent(in)            :: N     ! Size of grid
       integer(kind=C_INT)            :: N_c   ! As C integer
       real(kind=dp), intent(in)      :: real_field(0:N-1, 0:N-1, 0:N-1)
                                               ! Scalar field
       complex(kind=C_DOUBLE_COMPLEX), intent(out) :: complex_field(0:N-1, 0:N-1, 0:N-1)
                                               ! Complex field output
       
       type(C_PTR)                    :: plan  ! FFTW plan pointer

       N_c = N
       
       complex_field = cmplx(real_field, kind=C_DOUBLE_COMPLEX)

       plan = fftw_plan_dft_3d(N_c, N_c, N_c, complex_field, complex_field,&
                              &FFTW_FORWARD, FFTW_ESTIMATE)

       !write (6,*) "Performing FFT..."
       call fftw_execute_dft(plan, complex_field, complex_field)
       complex_field = complex_field / sqrt(real(N**3, dp))
   
       call fftw_destroy_plan(plan)

       return
    end subroutine DFT_3D

    subroutine FFT_1D(complex_field_in, complex_field_out, N)
       integer, intent(in)            :: N        ! Size of grid
       complex(kind=C_DOUBLE_COMPLEX), intent(inout) :: complex_field_in(0:N-1)
                                                  ! Complex field input
       complex(kind=C_DOUBLE_COMPLEX), intent(out)   :: complex_field_out(0:N-1)
                                                  ! Complex field output
       call FFT_1D_int(complex_field_in, complex_field_out, N, .TRUE.)
    end subroutine FFT_1D

    subroutine IFT_1D(complex_field_in, complex_field_out, N)
       integer, intent(in)            :: N        ! Size of grid
       complex(kind=C_DOUBLE_COMPLEX), intent(inout) :: complex_field_in(0:N-1)
                                                  ! Complex field input
       complex(kind=C_DOUBLE_COMPLEX), intent(out)   :: complex_field_out(0:N-1)
                                                  ! Complex field output
       call FFT_1D_int(complex_field_in, complex_field_out, N, .FALSE.)
    end subroutine IFT_1D

    subroutine FFT_1D_int(complex_field_in, complex_field_out, N, forwards)
       ! Transform complex field into complex field

       integer, intent(in)            :: N        ! Size of grid
       complex(kind=C_DOUBLE_COMPLEX), intent(inout) :: complex_field_in(0:N-1)
                                                  ! Complex field input
       complex(kind=C_DOUBLE_COMPLEX), intent(out)   :: complex_field_out(0:N-1)
                                                  ! Complex field output
       logical, intent(in)            :: forwards ! Forwards transform?
       
       integer(kind=C_INT)            :: N_c      ! As C integer
       type(C_PTR)                    :: plan  ! FFTW plan pointer

       N_c = N

       if (forwards) then
           plan = fftw_plan_dft_1d(N_c, complex_field_in, complex_field_out,&
                                  &FFTW_FORWARD, FFTW_ESTIMATE)
       else
           plan = fftw_plan_dft_1d(N_c, complex_field_in, complex_field_out,&
                                  &FFTW_BACKWARD, FFTW_ESTIMATE)
       end if

       !write (6,*) "Performing FFT..."
       call fftw_execute_dft(plan, complex_field_in, complex_field_out)
       complex_field_out = complex_field_out / sqrt(real(N, dp))

       call fftw_destroy_plan(plan)

       return
    end subroutine FFT_1D_int

    subroutine FFT_2D(complex_field_in, complex_field_out, N)
       integer, intent(in)            :: N        ! Size of grid
       complex(kind=C_DOUBLE_COMPLEX), intent(inout) :: &
           &complex_field_in(0:N-1, 0:N-1)        ! Complex field input
       complex(kind=C_DOUBLE_COMPLEX), intent(out)   :: &
           &complex_field_out(0:N-1, 0:N-1)       ! Complex field output
       call FFT_2D_int(complex_field_in, complex_field_out, N, .TRUE.)
    end subroutine FFT_2D

    subroutine IFT_2D(complex_field_in, complex_field_out, N)
       integer, intent(in)            :: N        ! Size of grid
       complex(kind=C_DOUBLE_COMPLEX), intent(inout) :: &
           &complex_field_in(0:N-1, 0:N-1)        ! Complex field input
       complex(kind=C_DOUBLE_COMPLEX), intent(out)   :: &
           &complex_field_out(0:N-1, 0:N-1)       ! Complex field output
       call FFT_2D_int(complex_field_in, complex_field_out, N, .FALSE.)
    end subroutine IFT_2D

    subroutine FFT_2D_int(complex_field_in, complex_field_out, N, forwards)
       ! Transform complex field into complex field

       integer, intent(in)            :: N        ! Size of grid
       complex(kind=C_DOUBLE_COMPLEX), intent(inout) :: &
           &complex_field_in(0:N-1, 0:N-1)        ! Complex field input
       complex(kind=C_DOUBLE_COMPLEX), intent(out)   :: &
           &complex_field_out(0:N-1, 0:N-1)       ! Complex field output
       logical, intent(in)            :: forwards ! Forwards transform?
       
       integer(kind=C_INT)            :: N_c      ! As C integer
       type(C_PTR)                    :: plan  ! FFTW plan pointer

       N_c = N

       if (forwards) then
           plan = fftw_plan_dft_2d(N_c, N_c, complex_field_in, complex_field_out,&
                                  &FFTW_FORWARD, FFTW_ESTIMATE)
       else
           plan = fftw_plan_dft_2d(N_c, N_c, complex_field_in, complex_field_out,&
                                  &FFTW_BACKWARD, FFTW_ESTIMATE)
       end if

       !write (6,*) "Performing FFT..."
       call fftw_execute_dft(plan, complex_field_in, complex_field_out)
       complex_field_out = complex_field_out / sqrt(real(N**2, dp))

       call fftw_destroy_plan(plan)

       return
    end subroutine FFT_2D_int

    subroutine FFT_3D(complex_field_in, complex_field_out, N)
       integer, intent(in)            :: N        ! Size of grid
       complex(kind=C_DOUBLE_COMPLEX), intent(inout) :: &
           &complex_field_in(0:N-1, 0:N-1)        ! Complex field input
       complex(kind=C_DOUBLE_COMPLEX), intent(out)   :: &
           &complex_field_out(0:N-1, 0:N-1)       ! Complex field output
       call FFT_3D_int(complex_field_in, complex_field_out, N, .TRUE.)
    end subroutine FFT_3D

    subroutine IFT_3D(complex_field_in, complex_field_out, N)
       integer, intent(in)            :: N          ! Size of grid
       complex(kind=C_DOUBLE_COMPLEX), intent(inout) :: &
           &complex_field_in(0:N-1, 0:N-1, 0:N-1)   ! Complex field input
       complex(kind=C_DOUBLE_COMPLEX), intent(out)   :: &
           &complex_field_out(0:N-1, 0:N-1, 0:N-1)  ! Complex field output
       call FFT_3D_int(complex_field_in, complex_field_out, N, .FALSE.)
    end subroutine IFT_3D

    subroutine FFT_3D_int(complex_field_in, complex_field_out, N, forwards)
       ! Transform complex field into complex field

       integer, intent(in)            :: N          ! Size of grid
       complex(kind=C_DOUBLE_COMPLEX), intent(inout) :: &
           &complex_field_in(0:N-1, 0:N-1, 0:N-1)   ! Complex field input
       complex(kind=C_DOUBLE_COMPLEX), intent(out)   :: &
           &complex_field_out(0:N-1, 0:N-1, 0:N-1)  ! Complex field output
       logical, intent(in)            :: forwards ! Forwards transform?
       
       integer(kind=C_INT)            :: N_c        ! As C integer
       type(C_PTR)                    :: plan  ! FFTW plan pointer

       N_c = N

       if (forwards) then
           plan = fftw_plan_dft_3d(N_c, N_c, N_c, complex_field_in,&
                                  &complex_field_out, FFTW_FORWARD,&
                                  &FFTW_ESTIMATE)
       else
           plan = fftw_plan_dft_3d(N_c, N_c, N_c, complex_field_in,&
                                  &complex_field_out, FFTW_BACKWARD,&
                                  &FFTW_ESTIMATE)
       end if

       !write (6,*) "Performing FFT..."
       call fftw_execute_dft(plan, complex_field_in, complex_field_out)
       complex_field_out = complex_field_out / sqrt(real(N**3, dp))

       call fftw_destroy_plan(plan)

       return
    end subroutine FFT_3D_int

end module fftw_module