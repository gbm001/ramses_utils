program convert_to_single_program
    use amr_utils
    
    character (LEN=:), allocatable :: output_basedir
    character (LEN=:), allocatable :: ioutput_str
    integer                        :: ioutput
    integer                        :: length
    
    lowmem_tables = .TRUE.
    
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
    
    call convert_to_single(output_basedir, ioutput)
    
end program convert_to_single_program
