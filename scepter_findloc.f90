!************************************************************************
! Name: scepter_findloc
! Purpose: Find location of a value in an array
!************************************************************************

module scepter_findloc
    use scepter_constants
    use scepter_variables
    implicit none

    private
    public :: findloc

    contains 
    
    !-----------------------------------------------------------------------
    ! Custom implementation of findloc for compatibility
    ! Searches for a specific string in an array of strings
    !-----------------------------------------------------------------------
#ifdef no_intr_findloc
        function findloc(chrlist_in,chrspecific,dim)
            implicit none
            character(*),intent(in)::chrlist_in(:),chrspecific
            integer,intent(in)::dim
            integer findloc,i

            findloc = 0
            do i=1, size(chrlist_in,dim=dim)
                if (trim(adjustl(chrspecific)) == trim(adjustl(chrlist_in(i)))) then
                    findloc = i
                    return
                endif 
            enddo 
        endfunction findloc
#endif 
    
endmodule scepter_findloc