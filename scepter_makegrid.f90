!************************************************************************
! Name: scepter_makegrid
! Purpose: Create computational grid
!************************************************************************

module scepter_makegrid
    use scepter_constants
    use scepter_variables
    implicit none

    private
    public :: make_grid

    contains
    !-----------------------------------------------------------------------
    ! Subroutine: make_grid
    ! Purpose: create computational grid, after Hoffmann & Chiang, 2000
    ! Creates either regular or non-uniform grid based on beta parameter
    ! For non-uniform grid, uses transformation to concentrate points near surface     
    !-----------------------------------------------------------------------
    subroutine makegrid(beta,nz,ztot,dz,z,regular_grid)  
        implicit none
        integer(kind=4),intent(in) :: nz          ! Number of grid points
        logical,intent(in)::regular_grid          ! Flag for regular vs. non-uniform grid
        real(kind=8),intent(in)::beta,ztot       ! Grid transformation parameter and total depth
        real(kind=8),intent(out)::dz(nz),z(nz)   ! Grid spacing and node positions
        integer(kind=4) iz

        do iz = 1, nz 
            z(iz) = iz*ztot/nz  ! regular grid 
            if (iz==1) then
                dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0))
            endif
            if (iz/=1) then 
                dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0)) - sum(dz(:iz-1))
            endif
        enddo

        ! Override with regular grid if specified
        if (regular_grid) then 
            dz = ztot/nz  ! when implementing regular grid
        endif 

        do iz=1,nz  ! depth is defined at the middle of individual layers 
            if (iz==1) z(iz)=dz(iz)*0.5d0  
            if (iz/=1) z(iz) = z(iz-1)+dz(iz-1)*0.5d0 + 0.5d0*dz(iz)
        enddo

    endsubroutine makegrid

endmodule scepter_makegrid