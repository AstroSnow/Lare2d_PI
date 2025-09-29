!Module for the two-fluid version of Lare2D
!Experimental version. Not for public use.
!Written by Ben Snow (Exeter University) - 2025


MODULE two_fluid

  USE shared_data
  !USE boundary

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: setup_two_fluid

  REAL(num), DIMENSION(:,:), ALLOCATABLE :: visc_heat, rho_v, cv_v
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: alpha1, alpha2
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: pressure_n
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: qx,qy,qz
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: fx_visc,fy_visc,fz_visc
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: flux_x,flux_y,flux_z

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE setup_two_fluid
        !Reserved for any potential setup requirements
        write(*,*) 'setup for two fluid'

        ALLOCATE(xi_n(-1:nx+2, -1:ny+2))
        xi_n = 0.0_num
        
    END SUBROUTINE setup_two_fluid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
END MODULE two_fluid
