!Module for the two-fluid version of Lare2D
!Experimental version. Not for public use.
!Written by Ben Snow (Exeter University) - 2025


MODULE two_fluid

  USE shared_data
  !USE boundary

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: setup_two_fluid


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
