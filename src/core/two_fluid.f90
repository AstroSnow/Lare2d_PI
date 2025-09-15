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

    END SUBROUTINE setup_two_fluid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE two_fluid
