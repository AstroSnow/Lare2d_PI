  ! Copyright 2020 University of Warwick

  ! Licensed under the Apache License, Version 2.0 (the "License");
  ! you may not use this file except in compliance with the License.
  ! You may obtain a copy of the License at

  !    http://www.apache.org/licenses/LICENSE-2.0

  ! Unless required by applicable law or agreed to in writing, software
  ! distributed under the License is distributed on an "AS IS" BASIS,
  ! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ! See the License for the specific language governing permissions and
  ! limitations under the License.
  
MODULE remap

  USE shared_data,ONLY : num,rke,dxb,dyb,xpass,ypass,nx,ny,delta_ke,ix,iy
  !USE boundary
  !USE xremap
  !USE yremap
  !USE zremap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eulerian_remap

CONTAINS

  ! Strang splitting

  SUBROUTINE eulerian_remap(i,bx,by,bz,cv,cv1)

    REAL(num),intent(inout):: bx(-1:nx+2,-1:ny+2),by(-1:nx+2,-1:ny+2),bz(-1:nx+2,-1:ny+2)
    REAL(num),intent(inout):: cv(-1:nx+2,-1:ny+2),cv1(-1:nx+2,-1:ny+2)
    INTEGER, INTENT(IN) :: i
    INTEGER :: case_test

    IF (rke) delta_ke = 0.0_num
    xpass = 1
    ypass = 1

    DO ix = -1, nx + 2
      DO iy = -1, ny + 2
        bx(ix,iy) = bx(ix,iy) * dyb(iy)
        by(ix,iy) = by(ix,iy) * dxb(ix)
        bz(ix,iy) = bz(ix,iy) * cv(ix,iy)
      END DO
    END DO

    DO iy = -1, ny + 2
      bx(-2,iy) = bx(-2,iy) * dyb(iy)
    END DO

    DO ix = -1, nx + 2
      by(ix,-2) = by(ix,-2) * dxb(ix)
    END DO

    case_test = MODULO(i, 6)

    ! Strang ordering
    SELECT CASE(case_test)
    CASE (0)
      CALL remap_x(cv1)
      CALL remap_y(cv1)
      CALL remap_z(cv1)
    CASE (1)
      CALL remap_y(cv1)
      CALL remap_z(cv1)
      CALL remap_x(cv1)
    CASE (2)
      CALL remap_z(cv1)
      CALL remap_x(cv1)
      CALL remap_y(cv1)
    CASE (3)
      CALL remap_x(cv1)
      CALL remap_z(cv1)
      CALL remap_y(cv1)
    CASE (4)
      CALL remap_z(cv1)
      CALL remap_y(cv1)
      CALL remap_x(cv1)
    CASE (5)
      CALL remap_y(cv1)
      CALL remap_x(cv1)
      CALL remap_z(cv1)
    END SELECT

    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        bx(ix,iy) = bx(ix,iy) / dyb(iy)
        by(ix,iy) = by(ix,iy) / dxb(ix)
        bz(ix,iy) = bz(ix,iy) / cv(ix,iy)
      END DO
    END DO

    DO iy = -1, ny + 2
      bx(-2,iy) = bx(-2,iy) / dyb(iy)
    END DO

    DO ix = -1, nx + 2
      by(ix,-2) = by(ix,-2) / dxb(ix)
    END DO

    CALL bfield_bcs

  END SUBROUTINE eulerian_remap

END MODULE remap
