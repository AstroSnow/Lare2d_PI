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
  
!******************************************************************************
! Lagrangian step routines
!******************************************************************************

MODULE lagran

  USE shared_data, only : num,nx,ny,cv1_plasma,sixth,cowling_resistivity,none_zero,boris,comm,errcode,&
            dt,dt_factor,dt_from_restart,dt_multiplier,dt_previous,dt_snapshots,gamma,hall_mhd,&
            largest_number,MPI_MIN,mpireal,restart,time,predictor_step,step,va_max2,any_open,&
            visc1,visc2,ionise_pot,dxc,dyc
  !USE boundary, only : energy_bcs
  !USE neutral
  !USE conduct
  !USE radiative
  !USE openboundary
  !USE remap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lagrangian_step, eta_calc

  ! Only used inside lagran.f90
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: alpha1, alpha2
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: visc_heat, pressure, rho_v, cv_v
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: flux_x, flux_y, flux_z, curlb
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: fx_visc, fy_visc, fz_visc
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: bx0, by0, bz0
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: qx, qy, qz
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: energy0, delta_energy
  
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: rho_temp, en_temp
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: vx_temp,vy_temp,vz_temp
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: bx_temp,by_temp,bz_temp
  
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: bx1, by1, bz1
  
  REAL(num)::dt2 !This is the timestep/2 for the predictor step
  
  !I am not sure what these variables are but they don't seem to be globally passed.
  REAL(num)::w1,w2
  
  INTEGER :: ix,iy,ixm,iym,ixp,iyp

CONTAINS

  !****************************************************************************
  ! This subroutine manages the progress of the lagrangian step
  !****************************************************************************

  SUBROUTINE lagrangian_step(rho,vx,vy,vz,energy,bx,by,bz)

    INTEGER :: substeps, subcycle
    REAL(num) :: actual_dt, dt_sub
    REAL(num),intent(inout):: rho(-1:nx+1,-1:ny+1)
    REAL(num),intent(inout):: vx(-2:nx+2, -2:ny+2),vy(-2:nx+2, -2:ny+2),vz(-2:nx+2, -2:ny+2)
    REAL(num),intent(inout):: energy(-1:nx+2,-1:ny+2)
    REAL(num),intent(inout):: bx(-1:nx+2,-1:ny+2),by(-1:nx+2,-1:ny+2),bz(-1:nx+2,-1:ny+2)

#ifndef CAUCHY
    ALLOCATE(bx0(-2:nx+2,-1:ny+2))
    ALLOCATE(by0(-1:nx+2,-2:ny+2))
    ALLOCATE(bz0(-1:nx+2,-1:ny+2))
#endif
    ALLOCATE(bx1(-1:nx+2,-1:ny+2))
    ALLOCATE(by1(-1:nx+2,-1:ny+2))
    ALLOCATE(bz1(-1:nx+2,-1:ny+2))
    ALLOCATE(alpha1(0:nx+1,0:ny+2))
    ALLOCATE(alpha2(-1:nx+1,0:ny+1))
    ALLOCATE(visc_heat(0:nx+1,0:ny+1))
    ALLOCATE(pressure(-1:nx+2,-1:ny+2))
    ALLOCATE(rho_v(-1:nx+1,-1:ny+1))
    ALLOCATE(cv_v(-1:nx+1,-1:ny+1))
    ALLOCATE(qx(0:nx+1,0:ny+1))
    ALLOCATE(qy(0:nx+1,0:ny+1))
    ALLOCATE(qz(0:nx+1,0:ny+1))
    ALLOCATE(fx_visc(0:nx,0:ny))
    ALLOCATE(fy_visc(0:nx,0:ny))
    ALLOCATE(fz_visc(0:nx,0:ny))
    ALLOCATE(flux_x(0:nx,0:ny))
    ALLOCATE(flux_y(0:nx,0:ny))
    ALLOCATE(flux_z(0:nx,0:ny))
    ALLOCATE(curlb (0:nx,0:ny))
    ALLOCATE(energy0(-1:nx+2,-1:ny+2))
    ALLOCATE(delta_energy(-1:nx+2,-1:ny+2))

    DO iy = -1, ny + 2
      iym = iy - 1
      DO ix = -1, nx + 2
        ixm = ix - 1
        bx1(ix,iy) = (bx(ix,iy) + bx(ixm,iy )) * 0.5_num
        by1(ix,iy) = (by(ix,iy) + by(ix ,iym)) * 0.5_num
        bz1(ix,iy) = bz(ix,iy)

        pressure(ix,iy) = (gamma - 1.0_num) * rho(ix,iy) &
            * (energy(ix,iy) - (1.0_num - xi_n(ix,iy)) * ionise_pot)
      END DO
    END DO

    DO iy = -1, ny + 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ixp = ix + 1
        rho_v(ix,iy) = rho(ix,iy) * cv(ix,iy) + rho(ixp,iy) * cv(ixp,iy) &
            +   rho(ix,iyp) * cv(ix,iyp) + rho(ixp,iyp) * cv(ixp,iyp)
        cv_v(ix,iy) = cv(ix,iy) + cv(ixp,iy) + cv(ix,iyp) + cv(ixp,iyp)
        rho_v(ix,iy) = rho_v(ix,iy) / cv_v(ix,iy)

        cv_v(ix,iy) = 0.25_num * cv_v(ix,iy)
      END DO
    END DO

    CALL edge_shock_viscosity
    CALL set_dt
    dt2 = dt * 0.5_num

    energy(:,:) = MAX(energy(:,:), 0.0_num)

    CALL predictor_corrector_step(energy,bx,by,bz,cv1_plasma)

    DEALLOCATE(bx1, by1, bz1, alpha1, alpha2)
    DEALLOCATE(visc_heat, pressure, rho_v, cv_v, flux_x, flux_y, flux_z, curlb)
    DEALLOCATE(qx, qy, qz)
    DEALLOCATE(fx_visc, fy_visc, fz_visc)
    DEALLOCATE(energy0, delta_energy)
#ifndef CAUCHY
    DEALLOCATE(bx0, by0, bz0)
#endif

    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

  END SUBROUTINE lagrangian_step



  !****************************************************************************
  ! The main predictor / corrector step which advances the momentum equation
  !****************************************************************************

  SUBROUTINE predictor_corrector_step(energy,bx,by,bz,cv1)

    REAL(num),intent(inout):: energy(-1:nx+2,-1:ny+2)
    REAL(num),intent(inout):: bx(-1:nx+2,-1:ny+2),by(-1:nx+2,-1:ny+2),bz(-1:nx+2,-1:ny+2)
    REAL(num),intent(inout)::cv1(-1:nx+2, -1:ny+2)
    REAL(num) :: pp, ppx, ppy, ppxy
    REAL(num) :: e1
    REAL(num) :: vxb, vxbm, vyb, vybm
    REAL(num) :: bxv, byv, bzv, jx, jy, jz
    REAL(num) :: dv
    REAL(num) :: fx, fy, fz
#ifdef CAUCHY
    REAL(num) :: cvx, cvxp, cvy, cvyp
#endif

#ifndef CAUCHY
    bx0(:,:) = bx(:,:) 
    by0(:,:) = by(:,:) 
    bz0(:,:) = bz(:,:) 
#endif

    CALL b_field_and_cv1_update

    bx1(:,:) = bx1(:,:) * cv1(:,:)
    by1(:,:) = by1(:,:) * cv1(:,:)
    bz1(:,:) = bz1(:,:) * cv1(:,:)

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        dv = cv1(ix,iy) / cv(ix,iy) - 1.0_num
        ! Predictor energy
        e1 = energy(ix,iy) - pressure(ix,iy) * dv / rho(ix,iy)
        e1 = e1 + visc_heat(ix,iy) * dt2 / rho(ix,iy)

        ! Now define the predictor step pressures
        pressure(ix,iy) = (e1 - (1.0_num - xi_n(ix,iy)) * ionise_pot) &
            * (gamma - 1.0_num) * rho(ix,iy) * cv(ix,iy) / cv1(ix,iy)
      END DO
    END DO

    DO iy = 0, ny
      iyp = iy + 1
      iym = iy - 1
      DO ix = 0, nx
        ixp = ix + 1
        ixm = ix - 1

        pp    = pressure(ix ,iy )
        ppx   = pressure(ixp,iy )
        ppy   = pressure(ix ,iyp)
        ppxy  = pressure(ixp,iyp)

        ! P total at Ex(i,j)
        w1 = (pp + ppy) * 0.5_num
        ! P total at Ex(i+1,j)
        w2 = (ppx + ppxy) * 0.5_num
        fx = -(w2 - w1) / dxc(ix)

        ! P total at Ey(i,j)
        w1 = (pp + ppx) * 0.5_num
        ! P total at Ey(i,j+1)
        w2 = (ppy + ppxy) * 0.5_num
        fy = -(w2 - w1) / dyc(iy)

        fz = 0.0_num

#ifdef CAUCHY
        cvx  = cv1(ix ,iy ) + cv1(ix ,iyp)
        cvxp = cv1(ixp,iy ) + cv1(ixp,iyp)
        cvy  = cv1(ix ,iy ) + cv1(ixp,iy )
        cvyp = cv1(ix ,iyp) + cv1(ixp,iyp)

        w1 = (bz1(ix ,iy ) + bz1(ixp,iy )) / cvy
        w2 = (bz1(ix ,iyp) + bz1(ixp,iyp)) / cvyp
        jx = (w2 - w1) / dyc(iy)

        w1 = (bz1(ix ,iy ) + bz1(ix ,iyp)) / cvx
        w2 = (bz1(ixp,iy ) + bz1(ixp,iyp)) / cvxp
        jy = -(w2 - w1) / dxc(ix)

        w1 = (by1(ix ,iy ) + by1(ix ,iyp)) / cvx
        w2 = (by1(ixp,iy ) + by1(ixp,iyp)) / cvxp
        jz = (w2 - w1) / dxc(ix)

        w1 = (bx1(ix ,iy ) + bx1(ixp,iy )) / cvy
        w2 = (bx1(ix ,iyp) + bx1(ixp,iyp)) / cvyp
        jz = jz - (w2 - w1) / dyc(iy)

        bxv = (bx1(ix,iy ) + bx1(ixp,iy ) + bx1(ix,iyp) + bx1(ixp,iyp)) &
            / (cvx + cvxp)

        byv = (by1(ix,iy ) + by1(ixp,iy ) + by1(ix,iyp) + by1(ixp,iyp)) &
            / (cvx + cvxp)

        bzv = (bz1(ix,iy ) + bz1(ixp,iy ) + bz1(ix,iyp) + bz1(ixp,iyp)) &
            / (cvx + cvxp)
#else
        bxv = 0.5_num * (bx(ix,iy) + bx(ix,iyp))   
        byv = 0.5_num * (by(ix,iy) + by(ixp,iy))
        bzv = 0.25_num * (bz(ix,iy ) + bz(ixp,iy ) + bz(ix,iyp) + bz(ixp,iyp)) 
        
        jx = 0.5_num * (bz(ix,iyp) + bz(ixp,iyp) - bz(ix,iy) - bz(ixp,iy)) / dyc(iy)
        jy = - 0.5_num * (bz(ixp,iy) + bz(ixp,iyp)- bz(ix,iy) - bz(ix,iyp)) / dxc(ix)
        jz = (by(ixp,iy) - by(ix,iy)) / dxc(ix) - (bx(ix,iyp) - bx(ix,iy)) / dyc(iy)
#endif
        fx = fx + gamma_boris(ix,iy) * (jy * bzv - jz * byv)
        fy = fy + gamma_boris(ix,iy) * (jz * bxv - jx * bzv)
        fz = fz + gamma_boris(ix,iy) * (jx * byv - jy * bxv)

        fy = fy - rho_v(ix,iy) * grav(iy)

        ! Find half step velocity needed for remap
        vx1(ix,iy) = vx(ix,iy) + dt2 * (fx_visc(ix,iy) + fx) / rho_v(ix,iy)
        vy1(ix,iy) = vy(ix,iy) + dt2 * (fy_visc(ix,iy) + fy) / rho_v(ix,iy)
        vz1(ix,iy) = vz(ix,iy) + dt2 * (fz_visc(ix,iy) + fz) / rho_v(ix,iy)
      END DO
    END DO

#ifndef CAUCHY
    bx(:,:) = bx0(:,:) 
    by(:,:) = by0(:,:) 
    bz(:,:) = bz0(:,:) 
#endif
    
    CALL remap_v_bcs

    CALL edge_shock_heating

    DO iy = 0, ny
      DO ix = 0, nx
        ! Velocity at the end of the Lagrangian step
        vx(ix,iy) = 2.0_num * vx1(ix,iy) - vx(ix,iy) 
        vy(ix,iy) = 2.0_num * vy1(ix,iy) - vy(ix,iy) 
        vz(ix,iy) = 2.0_num * vz1(ix,iy) - vz(ix,iy) 
      END DO
    END DO

    CALL velocity_bcs
    IF (any_open) THEN
      CALL open_bcs         
    END IF 

    ! Finally correct density and energy to final values
    DO iy = 1, ny
      iym = iy - 1
      DO ix = 1, nx
        ixm = ix - 1

        ! vx1 at Bx(i,j)
        vxb  = (vx1(ix ,iy ) + vx1(ix ,iym)) * 0.5_num
        ! vx1 at Bx(i-1,j)
        vxbm = (vx1(ixm,iy ) + vx1(ixm,iym)) * 0.5_num
        ! vy1 at By(i,j)
        vyb  = (vy1(ix ,iy ) + vy1(ixm,iy )) * 0.5_num
        ! vy1 at By(i,j-1)
        vybm = (vy1(ix ,iym) + vy1(ixm,iym)) * 0.5_num

        dv = ((vxb - vxbm) / dxb(ix) + (vyb - vybm) / dyb(iy)) * dt

        cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)

        ! Energy at end of Lagrangian step
        energy(ix,iy) = energy(ix,iy) &
            + (dt * visc_heat(ix,iy) - dv * pressure(ix,iy)) &
            / rho(ix,iy)

        visc_dep(ix,iy) = visc_dep(ix,iy) + dt * visc_heat(ix,iy)

        rho(ix,iy) = rho(ix,iy) / (1.0_num + dv)

        !total_visc_heating = total_visc_heating &
        !    + dt * visc_heat(ix,iy) * cv(ix,iy)

      END DO
    END DO

  END SUBROUTINE predictor_corrector_step



  !****************************************************************************
  ! This subroutine calculates the viscous effects and updates the
  ! magnetic field
  !****************************************************************************

  SUBROUTINE edge_shock_viscosity

    REAL(num) :: dvdots, dx, dxm, dxp
    REAL(num) :: b2, rmin
    REAL(num) :: a1, a2, a3, a4
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: cs, cs_v
    REAL(num),DIMENSION(:,:),ALLOCATABLE :: p_visc
    INTEGER :: i0, i1, i2, i3, j0, j1, j2, j3
    LOGICAL, SAVE :: first_call = .TRUE.
    
    REAL(num)::visc2_norm

    ALLOCATE(p_visc(-1:nx+2, -1:ny+2))
    ALLOCATE(cs(-1:nx+2,-1:ny+2), cs_v(-1:nx+1,-1:ny+1))

    IF (first_call) THEN
      first_call = .FALSE.
      visc2_norm = 0.25_num * (gamma + 1.0_num) * visc2
    END IF

    p_visc = 0.0_num
    visc_heat = 0.0_num

    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        rmin = MAX(rho(ix,iy), none_zero)
        b2 = bx1(ix,iy)**2 + by1(ix,iy)**2 + bz1(ix,iy)**2
        cs(ix,iy) = SQRT((gamma * pressure(ix,iy) + gamma_boris(ix,iy) * b2) / rmin)
      END DO
    END DO

    DO iy = -1, ny + 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ixp = ix + 1
        cs_v(ix,iy) = cs(ix,iy) * cv(ix,iy) + cs(ixp,iy) * cv(ixp,iy) &
            +   cs(ix,iyp) * cv(ix,iyp) + cs(ixp,iyp) * cv(ixp,iyp)
        cs_v(ix,iy) = 0.25_num * cs_v(ix,iy) / cv_v(ix,iy)
      END DO
    END DO

    DO iy = 0, ny + 2
      iym = iy - 1
      iyp = iy + 1
      DO ix = 0, nx + 1
        ixm = ix - 1
        ixp = ix + 1

        ! Edge viscosities from Caramana
        ! Triangles numbered as in Goffrey thesis

        ! Edge viscosity for triangle 1
        i1 = ixm
        j1 = iym
        i2 = ix
        j2 = iym
        i0 = i1 - 1
        j0 = j1
        i3 = i2 + 1
        j3 = j2
        dx = dxb(ix)
        dxp = dxb(ixp)
        dxm = dxb(ixm)
        ! dv in direction of dS, i.e. dv.dS / abs(dS)
        dvdots = - (vx(i1,j1) - vx(i2,j2))
        ! Force on node is alpha*dv*ds but store only alpha and convert to force
        ! when needed.  
        alpha1(ix,iy) = edge_viscosity()
      END DO
    END DO

    ! Edge viscosity for triangle 2
    DO iy = 0, ny + 1
      iym = iy - 1
      iyp = iy + 1
      DO ix = -1, nx + 1
        ixm = ix - 1
        ixp = ix + 1

        i1 = ix
        j1 = iym
        i2 = ix
        j2 = iy
        i0 = i1
        j0 = j1 - 1
        i3 = i2
        j3 = j2 + 1
        dx = dyb(iy)
        dxp = dyb(iyp)
        dxm = dyb(iym)
        dvdots = - (vy(i1,j1) - vy(i2,j2))
        alpha2(ix,iy) = edge_viscosity()
      END DO
    END DO

    DO iy = 0, ny + 1 
      iym = iy - 1
      iyp = iy + 1
      DO ix = 0, nx + 1 
        ixm = ix - 1
        ixp = ix + 1
        ! Estimate p_visc based on alpha * dv, for timestep control
        a1 = ((vx(ixm,iym) - vx(ix ,iym))**2  &
              + (vy(ixm,iym) - vy(ix ,iym))**2 + (vz(ixm,iym) - vz(ix ,iym))**2) 
        a2 = ((vx(ix ,iym) - vx(ix ,iy ))**2  &
              + (vy(ix ,iym) - vy(ix ,iy ))**2 + (vz(ix ,iym) - vz(ix ,iy ))**2)
        a3 = ((vx(ix ,iy ) - vx(ixm,iy ))**2  &
              + (vy(ix ,iy ) - vy(ixm,iy ))**2 + (vz(ix ,iy ) - vz(ixm,iy ))**2) 
        a4 = ((vx(ixm,iy ) - vx(ixm,iym))**2  &
              + (vy(ixm,iy ) - vy(ixm,iym))**2 + (vz(ixm,iy ) - vz(ixm,iym))**2)

        p_visc(ix,iy) = MAX(p_visc(ix,iy), - alpha1(ix,iy)*SQRT(a1)) 
        p_visc(ix,iy) = MAX(p_visc(ix,iy), - alpha2(ix,iy)*SQRT(a2)) 

        visc_heat(ix,iy) = &
            - 0.5_num * dyb(iy) * alpha1(ix ,iy ) * a1 &
            - 0.5_num * dxb(ix) * alpha2(ix ,iy ) * a2 &
            - 0.5_num * dyb(iy) * alpha1(ix ,iyp) * a3 &
            - 0.5_num * dxb(ix) * alpha2(ixm,iy ) * a4

        visc_heat(ix,iy) = visc_heat(ix,iy) / cv(ix,iy)
      END DO
    END DO

    fx_visc = 0.0_num
    fy_visc = 0.0_num
    fz_visc = 0.0_num
    DO iy = 0, ny
      iym = iy - 1
      iyp = iy + 1
      DO ix = 0, nx
        ixm = ix - 1
        ixp = ix + 1

        a1 = alpha1(ix ,iyp) * dyc(iy)
        a2 = alpha1(ixp,iyp) * dyc(iy)
        a3 = alpha2(ix ,iy ) * dxc(ix)
        a4 = alpha2(ix ,iyp) * dxc(ix)

        fx_visc(ix,iy) = (a1 * (vx(ix,iy) - vx(ixm,iy )) &
                        + a2 * (vx(ix,iy) - vx(ixp,iy )) &
                        + a3 * (vx(ix,iy) - vx(ix ,iym)) &
                        + a4 * (vx(ix,iy) - vx(ix ,iyp)) ) / cv_v(ix,iy)

        fy_visc(ix,iy) = (a1 * (vy(ix,iy) - vy(ixm,iy )) &
                        + a2 * (vy(ix,iy) - vy(ixp,iy )) &
                        + a3 * (vy(ix,iy) - vy(ix ,iym)) &
                        + a4 * (vy(ix,iy) - vy(ix ,iyp)) ) / cv_v(ix,iy)

        fz_visc(ix,iy) = (a1 * (vz(ix,iy) - vz(ixm,iy )) &
                        + a2 * (vz(ix,iy) - vz(ixp,iy )) &
                        + a3 * (vz(ix,iy) - vz(ix ,iym)) &
                        + a4 * (vz(ix,iy) - vz(ix ,iyp)) ) / cv_v(ix,iy)

      END DO
    END DO

    DEALLOCATE(cs, cs_v)

  CONTAINS

    DOUBLE PRECISION FUNCTION edge_viscosity()

      ! Actually returns q_k_bar = q_kur*(1-psi) / abs(dv)
      ! Other symbols follow notation in Caramana

      REAL(num) :: dvx, dvy, dvz, dv, dv2
#ifdef SHOCKLIMITER      
      REAL(num) :: dvxm, dvxp, dvym, dvyp, dvzm, dvzp
      REAL(num) :: rl, rr
#endif
      REAL(num) :: psi, rho_edge, cs_edge, q_k_bar

#ifdef SHOCKEXPANSION
      ! Allow shock viscoity on expanding edge
      dvdots = -ABS(dvdots) 
#else
      ! Turn off shock viscosity if cell edge expanding
      dvdots = MIN(0.0_num, dvdots)
#endif

      rho_edge = 2.0_num * rho_v(i1,j1) * rho_v(i2,j2) &
          / (rho_v(i1,j1) + rho_v(i2,j2))
      cs_edge = MIN(cs_v(i1,j1), cs_v(i2,j2))

      dvx = vx(i1,j1) - vx(i2,j2)
      dvy = vy(i1,j1) - vy(i2,j2)
      dvz = vz(i1,j1) - vz(i2,j2)
      dv2 = dvx**2 + dvy**2 + dvz**2
      dv = SQRT(dv2)
      psi = 0.0_num
      IF (dv * dt / dx < 1.e-14_num) THEN
        dvdots = 0.0_num
      ELSE
        dvdots = dvdots / dv
      END IF

#ifdef SHOCKLIMITER
      dvxm = vx(i0,j0) - vx(i1,j1)
      dvxp = vx(i2,j2) - vx(i3,j3)
      dvym = vy(i0,j0) - vy(i1,j1)
      dvyp = vy(i2,j2) - vy(i3,j3)
      dvzm = vz(i0,j0) - vz(i1,j1)
      dvzp = vz(i2,j2) - vz(i3,j3)
      IF (dv * dt / dx < 1.e-14_num) THEN
        rl = 1.0_num
        rr = 1.0_num
      ELSE
        rl = (dvxp * dvx + dvyp * dvy + dvzp * dvz) * dx / (dxp * dv2)
        rr = (dvxm * dvx + dvym * dvy + dvzm * dvz) * dx / (dxm * dv2)
      END IF
      psi = MIN(0.5_num * (rr + rl), 2.0_num * rl, 2.0_num * rr, 1.0_num)
      psi = MAX(0.0_num, psi)
#endif

      ! Find q_kur / abs(dv)
      q_k_bar = rho_edge &
          * (visc2_norm * dv + SQRT(visc2_norm**2 * dv2 + (visc1 * cs_edge)**2))

      edge_viscosity = q_k_bar * (1.0_num - psi) * dvdots

    END FUNCTION edge_viscosity

  END SUBROUTINE edge_shock_viscosity



  SUBROUTINE edge_shock_heating

    REAL(num) :: a1, a2, a3, a4

    visc_heat = 0.0_num

    DO iy = 0, ny + 1 
      iym = iy - 1
      iyp = iy + 1
      DO ix = 0, nx + 1 
        ixm = ix - 1
        ixp = ix + 1

        a1 =  (vx(ixm,iym) - vx(ix ,iym))*(vx1(ixm,iym) - vx1(ix ,iym)) &
            + (vy(ixm,iym) - vy(ix ,iym))*(vy1(ixm,iym) - vy1(ix ,iym)) &
            + (vz(ixm,iym) - vz(ix ,iym))*(vz1(ixm,iym) - vz1(ix ,iym)) 
        a2 =  (vx(ix ,iym) - vx(ix ,iy ))*(vx1(ix ,iym) - vx1(ix ,iy )) &
            + (vy(ix ,iym) - vy(ix ,iy ))*(vy1(ix ,iym) - vy1(ix ,iy )) &
            + (vz(ix ,iym) - vz(ix ,iy ))*(vz1(ix ,iym) - vz1(ix ,iy ))
        a3 =  (vx(ix ,iy ) - vx(ixm,iy ))*(vx1(ix ,iy ) - vx1(ixm,iy )) &
            + (vy(ix ,iy ) - vy(ixm,iy ))*(vy1(ix ,iy ) - vy1(ixm,iy )) &
            + (vz(ix ,iy ) - vz(ixm,iy ))*(vz1(ix ,iy ) - vz1(ixm,iy ))
        a4 =  (vx(ixm,iy ) - vx(ixm,iym))*(vx1(ixm,iy ) - vx1(ixm,iym)) &
            + (vy(ixm,iy ) - vy(ixm,iym))*(vy1(ixm,iy ) - vy1(ixm,iym)) &
            + (vz(ixm,iy ) - vz(ixm,iym))*(vz1(ixm,iy ) - vz1(ixm,iym))

        visc_heat(ix,iy) = &
            - 0.5_num * dyb(iy) * alpha1(ix,iy) * a1 &
            - 0.5_num * dxb(ix) * alpha2(ix,iy) * a2 &
            - 0.5_num * dyb(iy) * alpha1(ix,iyp) * a3 &
            - 0.5_num * dxb(ix) * alpha2(ixm,iy) * a4

        visc_heat(ix,iy) = visc_heat(ix,iy) / cv(ix,iy)
      END DO
    END DO

    visc_heat = MAX(visc_heat, 0.0_num)

  END SUBROUTINE edge_shock_heating



  SUBROUTINE b_field_and_cv1_update

    REAL(num) :: vxb, vxbm, vyb, vybm, dvxdx, dvydy, dv
#ifdef CAUCHY
    REAL(num) :: dvxdy, dvydx, dvzdx, vzb, vzbm, dvzdy
#endif

    DO iy = -1, ny + 2
      iym = iy - 1
      DO ix = -1, nx + 2
        ixm = ix - 1

        ! vx at Bx(i,j)
        vxb  = (vx(ix ,iy ) + vx(ix ,iym)) * 0.5_num
        ! vx at Bx(i-1,j)
        vxbm = (vx(ixm,iy ) + vx(ixm,iym)) * 0.5_num
        ! vy at By(i,j)
        vyb  = (vy(ix ,iy ) + vy(ixm,iy )) * 0.5_num
        ! vy at By(i,j-1)
        vybm = (vy(ix ,iym) + vy(ixm,iym)) * 0.5_num

        dvxdx = (vxb - vxbm) / dxb(ix)
        dvydy = (vyb - vybm) / dyb(iy)

        dv = (dvxdx + dvydy) * dt2
        cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)
#ifdef CAUCHY
        ! vx at By(i,j)
        vxb  = (vx(ix ,iy ) + vx(ixm,iy )) * 0.5_num
        ! vx at By(i,j-1)
        vxbm = (vx(ix ,iym) + vx(ixm,iym)) * 0.5_num
        ! vy at Bx(i,j)
        vyb  = (vy(ix ,iy ) + vy(ix ,iym)) * 0.5_num
        ! vy at Bx(i-1,j)
        vybm = (vy(ixm,iy ) + vy(ixm,iym)) * 0.5_num

        dvxdy = (vxb - vxbm) / dyb(iy)
        dvydx = (vyb - vybm) / dxb(ix)

        ! vz at Bx(i,j)
        vzb  = (vz(ix ,iy ) + vz(ix ,iym)) * 0.5_num
        ! vz at Bx(i-1,j)
        vzbm = (vz(ixm,iy ) + vz(ixm,iym)) * 0.5_num
        ! vz at By(i,j)
        dvzdx = (vzb - vzbm) / dxb(ix)

        vzb  = (vz(ix ,iy ) + vz(ixm,iy )) * 0.5_num
        ! vz at By(i,j-1)
        vzbm = (vz(ix ,iym) + vz(ixm,iym)) * 0.5_num
        dvzdy = (vzb - vzbm) / dyb(iy)

        w3 =  bx1(ix,iy) * dvxdx + by1(ix,iy) * dvxdy
        w4 =  bx1(ix,iy) * dvydx + by1(ix,iy) * dvydy
        w5 =  bx1(ix,iy) * dvzdx + by1(ix,iy) * dvzdy

        bx1(ix,iy) = (bx1(ix,iy) + w3 * dt2) / (1.0_num + dv)
        by1(ix,iy) = (by1(ix,iy) + w4 * dt2) / (1.0_num + dv)
        bz1(ix,iy) = (bz1(ix,iy) + w5 * dt2) / (1.0_num + dv)
#endif
      END DO
    END DO

#ifndef CAUCHY
    vx1(:,:) = vx(:,:)
    vy1(:,:) = vy(:,:)
    vz1(:,:) = vz(:,:)
    
    predictor_step = .TRUE.
    dt = 0.5_num * dt
    CALL eulerian_remap(step)
    dt = 2.0_num * dt
    predictor_step = .FALSE. 
#endif

  END SUBROUTINE b_field_and_cv1_update



  !****************************************************************************
  ! Sets CFL limited step
  !****************************************************************************

  SUBROUTINE set_dt

    ! Assumes all variables are defined at the same point. Be careful with
    ! setting 'dt_multiplier' if you expect massive changes across cells.

    REAL(num) :: cs2, c_visc2, rho0, length
    REAL(num) :: dxlocal, dt_local, dtr, dtr_local, dth, dth_local
    REAL(num) :: dt1, dt2, dt3, dt4, ss_reduct_fac
    REAL(num) :: dt_locals(3), dt_min(3)
    REAL(num) :: dt0, time_dump, time_rem
    REAL(num) :: dt_fudge = 1e-4_num
    CHARACTER(LEN=1) :: dt_reason
    LOGICAL :: is_restart = .FALSE.
    LOGICAL, SAVE :: first = .TRUE.
    ! Maximum number of super-steps
    INTEGER, PARAMETER :: ss_limit = 60

    IF (first) THEN
      first = .FALSE.
      IF (restart) THEN
        dt = dt_from_restart
        RETURN
      END IF
    END IF

    dt_local = largest_number
    dtr_local = largest_number
    dth_local = largest_number

    gamma_boris(:,:) = 1.0_num

    DO iy = 0, ny
      iym = iy - 1
      DO ix = 0, nx
        ixm = ix - 1

        ! Fix dt for Lagrangian step
        w1 = bx1(ix,iy)**2 + by1(ix,iy)**2 + bz1(ix,iy)**2
        ! Sound speed squared
        rho0 = MAX(rho(ix,iy), none_zero)
        w2 = w1 / rho0
        IF (boris .AND. (w2 .GE. va_max2)) THEN
          gamma_boris(ix,iy) = 1.0_num / (1.0_num + w2 / va_max2)
        END IF
        cs2 =  (gamma * pressure(ix,iy) +  gamma_boris(ix,iy) * w1) / rho0

        !effective speed from viscous pressure
        c_visc2 = p_visc(ix,iy) / rho0

        ! length based on simple DYNA2D estimates
        length = dxb(ix) * dyb(iy) / SQRT(dxb(ix)**2 + dyb(iy)**2)

        ! Find ideal MHD CFL limit for Lagrangian step
        dt1 = length / (SQRT(c_visc2) + SQRT(cs2 + c_visc2))
        dt_local = MIN(dt_local, dt1)

        ! Check no node moves more than one cell to allow remap
        dt2 = MIN(dxb(ix) / MAX(ABS(vx(ix,iy)),none_zero), &
              dyb(iy) / MAX(ABS(vy(ix,iy)),none_zero))
        dt_local = MIN(dt_local, dt2)

        ! Note resistive limits assumes uniform resistivity hence cautious
        ! factor 0.2
        dxlocal = 1.0_num / (1.0_num / dxb(ix)**2 + 1.0_num / dyb(iy)**2)

        IF (cowling_resistivity) THEN
          dt3 = 0.2_num * dxlocal &
              / MAX(MAX(eta(ix,iy), eta_perp(ix,iy)), none_zero)
        ELSE
          dt3 = 0.2_num * dxlocal / MAX(eta(ix,iy), none_zero)
        END IF

        ! Adjust to accomodate resistive effects
        dtr_local = MIN(dtr_local, dt3)

        ! Hall MHD CFL limit
        IF (hall_mhd) THEN
          dt4 = 0.75_num * rho(ix,iy) * MIN(dxb(ix), dyb(iy))**2 &
            / MAX(lambda_i(ix,iy) * SQRT(w1), none_zero)
          dth_local = MIN(dth_local, dt4)
        END IF
        
      END DO
    END DO

    dt_locals(1) = dt_local
    dt_locals(2) = dtr_local
    dt_locals(3) = dth_local

    CALL MPI_ALLREDUCE(dt_locals, dt_min, 3, mpireal, MPI_MIN, comm, errcode)

    dt  = dt_multiplier * dt_min(1)
    dtr = dt_multiplier * dt_min(2)
    dth = dt_multiplier * dt_min(3)


    time = time + dt

  END SUBROUTINE set_dt


END MODULE lagran
