!Module for the two-fluid version of Lare2D
!Experimental version. Not for public use.
!Written by Ben Snow (Exeter University) - 2025


MODULE two_fluid

  USE shared_data
  !USE boundary
  USE xremap_neutral
  USE yremap_neutral
  USE zremap_neutral

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: setup_two_fluid, lagrangian_neutral

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE lagrangian_neutral
        !Neutral lagrangian step

        ALLOCATE(alpha1(0:nx+1,0:ny+2))
        ALLOCATE(alpha2(-1:nx+1,0:ny+1))
        ALLOCATE(visc_heat(0:nx+1,0:ny+1))
        ALLOCATE(pressure_n(-1:nx+2,-1:ny+2))
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
        
        
        !Get the pressure from the energy
        DO iy = -1, ny + 2
          DO ix = -1, nx + 2
            pressure_n(ix,iy) = (gamma - 1.0_num) * rho_n(ix,iy) &
                * (energy_n(ix,iy))
          END DO
        END DO    
        
        !get the cell volume?
        DO iy = -1, ny + 1
          iyp = iy + 1
          DO ix = -1, nx + 1
            ixp = ix + 1
            rho_v(ix,iy) = rho_n(ix,iy) * cv(ix,iy) + rho_n(ixp,iy) * cv(ixp,iy) &
                +   rho_n(ix,iyp) * cv(ix,iyp) + rho_n(ixp,iyp) * cv(ixp,iyp)
            cv_v(ix,iy) = cv(ix,iy) + cv(ixp,iy) + cv(ix,iyp) + cv(ixp,iyp)
            rho_v(ix,iy) = rho_v(ix,iy) / cv_v(ix,iy)

            cv_v(ix,iy) = 0.25_num * cv_v(ix,iy)
          END DO
        END DO
        
        !Viscous heating in the neutral fluid
        CALL edge_shock_viscosity
        
        !!!!!!!!!!!!!!!!!!!!!!!!
        ! MISSING DT CALCULATION HERE
        !!!!!!!!!!!!!!!!!!!!!!!!
        
        !
        energy_n(:,:) = MAX(energy_n(:,:), 0.0_num)        
        
        CALL predictor_corrector_step

        DEALLOCATE(alpha1, alpha2)
        DEALLOCATE(visc_heat, pressure_n, rho_v, cv_v, flux_x, flux_y, flux_z)
        DEALLOCATE(qx, qy, qz)
        DEALLOCATE(fx_visc, fy_visc, fz_visc)

        CALL energy_bcs
        CALL density_bcs
        CALL velocity_bcs
    
    END SUBROUTINE lagrangian_neutral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE eulerian_remap_neutral(i)

    INTEGER, INTENT(IN) :: i
    INTEGER :: case_test

    IF (rke) delta_ke = 0.0_num
    xpass = 1
    ypass = 1

    case_test = MODULO(i, 6)

    ! Strang ordering
    SELECT CASE(case_test)
    CASE (0)
      CALL remap_x
      CALL remap_y
      CALL remap_z
    CASE (1)
      CALL remap_y
      CALL remap_z
      CALL remap_x
    CASE (2)
      CALL remap_z
      CALL remap_x
      CALL remap_y
    CASE (3)
      CALL remap_x
      CALL remap_z
      CALL remap_y
    CASE (4)
      CALL remap_z
      CALL remap_y
      CALL remap_x
    CASE (5)
      CALL remap_y
      CALL remap_x
      CALL remap_z
    END SELECT

  END SUBROUTINE eulerian_remap_neutral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    SUBROUTINE edge_shock_viscosity

        REAL(num) :: dvdots, dx, dxm, dxp
        REAL(num) :: rmin
        REAL(num) :: a1, a2, a3, a4
        REAL(num), DIMENSION(:,:), ALLOCATABLE :: cs, cs_v
        INTEGER :: i0, i1, i2, i3, j0, j1, j2, j3
        LOGICAL, SAVE :: first_call = .TRUE.

        ALLOCATE(cs(-1:nx+2,-1:ny+2), cs_v(-1:nx+1,-1:ny+1))

        IF (first_call) THEN
          first_call = .FALSE.
          visc2_norm = 0.25_num * (gamma + 1.0_num) * visc2
        END IF

        p_visc = 0.0_num
        visc_heat = 0.0_num

        DO iy = -1, ny + 2
          DO ix = -1, nx + 2
            rmin = MAX(rho_n(ix,iy), none_zero)
            cs(ix,iy) = SQRT((gamma * pressure_n(ix,iy)) / rmin)
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
          REAL(num) :: psi, rho_edge, cs_edge, q_k_bar

          ! Turn off shock viscosity if cell edge expanding
          dvdots = MIN(0.0_num, dvdots)

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

          ! Find q_kur / abs(dv)
          q_k_bar = rho_edge &
              * (visc2_norm * dv + SQRT(visc2_norm**2 * dv2 + (visc1 * cs_edge)**2))

          edge_viscosity = q_k_bar * (1.0_num - psi) * dvdots

        END FUNCTION edge_viscosity

    END SUBROUTINE edge_shock_viscosity
    
  !****************************************************************************
  ! The main predictor / corrector step which advances the momentum equation
  !****************************************************************************

  SUBROUTINE predictor_corrector_step

    REAL(num) :: pp, ppx, ppy, ppxy
    REAL(num) :: e1
    REAL(num) :: vxb, vxbm, vyb, vybm
    REAL(num) :: bxv, byv, bzv, jx, jy, jz
    REAL(num) :: dv
    REAL(num) :: fx, fy, fz

    CALL cv1_update

    DO iy = 0, ny + 1
      DO ix = 0, nx + 1
        dv = cv1(ix,iy) / cv(ix,iy) - 1.0_num
        ! Predictor energy
        e1 = energy_n(ix,iy) - pressure_n(ix,iy) * dv / rho_n(ix,iy)
        e1 = e1 + visc_heat(ix,iy) * dt2 / rho_n(ix,iy)

        ! Now define the predictor step pressures
        pressure_n(ix,iy) = e1* (gamma - 1.0_num) * rho_n(ix,iy) * cv(ix,iy) / cv1(ix,iy)
      END DO
    END DO

    DO iy = 0, ny
      iyp = iy + 1
      iym = iy - 1
      DO ix = 0, nx
        ixp = ix + 1
        ixm = ix - 1

        pp    = pressure_n(ix ,iy )
        ppx   = pressure_n(ixp,iy )
        ppy   = pressure_n(ix ,iyp)
        ppxy  = pressure_n(ixp,iyp)

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

        fy = fy - rho_v(ix,iy) * grav(iy)

        ! Find half step velocity needed for remap
        vn_x1(ix,iy) = vn_x(ix,iy) + dt2 * (fx_visc(ix,iy) + fx) / rho_v(ix,iy)
        vn_y1(ix,iy) = vn_y(ix,iy) + dt2 * (fy_visc(ix,iy) + fy) / rho_v(ix,iy)
        vn_z1(ix,iy) = vn_z(ix,iy) + dt2 * (fz_visc(ix,iy) + fz) / rho_v(ix,iy)
      END DO
    END DO

    
    CALL remap_v_bcs

    CALL edge_shock_heating

    DO iy = 0, ny
      DO ix = 0, nx
        ! Velocity at the end of the Lagrangian step
        vn_x(ix,iy) = 2.0_num * vn_x1(ix,iy) - vn_x(ix,iy) 
        vn_y(ix,iy) = 2.0_num * vn_y1(ix,iy) - vn_y(ix,iy) 
        vn_z(ix,iy) = 2.0_num * vn_z1(ix,iy) - vn_z(ix,iy) 
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
        vxb  = (vn_x1(ix ,iy ) + vn_x1(ix ,iym)) * 0.5_num
        ! vx1 at Bx(i-1,j)
        vxbm = (vn_x1(ixm,iy ) + vn_x1(ixm,iym)) * 0.5_num
        ! vy1 at By(i,j)
        vyb  = (vn_y1(ix ,iy ) + vn_y1(ixm,iy )) * 0.5_num
        ! vy1 at By(i,j-1)
        vybm = (vn_y1(ix ,iym) + vn_y1(ixm,iym)) * 0.5_num

        dv = ((vxb - vxbm) / dxb(ix) + (vyb - vybm) / dyb(iy)) * dt

        cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)

        ! Energy at end of Lagrangian step
        energy_n(ix,iy) = energy_n(ix,iy) &
            + (dt * visc_heat(ix,iy) - dv * pressure_n(ix,iy)) &
            / rho_n(ix,iy)

        visc_dep(ix,iy) = visc_dep(ix,iy) + dt * visc_heat(ix,iy)

        rho_n(ix,iy) = rho_n(ix,iy) / (1.0_num + dv)

        total_visc_heating = total_visc_heating &
            + dt * visc_heat(ix,iy) * cv(ix,iy)

      END DO
    END DO

  END SUBROUTINE predictor_corrector_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  SUBROUTINE cv1_update

    REAL(num) :: vxb, vxbm, vyb, vybm, dvxdx, dvydy, dv

    DO iy = -1, ny + 2
      iym = iy - 1
      DO ix = -1, nx + 2
        ixm = ix - 1

        ! vx at Bx(i,j)
        vxb  = (vn_x(ix ,iy ) + vn_x(ix ,iym)) * 0.5_num
        ! vx at Bx(i-1,j)
        vxbm = (vx(ixm,iy ) + vn_x(ixm,iym)) * 0.5_num
        ! vy at By(i,j)
        vyb  = (vn_y(ix ,iy ) + vn_y(ixm,iy )) * 0.5_num
        ! vy at By(i,j-1)
        vybm = (vn_y(ix ,iym) + vn_y(ixm,iym)) * 0.5_num

        dvxdx = (vxb - vxbm) / dxb(ix)
        dvydy = (vyb - vybm) / dyb(iy)

        dv = (dvxdx + dvydy) * dt2
        cv1(ix,iy) = cv(ix,iy) * (1.0_num + dv)

      END DO
    END DO

    vn_x1(:,:) = vn_x(:,:)
    vn_y1(:,:) = vn_y(:,:)
    vn_z1(:,:) = vn_z(:,:)
    
    predictor_step = .TRUE.
    dt = 0.5_num * dt
    CALL eulerian_remap_neutral(step)
    dt = 2.0_num * dt
    predictor_step = .FALSE. 


  END SUBROUTINE cv1_update
    
END MODULE two_fluid
