!     Steady State function that builds pl curve
      subroutine steady_state_pl(Nt, times, x, rate, pl)
      Implicit None
      Integer*4 Nt,i
      Integer*4, parameter :: Nt_guess=1000, Nr_guess=30
      Real*8 rate(Nr_guess),x(3),dx(3),t_old,t_new
      Real*8 times(Nt_guess)
      Real*8, intent(out) :: pl(Nt_guess)

      dx = 0d0
      t_old = 0
      do i = 1, Nt ! Loop over experimental time steps
          t_new = times(i)
          dx = 0d0
          Call Kinetic(x,rate,dx) ! Get time derivatives from kinetic equations
          x = x + dx * (t_new - t_old) ! Update the concentrations of each species using Euler Integration
          pl(i) = x(1)
          t_old = t_new
      end do

      end subroutine steady_state_pl


!     Steady State function that builds rolloff curve
      subroutine steady_state_el(NJ, Jvals, x_i, rate, eqe, jout, fail)
      Implicit None
      Integer*4 i,J,NJ
      Integer*4, parameter :: NJ_guess=1000, Nr_guess=30
      Real*8 rate(Nr_guess),flux(3),x(3),dx(3),t,dt
      Real*8 jvals(NJ_guess),x_i(3),dx_old(3),x_old(3)
      Real*8, intent(out) :: eqe(NJ_guess)
      Real*8, intent(out) :: Jout(NJ_guess)
      Integer*4, intent(out) :: fail

      x = x_i
      flux = 0d0
      fail = 0
      do J = 1, NJ ! Loop over experimental current densities

          t=0d0; dt=1d-11 ! 10 ps timestep
          Flux(2) = Jvals(J) ! apply flux J to charges

          ! initialize variables
          x=0d0
          x=x_i
          dx=0d0
          dx=flux     ! apply flux to conc derivatives
          dx_old=1d3  ! set to arbitrary number
          x_old=1d3   ! set to arbitrary number

          ! Evolve system until steady state is reached
          i = 1
          do while (abs(dx_old(1)-dx(1))>1d-100)
              t = t + dt
              dx_old = dx
              dx=0d0
              dx=flux
              Call Kinetic(x,rate,dx) ! Get time derivatives from kinetic equations
              x = x + dx * dt ! Update the concentrations of each species using Euler Integration
              if (t>1d0) then ! Exit if steady state is not reached in 1 sec
                  exit
              end if
              i = i + 1
          end do
          if (t>1d0) then
              fail = 1
              print*,"t > 1 sec, steady state not found"
              exit
          end if
!          print*,"steady state reached at",i,"at t=",t,"iter=",J 

          Jout(J) = flux(2)
          eqe(J) = x(1) / (flux(2)*dt*1d3) ! Calculate EQE, converting back to units s
          flush(6)

      end do
      end subroutine steady_state_el


!     Kinetic equations function that calculates the time derivatives for EL
!     TODO: comment out equations not used in model, update rate indices accordingly
!     Note: DO NOT change the order of the equations
      Subroutine Kinetic(x, rate, dx)
      Implicit None
      Integer*4, parameter :: Nr_guess=30
      Real*8 rate(Nr_guess),x(3)
      Real*8, intent(out) :: dx(3)
      
!     dx => 1=Exciton; 2=Charge; 3=Dark exciton; 4=Defect
      
!     a K_abs: Exciton Absorption:   -> ex
      dx(1)=dx(1)+rate(1)*1d0

!     b K_rad: Exciton Radiative Decay (emission): ex ->
      dx(1)=dx(1)-rate(2)*x(1)

!     c K_ET: Exciton Nonradiative Decay: ex ->
      dx(1)=dx(1)-rate(3)*x(1)

!     d K_DT: Dark Exciton Nonradiative Decay: dx ->
      dx(3)=dx(3)-rate(4)*x(3)

!     e K_CT: Charge Leakage: el + ho -> 
      dx(2)=dx(2)-2d0*rate(5)*x(2)*x(2)

!     f K_BD: Exciton conversion, bright to dark: ex -> dx
      dx(1)=dx(1)-rate(6)*x(1)
      dx(3)=dx(3)+rate(6)*x(1)

!     g K_DB: Exciton conversion, dark to bright: dx -> ex
      dx(1)=dx(1)+rate(7)*x(3)
      dx(3)=dx(3)-rate(7)*x(3)

!     h K_sep_bright: Charge Separation: ex -> el + ho
      dx(1)=dx(1)-rate(8)*x(1)
      dx(2)=dx(2)+2d0*rate(8)*x(1)

!     i K_sep_dark: Charge Separation: dx -> el + ho
      dx(3)=dx(3)-rate(9)*x(3)
      dx(2)=dx(2)+2d0*rate(9)*x(3)

!     j K_rec_bright: Charge Recombination: el + ho -> ex
      dx(2)=dx(2)-2d0*rate(10)*x(2)*x(2)
      dx(1)=dx(1)+rate(10)*x(2)*x(2)
      
!     k K_rec_dark: Charge Recombination: el + ho -> dx
      dx(2)=dx(2)-2d0*rate(11)*x(2)*x(2)
      dx(3)=dx(3)+rate(11)*x(2)*x(2)

!     l K_TTA: Triplet-Triplet Annihilation
!     ex + ex -> ex
!     ex + dx -> ex
!     ex + dx -> dx
!     dx + dx -> dx
      dx(1)=dx(1)-rate(12)*x(1)*x(1)
      dx(3)=dx(3)-rate(12)*x(1)*x(3)
      dx(1)=dx(1)-rate(12)*x(1)*x(3)
      dx(3)=dx(3)-rate(12)*x(3)*x(3)
      
!     m K_TCA: Triplet-Charge Annihilation
!     ex + ch -> ch
!     dx + ch -> ch
      dx(1)=dx(1)-rate(13)*x(1)*x(2)
      dx(3)=dx(3)-rate(13)*x(3)*x(2)

      End Subroutine Kinetic
