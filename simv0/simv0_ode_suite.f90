!=============================================================================
!  simv0_ode_suite.f90  --  Unified ODE Integration Suite
!
!  Provides four integration methods behind a single abstract interface so
!  that any physics RHS can be passed in without code duplication.
!
!  Methods implemented:
!    euler_step   – Forward Euler          O(h^1)  1 RHS evaluation
!    heun_step    – Improved Euler / Heun  O(h^2)  2 RHS evaluations
!    rk4_step     – Classical RK4          O(h^4)  4 RHS evaluations
!    rk45_step    – Dormand-Prince RK45    O(h^5)  7 RHS evaluations
!                   Returns 5th-order solution + WRMS error norm for
!                   adaptive step-size control (step accepted if err_norm ≤ 1)
!
!  Abstract interface for the right-hand side function:
!
!    subroutine rhs(t, y, n, dydt)
!        real(8), intent(in)  :: t        !  current time (s)
!        integer, intent(in)  :: n        !  number of state variables
!        real(8), intent(in)  :: y(n)     !  state vector at t
!        real(8), intent(out) :: dydt(n)  !  dy/dt = f(t,y)
!    end subroutine
!
!  Dormand-Prince (DOPRI5) Butcher tableau (non-FSAL variant for clarity):
!    c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1]
!    5th-order solution uses stages k1..k6 (b weights = a71..a76, a72=0)
!    Error  = y5 - y4  = h * (e1*k1 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7)
!    where k7 is evaluated at (t+h, y5) – the FSAL starting stage.
!
!  Compile (after simv0_config):  gfortran -O2 -c simv0_ode_suite.f90
!=============================================================================
module simv0_ode_suite
    use simv0_config
    implicit none

    !--------------------------------------------------------------------------
    !  Abstract interface for the RHS subroutine dy/dt = f(t,y)
    !--------------------------------------------------------------------------
    abstract interface
        subroutine rhs_interface(t, y, n, dydt)
            real(8), intent(in)  :: t
            integer, intent(in)  :: n
            real(8), intent(in)  :: y(n)
            real(8), intent(out) :: dydt(n)
        end subroutine rhs_interface
    end interface

    !--------------------------------------------------------------------------
    !  Dormand-Prince RK45 Butcher tableau  –  internal lower-triangular coefficients
    !--------------------------------------------------------------------------
    !  Stage 2
    real(8), parameter :: dp_a21 = 1.0d0/5.0d0
    !  Stage 3
    real(8), parameter :: dp_a31 = 3.0d0/40.0d0
    real(8), parameter :: dp_a32 = 9.0d0/40.0d0
    !  Stage 4
    real(8), parameter :: dp_a41 =  44.0d0/45.0d0
    real(8), parameter :: dp_a42 = -56.0d0/15.0d0
    real(8), parameter :: dp_a43 =  32.0d0/9.0d0
    !  Stage 5
    real(8), parameter :: dp_a51 =  19372.0d0/6561.0d0
    real(8), parameter :: dp_a52 = -25360.0d0/2187.0d0
    real(8), parameter :: dp_a53 =  64448.0d0/6561.0d0
    real(8), parameter :: dp_a54 =   -212.0d0/729.0d0
    !  Stage 6
    real(8), parameter :: dp_a61 =   9017.0d0/3168.0d0
    real(8), parameter :: dp_a62 =   -355.0d0/33.0d0
    real(8), parameter :: dp_a63 =  46732.0d0/5247.0d0
    real(8), parameter :: dp_a64 =     49.0d0/176.0d0
    real(8), parameter :: dp_a65 =  -5103.0d0/18656.0d0
    !  5th-order output weights  b = [a71..a76]  (a72 = 0 by construction)
    real(8), parameter :: dp_b1  =    35.0d0/384.0d0
    real(8), parameter :: dp_b3  =   500.0d0/1113.0d0
    real(8), parameter :: dp_b4  =   125.0d0/192.0d0
    real(8), parameter :: dp_b5  = -2187.0d0/6784.0d0
    real(8), parameter :: dp_b6  =    11.0d0/84.0d0
    !  Error coefficients  e_i = b5_i - b4_i  (5th-order minus embedded 4th-order)
    real(8), parameter :: dp_e1  =    71.0d0/57600.0d0
    real(8), parameter :: dp_e3  =   -71.0d0/16695.0d0
    real(8), parameter :: dp_e4  =    71.0d0/1920.0d0
    real(8), parameter :: dp_e5  = -17253.0d0/339200.0d0
    real(8), parameter :: dp_e6  =    22.0d0/525.0d0
    real(8), parameter :: dp_e7  =    -1.0d0/40.0d0

contains

    !--------------------------------------------------------------------------
    !  euler_step  –  Forward Euler  (O(h^1) global truncation error)
    !
    !    y_new = y + h * f(t, y)
    !--------------------------------------------------------------------------
    subroutine euler_step(rhs, t, y, n, h, y_new)
        procedure(rhs_interface) :: rhs
        real(8), intent(in)  :: t, h
        integer, intent(in)  :: n
        real(8), intent(in)  :: y(n)
        real(8), intent(out) :: y_new(n)

        real(8) :: dydt(n)

        call rhs(t, y, n, dydt)
        y_new = y + h * dydt

    end subroutine euler_step

    !--------------------------------------------------------------------------
    !  heun_step  –  Improved Euler / Heun  (O(h^2) global truncation error)
    !
    !    k1    = f(t,     y)
    !    k2    = f(t + h, y + h*k1)
    !    y_new = y + (h/2) * (k1 + k2)
    !--------------------------------------------------------------------------
    subroutine heun_step(rhs, t, y, n, h, y_new)
        procedure(rhs_interface) :: rhs
        real(8), intent(in)  :: t, h
        integer, intent(in)  :: n
        real(8), intent(in)  :: y(n)
        real(8), intent(out) :: y_new(n)

        real(8) :: k1(n), k2(n)

        call rhs(t,     y,          n, k1)
        call rhs(t + h, y + h * k1, n, k2)
        y_new = y + (h / 2.0d0) * (k1 + k2)

    end subroutine heun_step

    !--------------------------------------------------------------------------
    !  rk4_step  –  Classical 4th-order Runge-Kutta  (O(h^4) global error)
    !
    !    k1    = f(t,       y)
    !    k2    = f(t + h/2, y + (h/2)*k1)
    !    k3    = f(t + h/2, y + (h/2)*k2)
    !    k4    = f(t + h,   y +  h   *k3)
    !    y_new = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    !--------------------------------------------------------------------------
    subroutine rk4_step(rhs, t, y, n, h, y_new)
        procedure(rhs_interface) :: rhs
        real(8), intent(in)  :: t, h
        integer, intent(in)  :: n
        real(8), intent(in)  :: y(n)
        real(8), intent(out) :: y_new(n)

        real(8) :: k1(n), k2(n), k3(n), k4(n)

        call rhs(t,           y,                   n, k1)
        call rhs(t + h/2.0d0, y + (h/2.0d0) * k1,  n, k2)
        call rhs(t + h/2.0d0, y + (h/2.0d0) * k2,  n, k3)
        call rhs(t + h,       y + h          * k3,  n, k4)
        y_new = y + (h / 6.0d0) * (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)

    end subroutine rk4_step

    !--------------------------------------------------------------------------
    !  rk45_step  –  Dormand-Prince single step (no adaptation here)
    !
    !  Computes the 5th-order solution y_new and the WRMS error norm used for
    !  adaptive step-size control.  The adaptive loop lives in the driver so
    !  that the driver can write output at every accepted step.
    !
    !  Step is accepted when  err_norm ≤ 1.0.
    !  New step size hint:  h_new = h * safety * (1 / err_norm)^(1/5)
    !
    !  Inputs:  rhs, t, y, n, h
    !  Outputs: y_new (5th-order), err_norm (WRMS, dimensionless)
    !--------------------------------------------------------------------------
    subroutine rk45_step(rhs, t, y, n, h, y_new, err_norm)
        procedure(rhs_interface) :: rhs
        real(8), intent(in)  :: t, h
        integer, intent(in)  :: n
        real(8), intent(in)  :: y(n)
        real(8), intent(out) :: y_new(n)
        real(8), intent(out) :: err_norm

        real(8) :: k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n)
        real(8) :: err(n), scale(n)
        integer :: i

        !  7 stage evaluations  (Dormand-Prince DOPRI5)
        call rhs(t,                                                    &
                 y,                                                    n, k1)
        call rhs(t + h/5.0d0,                                         &
                 y + h*(dp_a21*k1),                                    n, k2)
        call rhs(t + 3.0d0*h/10.0d0,                                  &
                 y + h*(dp_a31*k1 + dp_a32*k2),                        n, k3)
        call rhs(t + 4.0d0*h/5.0d0,                                   &
                 y + h*(dp_a41*k1 + dp_a42*k2 + dp_a43*k3),            n, k4)
        call rhs(t + 8.0d0*h/9.0d0,                                   &
                 y + h*(dp_a51*k1 + dp_a52*k2                         &
                      + dp_a53*k3 + dp_a54*k4),                        n, k5)
        call rhs(t + h,                                                &
                 y + h*(dp_a61*k1 + dp_a62*k2 + dp_a63*k3             &
                      + dp_a64*k4 + dp_a65*k5),                        n, k6)

        !  5th-order solution  (b weights; b2 = 0)
        y_new = y + h*(dp_b1*k1 + dp_b3*k3 + dp_b4*k4 + dp_b5*k5 + dp_b6*k6)

        !  7th stage evaluated at the new point (FSAL: first-same-as-last)
        call rhs(t + h, y_new, n, k7)

        !  Error vector  e = h * sum(e_i * k_i)  (5th minus embedded 4th order)
        err = h * (dp_e1*k1 + dp_e3*k3 + dp_e4*k4 + dp_e5*k5 + dp_e6*k6 + dp_e7*k7)

        !  WRMS (Weighted Root Mean Square) error norm
        !    scale_i = atol + rtol * max(|y_i|, |y_new_i|)
        do i = 1, n
            scale(i) = atol + rtol * max(abs(y(i)), abs(y_new(i)))
        end do
        err_norm = sqrt(sum((err / scale)**2) / real(n, 8))

    end subroutine rk45_step

end module simv0_ode_suite
