! Fortran implementation of Lorenz system with Euler integration
! Created on Tue Feb 6 15:05:51 2020
! @author: Tammas Loughran
! @email: t.loughran@lmu.de
program lorenz
    use ogpf
    implicit none
    ! Parameters
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: N = 10000
    double precision, parameter :: H = 0.01_dp
    double precision, parameter :: sigma = 10, rho = 28, beta = 8.0_dp/3.0_dp
    integer :: t
    ! Local variables
    double precision, dimension(N) :: xstate, ystate, zstate
    double precision :: dxdt, dydt, dzdt
    type(gpf) :: gp

    ! Initialise the state arrays with 1
    xstate(1) = 1
    ystate(1) = 1
    zstate(1) = 1

    ! Loop over all time steps
    step: do t = 2, N
        ! Calculate the tendencies
        dxdt = xtend(xstate(t-1), ystate(t-1))
        dydt = ytend(xstate(t-1), ystate(t-1), zstate(t-1))
        dzdt = ztend(xstate(t-1), ystate(t-1), zstate(t-1))
        ! Integrate forward in time
        xstate(t) = euler(xstate(t-1), dxdt)
        ystate(t) = euler(ystate(t-1), dydt)
        zstate(t) = euler(zstate(t-1), dzdt)
    end do step

    ! Write state variables to file
    open(unit=1, file='data.dat')
    output: do t = 1, N
        write(1,*)  xstate(t), ystate(t), zstate(t)
    end do output
    close(unit=1)

    ! Plot the results
    call gp%lplot(xstate, ystate, zstate)

    contains
        double precision function xtend(x, y)
            ! Calculate the tendency of the x variable
            double precision :: x, y
            xtend = sigma * (y - x)
        end function xtend

        double precision function ytend(x, y, z)
            ! Calculate the tendency of the y variable
            double precision :: x, y, z
            ytend = x * (rho - z) - y
        end function ytend

        double precision function ztend(x, y, z)
            ! Calculate the tendency of the z variable
            double precision :: x, y, z
            ztend = x * y - beta * z
        end function ztend

        double precision function euler(yn, dyndt)
            ! Integrate the state forward in time
            double precision :: yn, dyndt
            euler = yn + dyndt * H
        end function euler
end program

