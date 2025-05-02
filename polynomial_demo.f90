! Author: Vegard G. Jervell
! Demo program for the Polynomials module
! Build using
!   mkdir build && cd build
!   cmake .. 
!   make polynomial_demo
! Will place the demo program in the directory: build/thermopack/polynomial_demo
! Building with `make install` will place the demo program in the "installed" directory

program polynomial_demo
    use polynomials
    implicit none
    
    integer :: N, m, i
    integer, allocatable, dimension(:, :) :: partitions
    integer, allocatable, dimension(:) :: partition_sizes
    real :: Z, dx, v2, v3

    type(Polynomial) :: poly1, poly2
    type(PolyExp) :: pexp
    type(PolyFrac) :: pfrac
    real :: x

    print*, " ------------- PARTITION DEMO STARTS ------------- "
    N = 5
    m = -1

    call get_partitions(N, m, partitions, partition_sizes)
    
    print*, "Partitions of ", N
    do i = 1, size(partitions, 1)
        Z = partition_multiplicity(partitions(i, :), partition_sizes(i))
        print*, Z, " : ", partitions(i, :)
    enddo
    print*, " ------------- PARTITION DEMO ENDS ------------- "
    print*, ""
    print*, " ------------- POLYNOMIAL DEMO STARTS ------------- "

    ! Polynomial 1: 1 + 2*x + 3*x^2
    call poly1%init(k_min=0, k_max=2, k_step=1, coeff=[1.0, 2.0, 3.0])
    ! call poly1%init(k_min=0, k_max=2, k_step=1, coeff=[1.0, 0., 0.])
    print*, "Polynomial 1: "
    call poly1%print()
    print*, ""

    x = 1
    
    print*, "Evaluate polynomial 1 at x = ", x, " : ", poly1%eval(x)
    do N = 0, 4
        print*, "Derivative ", N, " : ", poly1%derivative(x, N)
    enddo
    print*, ""
    ! Polynomial 2: 1 * x^{-2} + 2 + 3 * x^2
    call poly2%init(k_min=-2, k_max=2, k_step=2, coeff=[-1.0, 2.0, -3.0])
    ! call poly2%init(k_min=0, k_max=2, k_step=1, coeff=[1.0, 0., 0.])
    print*, "Polynomial 2: "
    call poly2%print()
    print*, ""
    print*, "Evaluate polynomial 2 at x = ", x, " : ", poly2%eval(x)
    do N = 0, 4
        print*, "Derivative ", N, " : ", poly2%derivative(x, N)
    enddo
    print*, " ------------- POLYNOMIAL DEMO ENDS -------------"
    print*, ""
    print*, " ------------- POLYEXP DEMO STARTS -------------"

    ! PolyExp object : f(x) * exp(g(x)), where f and g are poly1 and poly2
    call pexp%init(pref=poly1, expo=poly2)
    print*, "PolyExp object : "
    call pexp%print()
    print*, ""
    print*, "Eval PolyExp at x = ", x, " : ", pexp%eval(x)
    dx = 1e-6
    do N = 0, 4
        if (N > 0) then
            v2 = pexp%derivative(x - dx, N - 1)
            v3 = pexp%derivative(x + dx, N - 1)
            print*, "Derivative ", N, " : ", pexp%derivative(x, N), " / ", (v3 - v2) / (2.0d0 * dx)
        else
            print*, "Derivative ", N, " : ", pexp%derivative(x, N)
        endif
        
    enddo
    print*, " ------------- POLYEXP DEMO ENDS -------------"
    print*, ""
    print*, " ------------- POLYFRAC DEMO STARTS -------------"
    ! PolyFrac object : f(x) / g(x), where f and g are poly1 and poly2
    call pfrac%init(num=poly1, denom=poly2)
    print*, "PolyFrac object : "
    call pfrac%print()
    dx = 1e-8
    x = 2.1457
    print*, "Eval PolyFrac at x = ", x, " : ", pfrac%eval(x)

    do N = 0, 4
        if (N > 0) then
            v2 = pfrac%derivative(x - dx, N - 1)
            v3 = pfrac%derivative(x + dx, N - 1)
            print*, "Derivative ", N, " : ", pfrac%derivative(x, N), " / ", (v3 - v2) / (2.0d0 * dx)
        else
            print*, "Derivative ", N, " : ", pfrac%derivative(x, N)
        endif

    enddo
    print*, " ------------- POLYFRAC DEMO ENDS -------------"

end program polynomial_demo
