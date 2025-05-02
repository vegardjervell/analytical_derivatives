! Author: Vegard G. Jervell
! Module for dealing with generic functions of the form
!       t(x) = f(x) * exp[g(x)]
! where f and g are "polynomials" that can contain both negative and positive integer powers of x
! Specifically: We want analytical derivatives up to arbitrary order for these so-called "t-functions".

! Note: This module is a port of parts of the "Factorial" C++ module in use in the KineticGas library
!       Notably, the core "Factorial" classes have not been ported, so be aware that you could get
!       overflow issues if going to too high derivative orders with this module.
!       For scripting purposes, the Factorial C++ module in KineticGas has been exposed to Python.
! Note: Extending this module to deal with non-integer exponents should be pretty straight forward,
!       I just haven't bothered to do it yet as I haven't had the need for them. Let me know if 
!       there's a pressing need for non-integer exponents and I can probably throw it together.
! Note: I am not proficient in Fortran, so this hacky mess could probably benefit from someone 
!       competent having a look at the details.

! Working principle: 
!   "Polynomials" (that can contain both positive and negative integer powers) are represented by the `Polynomial` class
!       - These simply know their coefficients, as well as the start- stop- and step of their exponents. 
!       - Coefficients can be any real number.
!
!   "t-functions" are represented by the `PolyExp` class.
!       - These just hold two `Polynomial` objects: The prefactor and the exponent
! 
!   In order to compute higher order derivatives, we need binomial coefficients and factorials. Note that these 
!   can cause some overflow issues when going to higher order derivatives.
! 
!   For details on the math and algorithms, see the memo in the KineticGas repo: https://thermotools.github.io/KineticGas/memo/index.html

module polynomials
implicit none
public

    ! The Polynomial type represents Laurent polynomials with integer exponents and constant spacing
    ! between exponents (k_step)
    type :: Polynomial
        real, allocatable :: coeff(:)        ! Array of coefficients
        integer :: k_min, k_max, k_step      ! Range and step for exponents (inclusive both k_min and k_max)
        integer :: max_k_negative            ! (k_max < 0) ? k_max : -1
    contains
        procedure :: init => init_polynomial
        procedure :: eval => eval_polynomial
        procedure :: derivative => derivative_polynomial
        procedure :: print => print_polynomial
    end type Polynomial

    ! The PolyExp type represents functions of the form f(x) * exp[g(x)] where f and g are Laurent polynomials,
    ! with f being the prefactor (pref) and g being the exponent (expo)
    type :: PolyExp
        type(Polynomial) :: pref
        type(Polynomial) :: expo
    contains
        procedure :: init => init_polyexp
        procedure :: eval => eval_polyexp
        procedure :: derivative => derivative_polyexp
        procedure :: get_Gk => polyexp_get_Gk
        procedure :: print => print_polyexp
    end type PolyExp

    ! The PolyFrac type represents Laurent polynomial fractions, f(x) / g(x), with f and g being Laurent polynomials.
    ! num is the numerator, denom is the denominator.
    type :: PolyFrac
      type(Polynomial) :: num
      type(Polynomial) :: denom
    contains
      procedure :: init       => init_polyfrac
      procedure :: eval       => eval_polyfrac
      procedure :: derivative => derivative_polyfrac
      procedure :: get_Hk     => polyfrac_get_Hk
      procedure :: print      => print_polyfrac
    end type PolyFrac

contains

!-----------------------------------------------------------------------------
!> Find all integer partitions of N_in, with largest value smaller than or equal to m_in
!> Note: Passing m_in < 0 (default) is equivalent to m_in = N_in 
!> \author VGJ, 2025-01-27
!-----------------------------------------------------------------------------
recursive subroutine get_partitions(N_in, m_in, partitions, partition_sizes)
    implicit none
    integer, intent(in) :: N_in ! Number to partition
    integer, intent(in) :: m_in ! max value (inclusive) in the partition
    integer, allocatable, dimension(:, :), intent(out) :: partitions   ! Organised as partitions(i, :) is partition nr. i
    integer, allocatable, dimension(:), intent(out) :: partition_sizes ! partition_sizes(i) is the number of values in partitions(i, :)

    integer :: i, k, m, n, min_k, num_partitions, num_subparts, sub_index
    integer, allocatable, dimension(:, :) :: sub_partitions, tmp_partitions
    integer, allocatable, dimension(:) :: sub_partition_sizes, tmp_partition_sizes

    m = m_in
    if (m < 0) m = N_in

    if (N_in == 0) then
        allocate(partitions(1, 0))   ! One partition with 0 elements
        allocate(partition_sizes(1))
        partition_sizes(1) = 0
        return
    endif

    if (N_in == 1) then
        allocate(partitions(1, 1))   ! One partition with one element (1)
        partitions(1, 1) = 1
        allocate(partition_sizes(1))
        partition_sizes(1) = 1
        return
    endif

    allocate(partitions(0, N_in))
    allocate(partition_sizes(0))
    num_partitions = 0
    min_k = max(0, N_in - m)

    do k = min_k, N_in - 1
        n = N_in - k
        call get_partitions(k, n, sub_partitions, sub_partition_sizes) ! Recursively get partitions of k
        num_subparts = size(sub_partition_sizes)

        ! Expand arrays and append sub-partitions
        allocate(tmp_partitions(num_partitions + num_subparts, N_in))
        tmp_partitions = 0
        if (num_partitions > 0) tmp_partitions(1:num_partitions, :) = partitions
        deallocate(partitions)
        partitions = tmp_partitions
        deallocate(tmp_partitions)

        allocate(tmp_partition_sizes(num_partitions + num_subparts))
        tmp_partition_sizes = 0
        if (num_partitions > 0) tmp_partition_sizes(1:num_partitions) = partition_sizes
        deallocate(partition_sizes)
        partition_sizes = tmp_partition_sizes
        deallocate(tmp_partition_sizes)

        do sub_index = 1, num_subparts
            partitions(num_partitions + sub_index, 1) = n
            partitions(num_partitions + sub_index, 2:size(sub_partitions, 2) + 1) = sub_partitions(sub_index, :)
            partition_sizes(num_partitions + sub_index) = sub_partition_sizes(sub_index) + 1
        enddo

        num_partitions = num_partitions + num_subparts
        deallocate(sub_partitions)
        deallocate(sub_partition_sizes)
    enddo
end subroutine get_partitions

!-----------------------------------------------------------------------------
!> Calculate the multiplicity of a partition
!> Used together with get_partitions, as
!>      call get_partitions(N, m, partitions, sizes)
!>      Z = partition_multiplicity(partitions(1), sizes(1)) ! Get multiplicity of the first partition
!> \author VGJ, 2025-01-27
!-----------------------------------------------------------------------------
function partition_multiplicity(partition_in, psize) result(multiplicity)
    integer, intent(in), dimension(:) :: partition_in
    integer, intent(in) :: psize
    integer, dimension(psize) :: partition
    real :: multiplicity

    integer :: n, p, count, i, j, unique_size, freq
    integer, allocatable, dimension(:) :: unique_vals
    real :: num, denom

    partition = partition_in(1:psize)
    n = sum(partition)
    num = factorial(n)
    denom = 1
    do i = 1, psize
        denom = denom * factorial(partition(i))
    enddo

    ! Count the frequency of each unique value
    allocate(unique_vals(size(partition)))
    unique_size = 0

    do i = 1, psize
        if (all(partition(i) /= unique_vals(:unique_size))) then
            unique_size = unique_size + 1
            unique_vals(unique_size) = partition(i)
            freq = count(partition == partition(i))
            denom = denom * factorial(freq)
        endif
    enddo

    multiplicity = num / denom
    deallocate(unique_vals)
end function partition_multiplicity

! Beware of overflow!
function factorial(n) result(fact)
    integer, intent(in) :: n
    real :: fact
    integer :: i

    fact = 1
    do i = 2, n
        fact = fact * i
    end do
end function factorial

! Binomial coefficient
function binom(n, k) result(bcoeff)
    implicit none

    integer, intent(in) :: n, k
    integer :: bcoeff
    integer :: i
    integer :: num, denom

    ! Handle special cases
    if (k < 0 .or. k > n) then
        bcoeff = 0
        return
    end if

    if (k == 0 .or. k == n) then
        bcoeff = 1
        return
    end if

    ! Compute the binomial coefficient iteratively
    num = 1
    denom = 1

    do i = 1, k
        num = num * (n - i + 1)
        denom = denom * i
    end do

    bcoeff = num / denom
end function binom

! Partial factorial (n! / k!)
function partialfactorial(start, stop) result(fac)
    integer, intent(in) :: start, stop
    real :: fac
    integer :: i

    fac = 1.0
    do i = start + 1, stop
        fac = fac * i
    enddo
end function partialfactorial

!-----------------------------------------------------------------------------
!> Initialize generalized polynomial
!> f(x) = \sum C_k x^{\nu_k}, 
!> where \nu_k = {k_min, k_min + k_step, ..., k_max - k_step, k_max}
!> \author VGJ, 2025-01-27
!-----------------------------------------------------------------------------
subroutine init_polynomial(self, k_min, k_max, k_step, coeff)
    class(Polynomial), intent(out) :: self
    integer, intent(in) :: k_min, k_max, k_step
    real, intent(in), dimension(:) :: coeff

    self%k_min = k_min
    self%k_max = k_max
    self%k_step = k_step

    self%max_k_negative = -1
    if (k_max < 0) then
        self%max_k_negative = k_max
    endif

    allocate(self%coeff(size(coeff)))
    self%coeff = coeff
end subroutine init_polynomial

function eval_polynomial(self, x) result(val)
    class(Polynomial), intent(in) :: self
    real, intent(in) :: x                   ! Point at which to evaluate
    real :: val                           ! Return val
    integer :: i, k                         ! i = term index, k = exponent

    val = 0.0
    do i = 1, size(self%coeff)
        k = self%k_min + (i - 1) * self%k_step
        val = val + self%coeff(i) * x**k
    enddo
end function eval_polynomial

function derivative_polynomial(self, x, n) result(val)
    class(Polynomial), intent(in) :: self
    real, intent(in) :: x           ! The position at which to evaluate
    integer, intent(in) :: n        ! The derivative order (0 reduces to function evaluation)
    real :: val                   ! Return val
    integer :: C_idx, k                 ! k = exponent, i = term index
    integer prefactor               ! For negative exponents

    prefactor = -1
    if (mod(n, 2) == 0) then
        prefactor = 1
    endif

    val = 0
    C_idx = 1
    k = self%k_min

    do while (k <= self%max_k_negative) 
        val = val + prefactor * partialfactorial(- k - 1, n - k - 1) * self%coeff(C_idx) * (x ** (k - n));
        k = k + self%k_step
        C_idx = C_idx + 1
    enddo

    if (self%k_max < n) return

    do while (k <= self%k_max)
        if (k < n) then
            C_idx = C_idx + 1
            k = k + self%k_step
            continue
        endif
        val = val + partialfactorial(k - n, k) * self%coeff(C_idx) * (x ** (k - n))
        k = k + self%k_step
        C_idx = C_idx + 1
    enddo
end function derivative_polynomial

subroutine print_polynomial(self)
    class(Polynomial), intent(in) :: self
    integer :: i, k

    write(*, "(A)", advance="no")
    do i = 1, size(self%coeff)
        k = self%k_min + (i - 1) * self%k_step
        if (i > 1) write(*, "(A)", advance="no") " + "
        write(*, "(F6.2,A,I0)", advance="no") self%coeff(i), "*x^", k
    enddo
end subroutine print_polynomial

!-----------------------------------------------------------------------------
!> Initialize "t-function",
!>      t(x) = f(x) * exp[g(x)],
!> where f and g are `Polynomial` objects.
!> \author VGJ, 2025-01-27
!-----------------------------------------------------------------------------
subroutine init_polyexp(self, pref, expo)
    class(PolyExp), intent(out) :: self
    type(Polynomial), intent(in) :: pref, expo

    self%pref = pref
    self%expo = expo
end subroutine init_polyexp

function eval_polyexp(self, x) result(val)
    class(PolyExp), intent(in) :: self
    real, intent(in) :: x
    real :: val

    val = self%pref%eval(x) * exp(self%expo%eval(x))
end function eval_polyexp

function derivative_polyexp(self, x, n) result(val)
    class(PolyExp), intent(in) :: self
    real, intent(in) :: x
    integer, intent(in) :: n
    real :: val
    integer :: k
    real :: df(n + 1), dg(n + 1)

    if ( (self%expo%k_max == 0) .and. (self%expo%k_min == 0) ) then
        val = self%pref%derivative(x, n) * exp(self%expo%coeff(1))
        return
    endif
    
    do k = 0, n
        df(k + 1) = self%pref%derivative(x, k)
        dg(k + 1) = self%expo%derivative(x, k)
    end do

    val = 0
    do k = 0, n
        val = val + binom(n, k) * self%get_Gk(x, k, dg) * df(n - k + 1)
    end do

    val = val * exp(dg(1))
    
end function derivative_polyexp

!-----------------------------------------------------------------------------
!> Recursive factors G_k in the memo: https://thermotools.github.io/KineticGas/memo/index.html
!> \author VGJ, 2025-01-27
!-----------------------------------------------------------------------------
function polyexp_get_Gk(self, x, k, dg) result(Gk)
    class(PolyExp), intent(in) :: self
    real, intent(in) :: x
    integer, intent(in) :: k
    real, dimension(:), intent(in) :: dg
    real :: Gk
    integer, allocatable, dimension(:, :) :: partitions
    integer, allocatable, dimension(:) :: partition_sizes
    integer :: i, l
    real :: tmp, Z
    integer :: max_dg_order

    Gk = 1
    if (k == 0) return
    
    max_dg_order = k
    if ((self%expo%k_min >= 0) .and. (self%expo%k_max < k)) then
        max_dg_order = self%expo%k_max 
    endif

    call get_partitions(k, max_dg_order, partitions, partition_sizes)

    Gk = 0
    do i = 1, size(partitions, 1)
        tmp = 1
        do l = 1, partition_sizes(i)
            tmp = tmp * dg(partitions(i, l) + 1)
        enddo
        if (tmp /= 0) then
            Z = partition_multiplicity(partitions(i, :), partition_sizes(i))
            Gk = Gk + Z * tmp
        endif
    enddo
end function polyexp_get_Gk

subroutine print_polyexp(self)
    class(PolyExp), intent(in) :: self

    write(*, "(A)", advance="no") "[ "
    call self%pref%print()
    write(*, "(A)", advance="no") " ] * exp[ "
    call self%expo%print()
    write(*, "(A)", advance="no") " ]"
end subroutine print_polyexp

!-----------------------------------------------------------------------------
!> Initialize Laurent polynomial fraction,
!>      p(x) = f(x) / g(x),
!> where f and g are `Polynomial` objects.
!> \author VGJ, 2025-05-02
!-----------------------------------------------------------------------------
subroutine init_polyfrac(self, num, denom)
    class(PolyFrac), intent(out) :: self
    type(Polynomial), intent(in) :: num, denom

    self%num   = num
    self%denom = denom
end subroutine init_polyfrac

!-----------------------------------------------------------------------------
!> Evaluate Laurent polynomial fraction at x
!> \author VGJ, 2025-05-02
!-----------------------------------------------------------------------------
function eval_polyfrac(self, x) result(val)
    class(PolyFrac), intent(in) :: self
    real, intent(in)            :: x
    real                        :: val

    val = self%num%eval(x) / self%denom%eval(x)
end function eval_polyfrac

!-----------------------------------------------------------------------------
!> Evaluate polyfrac n'th derivative at x
!> \author VGJ, 2025-05-02
!-----------------------------------------------------------------------------
function derivative_polyfrac(self, x, n) result(val)
    class(PolyFrac), intent(in) :: self
    real, intent(in)            :: x
    integer, intent(in)         :: n
    real                        :: val
    integer                     :: k, k_start
    real                        :: binom_val
    real                        :: df(n + 1), dg(n + 1)

    do k = 0, n
      df(k + 1) = self%num%derivative(x, k)
      dg(k + 1) = self%denom%derivative(x, k)
    end do

    val = 0.0d0
    do k = 0, n
      val = val + binom(n, k) * self%get_Hk(x, k, dg) * df(n - k + 1)
    end do

    val = val / dg(1)
end function derivative_polyfrac

!-----------------------------------------------------------------------------
!> Recursive factors H_k in the memo: https://thermotools.github.io/KineticGas/memo/index.html
!> \author VGJ, 2025-05-02
!-----------------------------------------------------------------------------
function polyfrac_get_Hk(self, x, k, dg) result(Hk)
    class(PolyFrac), intent(in)    :: self
    real, intent(in)               :: x
    integer, intent(in)            :: k
    real, dimension(:), intent(in) :: dg
    real                           :: Hk, tmp, Z
    integer                        :: i, l, max_dg_order, psize
    integer, allocatable           :: partitions(:,:), partition_sizes(:)

    if (self%denom%k_min >= 0 .and. self%denom%k_max < k) then
      max_dg_order = self%denom%k_max
    else
      max_dg_order = k
    end if

    call get_partitions(k, max_dg_order, partitions, partition_sizes)
    Hk = 0.0

    do i = 1, size(partition_sizes)
      tmp = 1.0
      do l = 1, partition_sizes(i)
        if (partitions(i ,l) > max_dg_order) exit
        tmp = tmp * dg(partitions(i, l) + 1)
      end do
      if (tmp /= 0.0) then
        Z = partition_multiplicity(partitions(i,:), partition_sizes(i))
        Hk = Hk + Z * tmp * factorial(partition_sizes(i)) * ( (-1.0 / dg(1))**partition_sizes(i) )
      end if
    end do

    deallocate(partitions, partition_sizes)
end function polyfrac_get_Hk

subroutine print_polyfrac(self)
    class(PolyFrac), intent(in) :: self

    write(*, "(A)", advance="no") "( "
    call self%num%print()
    write(*, "(A)", advance="no") " ) / ( "
    call self%denom%print()
    write(*, "(A)")            " )"
end subroutine print_polyfrac

end module polynomials
