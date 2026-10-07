! Author : Emmanuel Osei Boakye
! University of Ghana Department of Physics
! The two root finders used in this project, written once so every program can use them:
!   bisection : the plain bisection method                       (tested in root_finding_tests/bisection.f90)
!   hybrid    : the hybrid of bisection and quadratic approximation (tested in root_finding_tests/hybrid.f90)
! Both take a function f and an interval [xl, xr] on which f changes sign.
! They return the root and set n to the number of iterations used.
! If log_unit is given, every iterate x3 is written to that file (used to plot convergence).

module roots
	implicit none

	integer, parameter :: maxit = 500

contains

	function bisection(f, xl_in, xr_in, tol, n, log_unit) result(x3)
		real(kind = 8), external :: f
		real(kind = 8), intent(in) :: xl_in, xr_in, tol
		integer, intent(out) :: n
		integer, intent(in), optional :: log_unit
		real(kind = 8) :: x3, xl, xr

		xl = xl_in
		xr = xr_in
		n = 0

		do
			n = n + 1
			x3 = (xl + xr)/2.0d0
			if (present(log_unit)) write(log_unit,*) n, x3

			if (f(xl)*f(x3) .le. 0) then
				xr = x3
			else
				xl = x3
			end if

			if (abs(xr - xl) .le. tol .or. n .ge. maxit) exit
		end do
	end function


	function hybrid(f, xl_in, xr_in, tol, n, log_unit) result(x3)
		real(kind = 8), external :: f
		real(kind = 8), intent(in) :: xl_in, xr_in, tol
		integer, intent(out) :: n
		integer, intent(in), optional :: log_unit
		real(kind = 8) :: x3, xl, xr, x0, x1, x2, a, b, c, D, delta
		logical :: quadratic_ok

		xl = xl_in
		xr = xr_in

		!choosing 3 guesses for quadratic approximation
		x0 = xl
		x1 = xr
		x2 = (xl + xr)/2.0d0
		n = 0

		do
			n = n + 1
			! a, b and c of the quadratic p(x) = a(x-x2)^2 + b(x-x2) + c through x0, x1, x2   (eqs 2.03-2.05)
			a = ( (x1 - x2)*(f(x0) - f(x2)) - (x0 - x2)*(f(x1) - f(x2)) ) / ( (x0 - x1)*(x0 - x2)*(x1 - x2) )
			b = ( (x0 - x2)**2*(f(x1) - f(x2)) - (x1 - x2)**2*(f(x0) - f(x2)) ) / ( (x0 - x1)*(x0 - x2)*(x1 - x2) )
			c = f(x2)

			!test if quadratic approximation will work: p(x) must have a real root that lies inside [xl, xr]
			quadratic_ok = .false.
			if (b**2 - 4*a*c .ge. 0) then
				if (b .ge. 0) then
					D = b + sqrt(b**2 - 4*a*c)
					if (D .ne. 0) delta = -2*c/D                        ! eq 2.09
				else
					D = -b + sqrt(b**2 - 4*a*c)
					if (D .ne. 0) delta = 2*c/D                         ! eq 2.11
				end if
				if (D .ne. 0) quadratic_ok = ( (x2 + delta - xl)*(x2 + delta - xr) .le. 0 )
			end if

			if (quadratic_ok) then
				x3 = x2 + delta                       !proceed with quadratic approximation
			else
				x3 = (xl + xr)/2.0d0                  !proceed with bisection   (eq 2.12)
				delta = xr - xl
			end if
			if (present(log_unit)) write(log_unit,*) n, x3

			!keep the root bracketed
			if (f(xl)*f(x3) .le. 0) then
				xr = x3
			else
				xl = x3
			end if

			!test for convergence
			if (abs(delta) .le. tol .or. abs(xr - xl) .le. tol .or. f(x3) .eq. 0 .or. n .ge. maxit) exit

			x0 = x1
			x1 = x2
			x2 = x3
		end do
	end function

end module roots
