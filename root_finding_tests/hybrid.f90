
! Author Emmanuel Osei Boakye
! Date 3rd May 2024 Friday
! A script the execute the hybrid of the bisection and the quadratic approximation method in finding the root f an equation
! Test equation: f(x) = 5x^2 + 12x + 7, which has roots x = -1.0 and x = -1.4

program hybrid
	implicit none

	!declaring variables
	real(kind = 8):: xl,xr,a,b,c,x0,x1,x2,x3,f,x, delta, root_tester, D
	real(kind = 8), parameter :: tol = 1e-12
	integer :: n
	logical :: quadratic_ok

	!selecting an the interval for the program (f must change sign between xl and xr)
	!f(-1.2) = -0.2 and f(-2.8) = 12.6, so only the root -1.4 lies inside
	xl =-1.2
	xr =-2.8

	!selecting three values of x that lies between the interval already selected(thus xl and xr)
	x0 = xl
	x1 = xr
	x2 = (xl + xr)/2.0

	n = 0

	Do
		n = n + 1
		! let evaluate a,b and c of the quadratic p(x) = a(x-x2)^2 + b(x-x2) + c
		a = ( (x1 - x2) * (f(x0) - f(x2)) - (x0 - x2)*(f(x1) - f(x2)) ) / ( (x0 - x1)*(x0 - x2)*(x1 - x2) )
		b = ( ((x0 - x2)**2.0 )*(f(x1) - f(x2)) - ( (x1 - x2)**2.0 )*(f(x0) - f(x2)) ) / ( (x0 - x1)*(x0 - x2)*(x1 - x2) )
		c = f(x2)

		!test if quadratic approximation will work: the quadratic must have a real root
		!and that root x3 must lie between xl and xr, i.e. (x3 - xl)*(x3 - xr) <= 0
		quadratic_ok = .false.
		if (b**2 - 4*a*c .ge. 0) then
			if (b .ge. 0) then
				D = b + sqrt(b**2 - 4*a*c)
				delta = -2.0*c / D
				root_tester = ((x2 - xl)*D - 2*c) * ((x2 - xr)*D - 2*c)
			else
				D = -b + sqrt(b**2 - 4*a*c)
				delta = 2.0*c / D
				root_tester = ((x2 - xl)*D + 2*c) * ((x2 - xr)*D + 2*c)
			end if
			if (root_tester .le. 0) quadratic_ok = .true.
		end if

		if (quadratic_ok) then
			!proceed with quadratic approximation
			x3 = x2 + delta
			print*, n, x3, '   quadratic'
		else
			!proceed with bisection method
			x3 = (xl + xr)/ 2.0
			delta = abs(xr - xl)
			print*, n, x3, '   bisection'
		end if

		!keep the root bracketed
		if (f(xl)*f(x3) .le. 0)  then
			xr = x3
		else
			xl = x3
		end if

		!test for convergence
		if (abs(delta) .le. tol .or. f(x3) .eq. 0) exit

		x0 = x1
		x1 = x2
		x2 = x3
	end do

	print*, 'root =', x3, ' found in', n, 'iterations'

end program hybrid


   function f(x)
		implicit none

		real(kind = 8):: f, x
		f = 5*x**2 + 12*x + 7

   end function
