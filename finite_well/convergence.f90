! Author Emmanuel Osei Boakye
! University of Ghana Department of Physics
! Compares the hybrid method with plain bisection on one energy level.
! The user gives an interval [xl, xr] that brackets the root, exactly as in the original thesis codes.
! usage:  echo "V0 a parity xl xr" | ./convergence > convergence.dat      (parity is even or odd)
! output: two blocks of (iteration, E):  block 0 = hybrid,  block 1 = bisection

	program convergence
		use well
		implicit none
		real(kind = 8) :: xl, xr, E
		character(len = 4) :: parity
		integer :: n, nb
		real(kind = 8), parameter :: tol = 1.0d-12

		read(*,*) V0, a, parity, xl, xr

		print*, '# hybrid  (V0 =', V0, ' a =', a, ' ', parity, ' state in [', xl, ',', xr, '])'
		if (parity .eq. 'even') then
			E = hybrid(f_even, xl, xr, tol, n, 6)
			print*
			print*
			print*, '# bisection'
			E = bisection(f_even, xl, xr, tol, nb, 6)
		else
			E = hybrid(f_odd, xl, xr, tol, n, 6)
			print*
			print*
			print*, '# bisection'
			E = bisection(f_odd, xl, xr, tol, nb, 6)
		end if
		print*
		print*
		print*, '# iterations: hybrid =', n, '  bisection =', nb
	end program
