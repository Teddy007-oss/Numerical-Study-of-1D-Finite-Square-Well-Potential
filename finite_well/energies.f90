! Author Emmanuel Osei Boakye
! University of Ghana Department of Physics
! Finds ALL the bound state energies (even and odd) of a finite square well
! with the hybrid of bisection and quadratic approximation.
! usage:  echo "V0 a" | ./energies > energies.dat
! output columns:  n   parity   E(eV)   iterations(hybrid)   iterations(bisection)   V0   a   P(outside the well)

	program energies
		use well
		implicit none
		integer :: nstates, n
		real(kind = 8) :: E(max_states)
		character(len = 4) :: parity(max_states)
		integer :: it_hybrid(max_states), it_bisection(max_states)

		read(*,*) V0, a

		call find_states(nstates, E, parity, it_hybrid, it_bisection)

		print '(a,f12.5,a,f10.5,a,f10.5)', '# V0 =', V0, '   a =', a, '   z0 =', z0()
		print '(a,i4,a,i4)', '# bound states found =', nstates, '   expected floor(2 z0/pi)+1 =', floor(2*z0()/pi) + 1
		print '(a)', '# n  parity         E (eV)        hybrid  bisection        V0          a     P(outside)'
		do n = 1, nstates
			print '(i4,2x,a4,f22.15,2i10,2es14.6,f10.4)', n, parity(n), E(n), it_hybrid(n), it_bisection(n), V0, a, &
				p_outside(E(n), parity(n))
		end do
	end program
