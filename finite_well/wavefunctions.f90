! Author Emmanuel Osei Boakye
! University of Ghana Department of Physics
! Evaluates the normalised wavefunctions of ALL the bound states (even and odd) on one grid of x,
! so they can be plotted together.
! usage:  echo "V0 a" | ./wavefunctions > wavefunctions.dat
! output columns:  x   V(x)   psi_1(x)   psi_2(x)  ...   (states in order of energy)

	program wavefunctions
		use well
		implicit none
		integer :: nstates, n, i
		real(kind = 8) :: E(max_states), x, L, Vx, psi(max_states)
		character(len = 4) :: parity(max_states)
		integer :: it_hybrid(max_states), it_bisection(max_states)
		integer, parameter :: npoints = 2000

		read(*,*) V0, a

		call find_states(nstates, E, parity, it_hybrid, it_bisection)

		!plot out to where the least bound state has decayed to about exp(-4) of its value at the wall
		L = a + 4.0d0/beta(E(nstates))

		print*, '# V0 =', V0, '  a =', a, '  states =', nstates
		print*, '# x (A)   V(x) (eV)   psi_1 ... psi_n (1/sqrt(A))'
		do i = 0, npoints
			x = -L + 2*L*i/npoints

			if (abs(x) .le. a) then
				Vx = 0.0d0
			else
				Vx = V0
			end if

			do n = 1, nstates
				if (parity(n) .eq. 'even') then
					psi(n) = psi_even(x, E(n))
				else
					psi(n) = psi_odd(x, E(n))
				end if
			end do

			print '(*(es16.7))', x, Vx, psi(1:nstates)
		end do
	end program
