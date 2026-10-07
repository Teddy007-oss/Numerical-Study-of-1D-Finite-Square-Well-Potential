! Author Emmanuel Osei Boakye
! University Of Ghana Department of Physics
! Tabulates the transcendental equations f_even(E) and f_odd(E) for 0 < E < V0,
! so they can be plotted to see where the roots (the energies) lie.
! usage:  echo "V0 a" | ./transcendental > table.dat
! output columns:  E   f_even(E)   f_odd(E)

	program transcendental
		use well
		implicit none
		real(kind = 8) :: E
		integer :: i
		integer, parameter :: npoints = 1000

		read(*,*) V0, a

		print*, '#  V0 =', V0, '  a =', a
		print*, '#  E (eV)      f_even      f_odd'
		do i = 1, npoints
			E = i*V0/npoints
			print*, E, f_even(E), f_odd(E)
		end do
	end program
