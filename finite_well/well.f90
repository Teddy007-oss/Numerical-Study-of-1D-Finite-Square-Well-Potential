! Author : Emmanuel Osei Boakye
! University of Ghana Department of Physics
! Everything about the finite square well  V(x) = 0 for |x| < a,  V(x) = V0 for |x| > a
! Units: energies in eV, lengths in Angstrom, mass in electron masses, hbar^2/m_e = 7.6199682 eV A^2
! V0 and a are set by the program that uses this module (they are read in, so any well can be studied).

module well
	use roots
	implicit none

	real(kind = 8) :: V0 = 10.0d0, a = 3.0d0
	real(kind = 8), parameter :: hbar_sq = 7.6199682d0, m = 1.0d0
	real(kind = 8), parameter :: pi = 3.14159265358979324d0
	integer, parameter :: nscan = 20000           ! number of steps used to scan 0 < E < V0 for sign changes
	integer, parameter :: max_states = 200

contains

	!wave number inside the well
	function alpha(E)
		real(kind = 8) :: alpha, E
		alpha = sqrt(2.0d0*m*E/hbar_sq)
	end function

	!decay constant outside the well
	function beta(E)
		real(kind = 8) :: beta, E
		beta = sqrt(2.0d0*m*(V0 - E)/hbar_sq)
	end function

	!transcendental equation for the even states, beta = alpha tan(alpha a), written without tan
	function f_even(E)
		real(kind = 8) :: f_even, E
		f_even = alpha(E)*sin(alpha(E)*a) - beta(E)*cos(alpha(E)*a)
	end function

	!transcendental equation for the odd states, -beta = alpha cot(alpha a), written without cot
	function f_odd(E)
		real(kind = 8) :: f_odd, E
		f_odd = alpha(E)*cos(alpha(E)*a) + beta(E)*sin(alpha(E)*a)
	end function

	!dimensionless well strength z0 = a sqrt(2 m V0)/hbar.  Number of bound states = floor(2 z0/pi) + 1
	function z0()
		real(kind = 8) :: z0
		z0 = a*sqrt(2.0d0*m*V0/hbar_sq)
	end function


	!normalised even wavefunction:  cos(alpha x) inside,  cos(alpha a) exp(-beta(|x|-a)) outside
	function psi_even(x, E)
		real(kind = 8) :: psi_even, x, E, k, q, norm
		k = alpha(E)
		q = beta(E)
		norm = sqrt(a + sin(2*k*a)/(2*k) + cos(k*a)**2/q)
		if (abs(x) .le. a) then
			psi_even = cos(k*x)/norm
		else
			psi_even = cos(k*a)*exp(-q*(abs(x) - a))/norm
		end if
	end function

	!normalised odd wavefunction:  sin(alpha x) inside,  sign(x) sin(alpha a) exp(-beta(|x|-a)) outside
	function psi_odd(x, E)
		real(kind = 8) :: psi_odd, x, E, k, q, norm
		k = alpha(E)
		q = beta(E)
		norm = sqrt(a - sin(2*k*a)/(2*k) + sin(k*a)**2/q)
		if (abs(x) .le. a) then
			psi_odd = sin(k*x)/norm
		else
			psi_odd = sign(1.0d0, x)*sin(k*a)*exp(-q*(abs(x) - a))/norm
		end if
	end function


	!probability of finding the particle OUTSIDE the well (in the classically forbidden region)
	function p_outside(E, parity)
		real(kind = 8) :: p_outside, E, k, q
		character(len = 4) :: parity
		k = alpha(E)
		q = beta(E)
		if (parity .eq. 'even') then
			p_outside = (cos(k*a)**2/q) / (a + sin(2*k*a)/(2*k) + cos(k*a)**2/q)
		else
			p_outside = (sin(k*a)**2/q) / (a - sin(2*k*a)/(2*k) + sin(k*a)**2/q)
		end if
	end function


	!Find every bound state of the well.
	!Step 1: scan 0 < E < V0 and note every interval where f_even or f_odd changes sign.
	!Step 2: solve for the root in each interval with the hybrid method (and with bisection, to compare).
	!The states come out sorted by energy, so state number n = 1, 2, 3, ...
	subroutine find_states(nstates, E, parity, it_hybrid, it_bisection)
		integer, intent(out) :: nstates
		real(kind = 8), intent(out) :: E(max_states)
		character(len = 4), intent(out) :: parity(max_states)
		integer, intent(out) :: it_hybrid(max_states), it_bisection(max_states)
		real(kind = 8) :: El, Er, dE, tol, Eb
		integer :: i

		dE = V0/nscan
		tol = 1.0d-12*V0
		nstates = 0

		do i = 0, nscan - 1
			El = i*dE
			Er = (i + 1)*dE

			if (f_even(El)*f_even(Er) .lt. 0) then
				nstates = nstates + 1
				E(nstates) = hybrid(f_even, El, Er, tol, it_hybrid(nstates))
				Eb = bisection(f_even, El, Er, tol, it_bisection(nstates))
				parity(nstates) = 'even'
			end if

			if (f_odd(El)*f_odd(Er) .lt. 0) then
				nstates = nstates + 1
				E(nstates) = hybrid(f_odd, El, Er, tol, it_hybrid(nstates))
				Eb = bisection(f_odd, El, Er, tol, it_bisection(nstates))
				parity(nstates) = 'odd'
			end if

			if (nstates .ge. max_states) exit
		end do
	end subroutine

end module well
