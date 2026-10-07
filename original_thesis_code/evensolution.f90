

!Name : Emmnanuel Osei Boakye
!Date : May 23rd 2024
! Department of Physics university of ghana
! this script contains codes the evaluates the wave function at points of x to be used to plot the wave functions for even states

	program evensolution
		implicit none
		real(kind = 8)::x,psy1,psy2,psy3,sol,E, alpha,beta, B,sol1,sol2,E1,E2
		real(kind = 8),parameter:: Vo = 10.0, a = 3.0, hbar_sq = 7.6199682, m = 1.0, C = 1.0
		
		x = -5.0
		E1 = 0.71545251406668120
		E2 = 6.14859905397473130
		
		!B = (C*exp(-beta(E)*a)/cos(alpha(E)*a))
		
		open(1, file = 'even_solution')
		
		Do while(x .le. 5)
			x = x + 0.001
			
			if (x .ge. -5.0 .and. x .le. -3.0) then
				sol1 = psy1(x,E1)
				sol2 = psy1(x, E2)
				
			else if (x .ge. -3.0 .and. x .le. 3.0) then
				sol1 = psy2(x,E1)
				sol2 = psy2(x, E2)
				
			else
			 	sol1 = psy3(x,E1)	
			 	sol2 = psy3(x,E2)	
			end if	
			
			print*, x,sol1,sol2
		        write(1,*) x,sol1,sol2
		end do
	end program


!!function for alpha!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	function alpha(E)
		implicit none
		real(kind = 8):: alpha, E
		real(kind = 8),parameter:: Vo = 10.0, a = 3.0, hbar_sq = 7.6199682, m = 1.0
		
		alpha = sqrt((2.0*m*E) / hbar_sq)
	end function
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX




!function for beta!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	function beta(E)
		implicit none 
		real(kind = 8):: beta, E
		real(kind = 8),parameter:: Vo = 10.0, a = 3.0, hbar_sq = 7.6199682, m = 1.0
		
		beta = sqrt((2.0*m*(Vo-E))/ hbar_sq)	
	end function
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx




!the equation for the wave function at position i. here the way decays as x approaches negative infinity.	
	function psy1(x,E)
		implicit none
		real(kind = 8):: psy1, x,E,beta,C
		
		C = 1.0
		
		psy1 = C*exp(beta(E)*x)
	end function
	
	
	
		
! equation for the wave function for inside the well.	
	function psy2(x,E)
		implicit none
		real(kind = 8):: psy2,x,E,alpha, B,C,beta,Q
		real(kind = 8),parameter::  a = 3.0
		!Q = E !0.71545251406668120
		Q = E
		C = 1.0
		
		B = (C*exp(-beta(Q)*a)/cos(alpha(Q)*a))
		
		psy2 = B*cos(alpha(E)*x)
	end function
	
	
	
	
	
!the equation for the wave function at position iii. here the way decays as x approaches infinity.	
	function psy3(x,E)
	 	implicit none
	 	real(kind = 8):: psy3,x,E,beta,C
	 	
	 	C=1.0
	 	
	 	psy3 = C*exp(-beta(E)*x)
	end function
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX.	
	

