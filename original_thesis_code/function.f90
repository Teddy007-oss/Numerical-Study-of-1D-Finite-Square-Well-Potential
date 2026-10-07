! Author Emmanuel Osei Boakye 
! Date : 20th May 2024 
! University Of Ghana Department of Physics
! this program helps us plot a graph for the energy transcendental equations to see where the roots lie


	program graphing
	
		implicit none
		real(kind = 8) :: ffo,ffe,fe,fo,E, alpha, beta,f
		real(kind = 8),parameter:: Vo = 10.0, a = 3.0, hbar_sq = 7.6199682, m = 1.0
		integer :: i 
			
		open(1, file = 'graph_table')
		
		E = 0.00
		
		do i = 1, 1000
		E = 0.01 + E
		
		fe = ffe(E)
		
		fo = ffo(E)
		
		print*, E, fe, fo
		write(1,*) E, fe, fo
		
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



!function for the even states!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	function ffe(E)
		implicit none 
		real(kind = 8)::ffe, E, alpha, beta
		real(kind = 8),parameter:: Vo = 10.0, a = 3.0, hbar_sq = 7.6199682, m = 1.0
		
		ffe = alpha(E)*sin(alpha(E)*a) - beta(E)*cos(alpha(E)*a)
	end function
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



!function for the odd states!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	function ffo(E)
		implicit none
		real(kind = 8):: ffo, E, alpha, beta
		real(kind = 8),parameter:: Vo = 10.0, a = 3.0, hbar_sq = 7.6199682, m = 1.0
		
		ffo = alpha(E)*cos(alpha(E)*a) + beta(E)*sin(alpha(E)*a)
	end function 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

