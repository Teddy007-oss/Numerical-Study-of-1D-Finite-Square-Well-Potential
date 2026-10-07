




!Emmanuel Osei Boakye
!22th May 2024 
!University of Ghana Department of Physics 
!This code uses a hybrid of bisection and quadratic approximation to evaluate the eigen energies of the evenstates for the finite square well.


	program evenstates
		implicit none
		
		real(kind = 8):: xl,xr,a,b,c,x0,x1,x2,x3,f,x, delta, root_tester,alpha,beta
		real(kind = 8), parameter :: tol = 1e-12, Vo = 10.0, z = 3.0, hbar_sq = 7.6199682, m = 1.0 
		integer :: n
		
		!selecting a boundary
		xl= 6.000
		xr= 6.4000

		!choosing 3 guesses for quadratic approximation
		x0 = xl
		x1 = xr
		x2 = (xl + xr)/2.0
		
		
		open(1, file ='evenstates')
		
		delta = 1.0
		n = 0
		
		!the algorithm will begin here!
		Do while (delta .gt. tol)
			! let evaluate a,b and c (this will be explained in a readme file)
			a = ( (x1 - x2) * (f(x0) - f(x2)) - (x0 - x2)*(f(x1) - f(x2)) ) / ( (x0 - x1)*(x0 - x2)*(x1 - x2) )
			b = ( ((x0 - x2)**2.0 )*(f(x1) - f(x2)) - ( (x1 - x2)**2.0 )*(f(x0) - f(x2)) ) / ( (x0 - x1)*(x0 - x2)*(x1 - x2) )
			c = f(x2)
			
			if (b .ge. 0) then
				delta = (- 2.0*c / (b + sqrt(b**2.0 - 4.0*a*c)))
				root_tester = ((x2 - xl)*(b+sqrt(b**2.0 - 4.0*a*c) - 2.0*c)) * ((x2 - xr)*(b + sqrt(b**2.0 - 4.0*a*C) - 2.0*c) )
			
				if (root_tester .le. 0) then
					x3 = x2 + delta
				end if	
			
			elseif (b .le. 0) then
			
				delta = (2.0*c / (-b + sqrt(b**2.0 - 4.0*a*c)))
				root_tester = ((xl - x2)*(-b + sqrt(b**2.0 - 4.0*a*c) - 2.0*c)) * ((xr - x2)*(-b + sqrt(b**2.0 - 4.0*a*C) - 2.0*c))
			
				if (root_tester .le. 0) then
					x3 = x2 + delta
					!print*,'r'
				end if
				
			else 
				x3 = (xl + xr)/ 2.0
				!print*, "BM"
				delta = abs(xr - xl)		
			end if 
			
			
			if (f(xl)*f(x3) .le. 0)  then
				xr = x3		
			else 
				xl = x3		
			end if	
			
			
			n = n + 1
			write(*,*) n, x3
		
			write(1,*) n, x3
			x0 = x1
			x1 = x2
			x2 = x3	
		end do 				
	end program


	!!!!!!! a function declaration for the function alpha
	function alpha(E)
		implicit none
		real(kind = 8):: alpha,E
		real(kind = 8),parameter :: Vo = 10.0, z = 3.0, hbar_sq = 7.6199682, m = 1.0 
		alpha = sqrt((2.0*m*E)/hbar_sq)
	end function

	!!!!!!! a function declaration for the function beta
	function beta(E)
		implicit none
		real(kind= 8):: beta,E
		real(kind = 8),parameter :: Vo = 10.0, z = 3.0, hbar_sq = 7.6199682, m = 1.0
		beta = sqrt((2.0*m*(Vo-E))/hbar_sq)
	end function	

	!the actual function for the even states		
	function f(E)
		implicit none
		real(kind = 8):: f,E,alpha,beta
		real(kind = 8),parameter :: Vo = 10.0, z = 3.0, hbar_sq = 7.6199682, m = 1.0
		f = alpha(E)*sin(alpha(E)*z) - beta(E)*cos(alpha(E)*z)	
	end function
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

