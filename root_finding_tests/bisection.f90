
!Author Emmanuel Osei Boakye
!Date 24th April 2024
! A program to implement bisection method

	program Bisection
		implicit none
		
		!Declaring variables
		real(kind = 8):: xr, xl, f, x3, error, error1  !where x3 happens to be the guessed root.
		integer:: n
		real(kind = 8), parameter :: tol = 1.0d-12

		! lets obtain an interval between which the root must exist and also do a test to confirm
		print*, 'choose a the left boundary'
			read(*,*)xl

		print*, 'choose a the right boundary'
			read(*,*)xr

		if (f(xl)*f(xr) .lt. 0) then
			print*, 'the root lies between that boundary'

		else
			print*, 'the root is not in the boundary hence start over'
			stop
		end if


		x3 = (xl + xr)/ 2.0
		!print*, x3

		error = abs(xr - xl)
		!print*, error
		n = 0

		Do while (error .gt. tol)
			n = n + 1

			if (f(xl)*f(x3) .le. 0)  then
				xr = x3
				
			else 
				xl = x3	
				
			end if	
			error1 = x3
			x3 = (xl + xr)/2.0
			
			error = abs(error1 - x3)
			
			print*, n, x3
		end do

		print*, 'root =', x3, ' found in', n, 'iterations'

	end program Bisection
		
	
	
	
	function f(x)
		implicit none 
		real(kind = 8):: f, x
		
		f = x**2 - 1
	end function
	
	
	
	

