! Author  Emmanuel Osei Boakye
! Date    Friday 26 April 2024
! This program finds the root of a given function by implementing the quadratic approximation methods......

program quadraticapproximation
implicit none


!defining variables
real(kind = 8) :: x0,x1,x2,x3, delta, x, f, a, b, c
real(kind = 8), parameter :: tol = 1.0d-12
integer :: n

!making three guesses 
x0 = 0.5
x1 = 4
x2 = 5
n = 0

do
n = n + 1

! lets evaluate a,b and c (this will be explained in a readme file)
a = ((x1 - x2) * (f(x0) - f(x2)) - (x0 - x2)*(f(x1) - f(x2))) / ((x0 - x1)*(x0 - x2)*(x1 - x2))
!print*, a

b = (((x0 - x2)**2.0)*(f(x1) - f(x2)) - ((x1 - x2)**2.0)*(f(x0) - f(x2))) / ((x0 - x1)*(x0 - x2)*(x1 - x2))
!print*, b

c = f(x2)
!print*, c


if (b .ge. 0.0) then

delta = (- 2.0*c / (b + sqrt(b**2.0 - 4*a*c)))
!print*, delta

x3 = x2 + delta

elseif (b .le. 0.0) then

delta = (2.0*c / (-b + sqrt(b**2.0 - 4*a*c)))
!print*, delta

x3 = x2 + delta

end if

print*, n, x3

x0 = x1
x1 = x2
x2 = x3

! stop once the step is smaller than the tolerance (or after 50 tries, in case it does not converge)
if (abs(delta) .le. tol .or. n .ge. 50) exit
end do

print*, 'root =', x3, ' found in', n, 'iterations'

end program quadraticapproximation



function f(x)
implicit none

real(kind = 8):: f, x
f = 2*x - 6

end function
