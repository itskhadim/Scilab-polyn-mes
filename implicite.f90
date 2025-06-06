module parametre
integer, parameter::N=10,M=5
end module parametre

module implicite
use parametre
contains
real*8 function  f(x,t)
real*8::x,t
f=2*t+x*(1-x)
end function 
subroutine chaleur(A,U,b)
real*8, dimension(N,N)::A,W
real*8, dimension(N+2,M+2)::U
real*8, dimension(N)::b
integer,parameter:: c=1 ,alpha=0 ,  beta=1, gamma=1
real*8::h,k
h=1.0/(N+1)
k=1.0/M
A=0
b=0
do i=1,N
A(i,i)=1/k+(2*c)/h**2
     if (i .lt. N ) then 
        A(i,i+1)=-c/h**2
        A(i+1,i)=-c/h**2
     end if
end do
W=DGEDI(A)
U=0
do j=1,M+2
    U(1,j)=alpha
    U(N+2,j)=beta
end do
do i=1,N+2
U(i,1)=gamma
end do
   do j=2,M+2
      b(1)=f(1*h,j*k)+(c/k)*U(2,j-1)+(c/h**2)*U(1,j)
      b(N)= f((N+1)*h,j*k)+(c/k)*U(N+1,j-1)+(c/h**2)*U(N+1,j)
      do i=2,N
      b(i)=f((i+1)*h,j*k)+(c/k)*U(i+1,j-1)
      end do
end do
  U(2:N+1,j)=matmul(W,b)
end subroutine 
end module implicite


