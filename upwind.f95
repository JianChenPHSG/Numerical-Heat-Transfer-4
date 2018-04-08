!by CXY
program main
  implicit none
  real:: aE,aW,aP,De,Dw,Fe,Fw,L,u,dx,f0,fL,ro,gamma
  integer::I,J,n
  !ряж╙
  f0 = 1.0
  fL = 0.0
  L = 1.0
  n = 5
  ro = 1.0
  gamma = 0.1
  dx = L/n

  print*,'input u='
  read(*,*)u

  De=gamma/dx
  Fe=ro*u
  Dw=gamma/dx
  Fw=ro*u

  !upwind
  aE=De+max(-Fe,0.0)
  aW=Dw+max(Fw,0.0)
  aP=aE+aW+(Fe-Fw)

  call solve(aE,aW,aP,n,f0,fL,Fe,Fw,De,Dw)

  end program main

  subroutine solve(aE,aW,aP,n,f0,fL,Fe,Fw,De,Dw)
    implicit none
    real::aW,aE,aP,f0,fL,Fe,Fw,De,Dw,A(5,5),b(5),c(4),d(5),x(5)
    integer::I,J,n

    A(5,5)=0

    if(Fe>0)then
		A(1,1)=Fe+De+2*Dw
		A(1,2)=-De
		b(1)=(A(1,1)+A(1,2)+(Fw-Fe))*f0
	else
		A(1,1)=De+2*Dw
		A(1,2)=-De+Fe
		b(1)=(A(1,1)+A(1,2)+(Fw-Fe))*f0
	end if

    do I=2,(n-1)
        A(I,I-1)=-aW
        A(I,I)=aP
        A(I,I+1)=-aE
        b(I)=0
    end do

    if(Fe>0)then
		A(n,n-1)=-(Fw+Dw)
		A(n,n)=2*De+Dw
		b(n)=(A(n,n-1)+A(n,n)+(Fw-Fe))*fL
	else
		A(n,n-1)=-Dw
		A(n,n)=Fw+2*De+Dw
		b(n)=(A(n,n-1)+A(n,n)+(Fw-Fe))*fL
	end if


    c(1)=A(1,2)/A(1,1)
    d(1)=b(1)/A(1,1)
    do I=2,(n-1)
        c(I)=A(I,I+1)/(A(I,I)-c(I-1)*A(I,I-1))
        d(I)=(b(I)-d(I-1)*A(I,I-1))/(A(I,I)-c(I-1)*A(I,I-1))
    end do
    d(n)=(b(n)-d(n-1)*A(n,n-1))/(A(n,n)-c(n-1)*A(n,n-1))
    x(n)=d(n)


    do I=1,(n-1)
        x(n-I)=d(n-I)-c(n-I)*x(n-I+1)
    end do

    print*,'upwind f=',x

end subroutine solve
