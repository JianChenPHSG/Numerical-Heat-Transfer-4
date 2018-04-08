
!by HXL
program main
  implicit none

  real:: aE,aW,aP,De,Dw,Fe,Fw,L,u,dx,f0,fL,ro,gamma
  integer::I,J,n

  f0 = 1.0
  fL = 0.0
  L = 1.0
  n = 5
  ro = 1.0
  gamma = 0.1
  dx = L/n

  print*,'ÊäÈëËÙ¶Èu='
  read(*,*)u

  De=gamma/dx
  Fe=ro*u
  Dw=gamma/dx
  Fw=ro*u

  aE=De-Fe/2
  aW=Dw+Fw/2
  aP=aE+aW+(Fe-Fw)

  call solve(aE,aW,aP,n,f0,fL)

  end program main

  subroutine solve(aE,aW,aP,n,f0,fL)
    implicit none
    real::aW,aE,aP,f0,fL,A(5,5),b(5),c(4),d(5),x(5)
    integer::I,J,n

    A(5,5)=0

    A(1,1)=0.5*Fe+3*De
    A(1,2)=-(De-0.5*Fe)
    b(1)=(2*De+Fe)*f0

    do I=2,(n-1)
        A(I,I-1)=-aW
        A(I,I)=aP
        A(I,I+1)=-aE
        b(I)=0
    end do

    A(n,n-1)=-(Dw+0.5*Fw)
    A(n,n)=3*Dw-0.5*Fw
    b(n)=(2*D-Fw)*fL


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

    print*,'ÖÐÐÄ²î·Ö·¨f=',x

end subroutine solve
