!By LYJ
program main
  implicit none
  !变量
  real(4):: aE,aW,aP,De,Dw,Fe,Fw,Pe,Pw,L,u,dx,f0,fL,ro,gama
  integer:: n
  !已知
  f0 = 1.0
  fL = 0.0
  L = 1.0
  ro = 1.0
  n = 5
  gama = 0.1
  dx = L/n
  !输入速度
  print*,'输入速度（以x增长方向为正）u=，'
  read(*,*)u
 ! print*,'输入网格数n=，'
!read(*,*)n

  De=gama/dx
  Fe=ro*u
  Dw=gama/dx
  Fw=ro*u
  Pe=Fe/De
  Pw=Fw/Dw

  aE=max(-Fe,De-Fe/2,0.0)
  aW=max(Fw,Dw+Fw/2,0.0)
  aP=aE+aW+(Fe-Fw)

  call solve(aE,aW,aP,n,f0,fL,Fe,Fw,De,Dw,Pe,Pw)

  end program main

  subroutine solve(aE,aW,aP,n,f0,fL,Fe,Fw,De,Dw,Pe,Pw)
    implicit none
    real(4)::aW,aE,aP,f0,fL,Fe,Fw,De,Dw,Pe,Pw,A(n,n),b(n),c(n-1),d(n),x(n)
    integer::I,J,n

    do I=1,n
        do J=1,n
            A(I,J)=0
        end do
        b(I)=0
    end do

	if(abs(Pe)<=2)then
		A(1,1)=0.5*Fe+De+2*Dw
		A(1,2)=0.5*Fe-De
		b(1)=(A(1,1)+A(1,2)+(Fw-Fe))*f0
	else if(Pe>2)then
		A(1,1)=Fe
		b(1)=Fw*f0
	else                                                    !(Pe<-2)
		A(1,1)=-Fw
		A(1,2)=Fe
		b(1)=0
	endif

    do I=2,(n-1)
        A(I,I-1)=-aW
        A(I,I)=aP
        A(I,I+1)=-aE
        b(I)=0
    end do

	if(abs(Pw)<=2)then
		A(n,n-1)=-(0.5*Fw+Dw)
		A(n,n)=-0.5*Fe+2*De+Dw
		b(n)=(A(n,n-1)+A(n,n)+(Fw-Fe))*fL
	else if(Pe>2)then
		A(n,n-1)=-Fw
		A(n,n)=Fe
	else                                                    !(Pw<-2)
		A(n,n)=-Fw
		b(n)=-Fe*fL
	endif

    !求解
    c(1)=A(1,2)/A(1,1)
    d(1)=b(1)/A(1,1)
    do I=2,(n-1)
        c(I)=A(I,I+1)/(A(I,I)-c(I-1)*A(I,I-1))
        d(I)=(b(I)-d(I-1)*A(I,I-1))/(A(I,I)-c(I-1)*A(I,I-1))
    end do
    d(n)=(b(n)-d(n-1)*A(n,n-1))/(A(n,n)-c(n-1)*A(n,n-1))
    x(n)=d(n)

    !结果
    do I=1,(n-1)
        x(n-I)=d(n-I)-c(n-I)*x(n-I+1)
    end do

    print*,'混合格式f=',x
    print*,'系数矩阵A=',A(1,2)

end subroutine solve
