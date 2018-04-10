!by CXY
program Work1
  implicit none

  real:: aE1,aE2,aE3,aW1,aW2,aW3,aP1,aP2,aP3,De,Dw,Fe,Fw,Pe,Pw,L,u,dx,f0,fL,ro,gama
  integer::I,J,n,k

!输入节点数
  print*,'INPUT NUMBER OF CELLS n='
  read(*,*)n

!输入速度
  print*,'INPUT INLET VELOCITY u='
  read(*,*)u


!选择方法
write(*,*)'输入选择:1.中心差分；2.上风格式；3.混合格式'
read(*,*)k
if(k==1)then
    CALL center(n,u)
else if(k==2)then
    CALL upwind(n,u)
else if(k==3)then
!if(u==2.5)then
    CALL mixed(n,u)
  !end if
end if

  end program Work1


!center
subroutine center(n,u)
  real:: aE1,aW1,aP1
  real::Dw,Fe,Fw,L,u,dx,f0,fL,ro,gama
  real:: A_1(n,n),b_1(n),c_1(n-1),d_1(n),x_1(n)
  integer::I,J,n


  f0 = 1.0
  fL = 0.0
  L = 1.0
  ro = 1.0
  gama = 0.1
  dx = L/n
  De=gama/dx
  Fe=ro*u
  Dw=gama/dx
  Fw=ro*u

  aE1=De-Fe/2
  aW1=Dw+Fw/2
  aP1=aE1+aW1+(Fe-Fw)

    do I=1,n
        do J=1,n
            A_1(I,J)=0
        end do
        b_1(I)=0
    end do

    A_1(1,1)=0.5*Fe+3*De
    A_1(1,2)=-(De-0.5*Fe)
    b_1(1)=(A_1(1,1)+A_1(1,2))*f0

    do I=2,(n-1)
        A_1(I,I-1)=-aW1
        A_1(I,I)=aP1
        A_1(I,I+1)=-aE1
        b_1(I)=0
    end do

    A_1(n,n-1)=-(Dw+0.5*Fw)
    A_1(n,n)=3*Dw-0.5*Fw
    b_1(n)=(2*D-Fw)*fL


    c_1(1)=A_1(1,2)/A_1(1,1)
    d_1(1)=b_1(1)/A_1(1,1)
    do I=2,(n-1)
        c_1(I)=A_1(I,I+1)/(A_1(I,I)-c_1(I-1)*A_1(I,I-1))
        d_1(I)=(b_1(I)-d_1(I-1)*A_1(I,I-1))/(A_1(I,I)-c_1(I-1)*A_1(I,I-1))
    end do
    d_1(n)=(b_1(n)-d_1(n-1)*A_1(n,n-1))/(A_1(n,n)-c_1(n-1)*A_1(n,n-1))
    x_1(n)=d_1(n)


    do I=1,(n-1)
        x_1(n-I)=d_1(n-I)-c_1(n-I)*x_1(n-I+1)
    end do

    print*,'center f=',x_1

    end subroutine


! upwind
subroutine upwind(n,u)
  real:: aE2,aW2,aP2
  real::Dw,Fe,Fw,L,u,dx,f0,fL,ro,gama
  real:: A_2(n,n),b_2(n),c_2(n-1),d_2(n),x_2(n)
  integer::I,J,n

  f0 = 1.0
  fL = 0.0
  L = 1.0
  ro = 1.0
  gama = 0.1
  dx = L/n
  De=gama/dx
  Fe=ro*u
  Dw=gama/dx
  Fw=ro*u

  aE2=De+max(-Fe,0.0)
  aW2=Dw+max(Fw,0.0)
  aP2=aE2+aW2+(Fe-Fw)

    do I=1,n
        do J=1,n
            A_2(I,J)=0
        end do
        b_2(I)=0
    end do


    if(Fe>0)then
		A_2(1,1)=Fe+De+2*Dw
		A_2(1,2)=-De
		b_2(1)=(A_2(1,1)+A_2(1,2)+(Fw-Fe))*f0
	else
		A_2(1,1)=De+2*Dw
		A_2(1,2)=-De+Fe
		b_2(1)=(A_2(1,1)+A_2(1,2)+(Fw-Fe))*f0
	end if

    do I=2,(n-1)
        A_2(I,I-1)=-aW2
        A_2(I,I)=aP2
        A_2(I,I+1)=-aE2
        b_2(I)=0
    end do

    if(Fe>0)then
		A_2(n,n-1)=-(Fw+Dw)
		A_2(n,n)=2*De+Dw
		b_2(n)=(A_2(n,n-1)+A_2(n,n)+(Fw-Fe))*fL
	else
		A_2(n,n-1)=-Dw
		A_2(n,n)=Fw+2*De+Dw
		b_2(n)=(A_2(n,n-1)+A_2(n,n)+(Fw-Fe))*fL
	end if


    c_2(1)=A_2(1,2)/A_2(1,1)
    d_2(1)=b_2(1)/A_2(1,1)
    do I=2,(n-1)
        c_2(I)=A_2(I,I+1)/(A_2(I,I)-c_2(I-1)*A_2(I,I-1))
        d_2(I)=(b_2(I)-d_2(I-1)*A_2(I,I-1))/(A_2(I,I)-c_2(I-1)*A_2(I,I-1))
    end do
    d_2(n)=(b_2(n)-d_2(n-1)*A_2(n,n-1))/(A_2(n,n)-c_2(n-1)*A_2(n,n-1))
    x_2(n)=d_2(n)


    do I=1,(n-1)
        x_2(n-I)=d_2(n-I)-c_2(n-I)*x_2(n-I+1)
    end do

    print*,'upwind f=',x_2

end subroutine


! mixed
subroutine mixed(n,u)
  real:: aE3,aW3,aP3
  real::Pe,Pw,Dw,Fe,Fw,L,u,dx,f0,fL,ro,gama
  real:: A_3(n,n),b_3(n),c_3(n-1),d_3(n),x_3(n)
  integer::I,J,n


  f0 = 1.0
  fL = 0.0
  L = 1.0
  ro = 1.0
  gama = 0.1
  dx = L/n
  De=gama/dx
  Fe=ro*u
  Dw=gama/dx
  Fw=ro*u
  Pe=Fe/De
  Pw=Fw/Dw


  aE3=max(-Fe,De-Fe/2,0.0)
  aW3=max(Fw,Dw+Fw/2,0.0)
  aP3=aE3+aW3+(Fe-Fw)

    do I=1,n
        do J=1,n
            A_3(I,J)=0
        end do
        b_3(I)=0
    end do



	if(abs(Pe)<=2)then
		A_3(1,1)=0.5*Fe+De+2*Dw
		A_3(1,2)=0.5*Fe-De
		b_3(1)=(A_3(1,1)+A_3(1,2)+(Fw-Fe))*f0
	else if(Pe>2)then
		A_3(1,1)=Fe
		b_3(1)=Fw*f0
	else                                                    !(Pe<-2)
		A_3(1,1)=-Fw
		A_3(1,2)=Fe
	endif

    do I=2,(n-1)
        A_3(I,I-1)=-aW3
        A_3(I,I)=aP3
        A_3(I,I+1)=-aE3
        b_3(I)=0
    end do

	if(abs(Pw)<=2)then
		A_3(n,n-1)=-(0.5*Fw+Dw)
		A_3(n,n)=-0.5*Fe+2*De+Dw
		b_3(n)=(A_3(n,n-1)+A_3(n,n)+(Fw-Fe))*fL
	else if(Pe>2)then
		A_3(n,n-1)=-Fw
		A_3(n,n)=Fe
	else                                                    !(Pw<-2)
		A_3(n,n-1)=-Fw
		b_3(n)=-Fe*fL
	endif

    c_3(1)=A_3(1,2)/A_3(1,1)
    d_3(1)=b_3(1)/A_3(1,1)
    do I=2,(n-1)
        c_3(I)=A_3(I,I+1)/(A_3(I,I)-c_3(I-1)*A_3(I,I-1))
        d_3(I)=(b_3(I)-d_3(I-1)*A_3(I,I-1))/(A_3(I,I)-c_3(I-1)*A_3(I,I-1))
    end do
    d_3(n)=(b_3(n)-d_3(n-1)*A_3(n,n-1))/(A_3(n,n)-c_3(n-1)*A_3(n,n-1))
    x_3(n)=d_3(n)


    do I=1,(n-1)
        x_3(n-I)=d_3(n-I)-c_3(n-I)*x_3(n-I+1)
    end do

    print*,'mixed f=',x_3

end subroutine
