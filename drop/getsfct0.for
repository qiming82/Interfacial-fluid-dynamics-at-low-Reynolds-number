      subroutine getsfct0(Gma,s,ds,r,f,un,ut,N,dt,Pe)
c------------------------------------------------------------------
c  Explicit Euler method solving mass conservation for surfactant
c------------------------------------------------------------------
      Implicit Double Precision (a-h,o-z)
	parameter(neq=1024*2)
	dimension Gma(neq),s(neq),r(neq),dgm(neq),gtmp(neq)
	dimension f(neq),un(neq),ut(neq),dut(neq)
	dimension ftmp1(neq),ftmp2(neq),ds(neq),Qx(neq)
      dimension dftmp1(neq),dftmp2(neq),ddftmp(neq)
c	dimension gtmp(neq),gtmpm(neq),gmam(neq)
      pi  = 3.14159 265358 97932384 D0

      if(Pe.eq.0.0d0)then
	 Dfs = 0.0d0
	else
	 Dfs = 1.0d0/Pe
	endif

c--------------------------------------------------------
c   compute derivatives
c--------------------------------------------------------
      Do i=2,N
	  dgm(i)=(gma(i+1)-gma(i-1))/(s(i+1)-s(i-1))
	  dut(i)=(ut(i)-ut(i-1))/(s(i)-s(i-1))  ! upwinding  
      Enddo
      dgm(1)=0.0d0
	dgm(N+1)=0.0d0
	dut(1)=ut(2)/(s(2)-s(1))  
      dut(N+1)=-ut(N)/(s(N+1)-s(N))

      Do i=1,N+1
	 ftmp1(i)=gma(i)*r(i)
	 ftmp2(i)=r(i)/ds(i)*dgm(i)
	enddo

      do i=2,N
	 dftmp1(i)=(ftmp1(i)-ftmp1(i-1))/(s(i)-s(i-1)) ! upwinding
	 ddftmp(i)=(ftmp2(i+1)-ftmp2(i-1))/(s(i+1)-s(i-1))
	enddo
      dftmp1(1)=0.0d0
	dftmp1(N+1)=gma(N+1)*(-pi)
      ddftmp(1)=ftmp2(2)/(s(2)-s(1))
      ddftmp(N+1)=-dgm(N)/(s(N+1)-s(N))  ! in computational sense
c---------------------------------------------------------
c  Euler advancing
c---------------------------------------------------------

      Do k=2,N
        Qx(k)= -(ut(k)*dftmp1(k)+ftmp1(k)*dut(k)-
     +	  Dfs*ddftmp(k))/r(k)/ds(k) - f(k)*gma(k)*un(k)
	Enddo
	Qx(1)= -gma(1)/ds(1)*dut(1) +Dfs/r(1)/ds(1)*ddftmp(1)
     +      -f(1)*gma(1)*un(1)
	Qx(N+1)= -2.0*gma(N+1)/ds(N+1)*dut(N+1) +
     +	  Dfs*2.0*ddftmp(N+1)/ds(N+1)**2 -f(N+1)*gma(N+1)*un(N+1)


      Do i=1,N+1
	  gtmp(i)=0.0d0
        gtmp(i)=gma(i)+dt*Qx(i)
	  gma(i)=gtmp(i)
      Enddo

	return
	end
