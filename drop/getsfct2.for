      subroutine getsfct2(Gma,s,ds,r,f,un,ut,N,dt,Pe)
c------------------------------------------------------------------
c  finite-volume method solving mass conservation for surfactant
c------------------------------------------------------------------
      Implicit Double Precision (a-h,o-z)
	parameter(neq=1024*2)
	dimension Gma(neq),s(neq),r(neq),dr(neq)
	dimension f(neq),un(neq),ut(neq),dut(neq)
	dimension gtmp(neq),ds(neq),sln(neq)
c	dimension so(neq),dso(neq),ro(neq)
      dimension at(neq),atp(neq),bt(neq),ct(neq),rhs(neq)
      dimension h(neq),RL(neq)
      dimension A1(neq),A2(neq),A3(neq)
	dimension ftmp(neq),dftmp(neq),ddftmp(neq)

      pi  = 3.14159 265358 97932384 D0

      if(Pe.eq.0.0d0)then
	 Dfs = 0.0d0
	else
	 Dfs = 1.0d0/Pe
	endif
      
c---------------
c compute h(i)
c---------------
	do i=1,N
	  h(i) = s(i+1)-s(i)
	enddo

      Do i=2,N
	  dr(i)=(r(i+1)-r(i-1))/(s(i+1)-s(i-1))  
	  dut(i)=(ut(i+1)-ut(i-1))/(s(i+1)-s(i-1))  
	  dftmp(i)=(r(i+1)/ds(i+1)-r(i-1)/ds(i-1))/(s(i+1)-s(i-1))  
      Enddo
	dr(1)=0.0d0
	dr(N+1)=-pi
	dut(1)=ut(2)/(s(2)-s(1))  
      dut(N+1)=-ut(N)/(s(N+1)-s(N))

      Do i=2,N
	 A1(i)=(ut(i)*dr(i)/r(i)+dut(i))/ds(i)+f(i)*un(i)
	 A2(i)=ut(i)/ds(i)-Dfs*dftmp(i)/r(i)/ds(i)
	 A3(i)=-Dfs/ds(i)**2
	Enddo
      A1(1)=dut(1)/ds(1)+f(1)*un(1)
	A1(N+1)=2.0*dut(N+1)/ds(N+1)+f(N+1)*un(N+1)
	dftmp(1)=(r(2)/ds(2)-r(1)/ds(1))/(s(2)-s(1))
      A2(1)=-Dfs/r(1)/ds(1)*dftmp(1)
	A2(N+1)=0.0d0
	A3(1)=-Dfs/ds(1)**2
	A3(N+1)=-Dfs*2.0/ds(N+1)**2

c	=====
c----------------------------------
c  generate the tridiagonal matrix
c----------------------------------
      Do k=2,N
	 ct(k)=(-A2(k)*h(k)+2.0*A3(k))*dt/h(k-1)/(h(k)+h(k-1))
	 at(k)=1.0+dt*A1(k)+dt*A2(k)*(h(k)-h(k-1))/h(k)/h(k-1)-
     +	 2.0*dt*A3(k)/(h(k)*h(k-1))
	 bt(k)=(A2(k)*h(k-1)+2.0*A3(k))*dt/h(k)/(h(k)+h(k-1))
	 rhs(k)=gma(k)
	Enddo
      
	at(1)=1.0+dt*A1(1)-2.0*dt*A3(1)/h(1)**2
	bt(1)=2.0*dt*A3(1)/h(1)**2
      rhs(1)=gma(1)
      at(N+1)=1.0+dt*A1(N+1)-2.0*dt*A3(N+1)/h(N)**2
      ct(N+1)=2.0*dt*A3(N+1)/h(N)**2
      rhs(N+1)=gma(N+1)

c---------------------
c solve N equations
c---------------------

      call thomas 
     +
     +   (N+1
     +   ,at,bt,ct
     +   ,rhs
     +   ,sln
     +   )

       gma(1:N+1) = sln

	return
	end
c---------------------------------
      subroutine getsfct2s(Gma,s,ds,r,f,un,ut,N,dt,Pe,Bi,Pk)
c------------------------------------------------------------------
c  finite-volume method solving mass conservation for surfactant
c------------------------------------------------------------------
      Implicit Double Precision (a-h,o-z)
	parameter(neq=1024)
	dimension Gma(neq),s(neq),r(neq),dr(neq)
	dimension f(neq),un(neq),ut(neq),dut(neq)
	dimension gtmp(neq),ds(neq),sln(neq)
c	dimension so(neq),dso(neq),ro(neq)
      dimension at(neq),atp(neq),bt(neq),ct(neq),rhs(neq)
      dimension h(neq),RL(neq)
      dimension A1(neq),A2(neq),A3(neq)
	dimension ftmp(neq),dftmp(neq),ddftmp(neq)

      pi  = 3.14159 265358 97932384 D0

      if(Pe.eq.0.0d0)then
	 Dfs = 0.0d0
	else
	 Dfs = 1.0d0/Pe
	endif
      
c---------------
c compute h(i)
c---------------
	do i=1,N
	  h(i) = s(i+1)-s(i)
	enddo

      Do i=2,N
	  dr(i)=(r(i+1)-r(i-1))/(s(i+1)-s(i-1))  
	  dut(i)=(ut(i+1)-ut(i-1))/(s(i+1)-s(i-1))  
	  dftmp(i)=(r(i+1)/ds(i+1)-r(i-1)/ds(i-1))/(s(i+1)-s(i-1))  
      Enddo
	dr(1)=0.0d0
	dr(N+1)=-pi
	dut(1)=ut(2)/(s(2)-s(1))  
      dut(N+1)=-ut(N)/(s(N+1)-s(N))

      Do i=2,N
	 A1(i)=(ut(i)*dr(i)/r(i)+dut(i))/ds(i)+f(i)*un(i)
	 A2(i)=ut(i)/ds(i)-Dfs*dftmp(i)/r(i)/ds(i)
	 A3(i)=-Dfs/ds(i)**2
	Enddo
      A1(1)=dut(1)/ds(1)+f(1)*un(1)
	A1(N+1)=2.0*dut(N+1)/ds(N+1)+f(N+1)*un(N+1)
	dftmp(1)=(r(2)/ds(2)-r(1)/ds(1))/(s(2)-s(1))
      A2(1)=-Dfs/r(1)/ds(1)*dftmp(1)
	A2(N+1)=0.0d0
	A3(1)=-Dfs/ds(1)**2
	A3(N+1)=-Dfs*2.0/ds(N+1)**2

c	=====
c----------------------------------
c  generate the tridiagonal matrix
c----------------------------------
      Do k=2,N
	 ct(k)=(-A2(k)*h(k)+2.0*A3(k))*dt/h(k-1)/(h(k)+h(k-1))
	 at(k)=1.0+dt*A1(k)+dt*A2(k)*(h(k)-h(k-1))/h(k)/h(k-1)-
     +	 2.0*dt*A3(k)/(h(k)*h(k-1))
	 bt(k)=(A2(k)*h(k-1)+2.0*A3(k))*dt/h(k)/(h(k)+h(k-1))
	 rhs(k)=gma(k)+dt*Bi*(1.d0+Pk)*(1.d0-gma(k))
	Enddo
      
	at(1)=1.0+dt*A1(1)-2.0*dt*A3(1)/h(1)**2
	bt(1)=2.0*dt*A3(1)/h(1)**2
      rhs(1)=gma(1)+dt*Bi*(1.d0+Pk)*(1.d0-gma(1))
      at(N+1)=1.0+dt*A1(N+1)-2.0*dt*A3(N+1)/h(N)**2
      ct(N+1)=2.0*dt*A3(N+1)/h(N)**2
      rhs(N+1)=gma(N+1)+dt*Bi*(1.d0+Pk)*(1.d0-gma(N+1))

c---------------------
c solve N equations
c---------------------

      call thomas 
     +
     +   (N+1
     +   ,at,bt,ct
     +   ,rhs
     +   ,sln
     +   )

       gma(1:N+1) = sln

	return
	end
