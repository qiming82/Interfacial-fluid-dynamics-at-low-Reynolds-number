      subroutine geten
     +
     +   (z,r
     +   ,Az,Bz,Cz
     +   ,Ar,Br,Cr
     +   ,N,d
     +   ,s
     +   ,y
     +   )

      Implicit Double Precision (a-h,o-z)
      parameter(neq=1024)
	dimension z(neq),r(neq),s(neq)
	dimension F(neq,neq),RHS(neq),x(neq),y(neq)
	dimension h1(neq),h2(neq),f1(neq),f2(neq)
	dimension Az(neq),Bz(neq),Cz(neq),Ar(neq),Br(neq),Cr(neq)
	dimension s1(neq),s2(neq)
      dimension ipiv(neq)
	

      do i=1,N+1
	  do k=1,N
           call int_pot_s
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,z(i),r(i),s(i)
     +   ,1,d
     +   ,Az(k),Bz(k),Cz(k)
     +   ,Ar(k),Br(k),Cr(k)
     +   ,f1(k)
     +   )

           call int_pot_s
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,z(i),r(i),s(i)
     +   ,2,d
     +   ,Az(k),Bz(k),Cz(k)
     +   ,Ar(k),Br(k),Cr(k)
     +   ,f2(k)
     +   )

           call int_pot_sc
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,z(i),r(i)
     +   ,1,d
     +   ,Az(k),Bz(k),Cz(k)
     +   ,Ar(k),Br(k),Cr(k)
     +   ,h1(k)
     +   )
           call int_pot_sc
     +
     +   (z(k),r(k)
     +   ,s(k),s(k+1)
     +   ,z(i),r(i)
     +   ,2,d
     +   ,Az(k),Bz(k),Cz(k)
     +   ,Ar(k),Br(k),Cr(k)
     +   ,h2(k)
     +   )

	      if(k.eq.1) then
               F(i,k)=f1(k)+h1(k)
	      else 
	         F(i,k)=f2(k-1)+h2(k-1)+f1(k)+h1(k)
	       end if
	  end do
	   
	   F(i,N+1)=f2(N)+h2(N)

	   RHS(i)=1.0D0
      end do

c       call gmriter(F,RHS,N+1,X,Y)
c     Gauss elimination
c         call gel
c     +
c     + (N+1      ! system size
c     + ,F      ! coefficient matrix
c     + ,RHS    ! right-hand side
c     + ,y      ! solution
c     + ,0
c     + ,1
c     + ,ml,mu
c     + ,det
c     + ,Istop
c     + )

         call dgesv(N+1, 1, F, neq, ipiv, RHS, neq, info)
      y(1:N+1)=RHS(1:N+1)

      return 
	end