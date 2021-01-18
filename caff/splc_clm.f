      subroutine splc_clm
     +
     +  (N
     +  ,xp,fp
     +  ,slope1
     +  ,slope2
     +  ,a,b,c
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c           C. Pozrikidis
c
c ``Numerical Computation in Science and Engineering''
c
c      Oxford University Press, 1998
c------------------------------------------------

c--------------------------------------
c  Cubic-spline fit of prescribed data 
c  with clamped-end conditions
c
c  SYMBOLS
c  -------
c
c  N  .... number of intervals
c  xp .... x coord. of prescribed data
c  fp .... y coord. of prescribed data
c	
c  a .... polynomial coefficient related to 3rd derivative
c  b .... polynomial coefficient related to 2nd derivative
c  c .... polynomial coefficient related to 1st derivative
c  h .... interval between prescribed data
c
c  at.... diagonal of tridiagonal matrix
c  bt.... superdiagonal of tridiagonal matrix
c  ct.... subdiagonal of tridiagonal matrix
c
c CAPACITY:
c --------
c
c N=512 max
c
c--------------------------------------

      Implicit Double Precision (a-h,o-z)
      parameter(neq=1024)
c      Dimension  xp(513),fp(513),b(513)
c      Dimension   h(512),dy(512),ddy(512)
c      Dimension   a(512),c(512)

c      Dimension at(512),bt(512),ct(512),rhs(512)
c      Dimension        sln(512)

      Dimension  xp(neq),fp(neq),b(neq)
      Dimension   h(neq),dy(neq),ddy(neq)
      Dimension   a(neq),c(neq)

      Dimension at(neq),bt(neq),ct(neq),rhs(neq)
      Dimension        sln(neq)

c--------
c prepare
c--------

      Na = N-1
      N1 = N+1

c-----------------------
c compute intervals h(i)
c-----------------------

      Do i=1,N
        h(i) = xp(i+1)-xp(i)
      End Do 

      h1 = h(1)
      hN = h(N)

c--------------------------------
C Generate the tridiagonal matrix
c--------------------------------
 
      at(1) = 2.0D0*(h(1)+h(2))- 0.50D0*h1         ! clamped end
      bt(1) = h(2)

      Do i=2,N-1
        i1 = i+1
        at(i) = 2.0D0*(h(i)+h(i1))
        bt(i) = h(i1)
        ct(i) = h(i)
      End Do 

      at(Na) = 2.0D0*(h(Na)+h(N)) - 0.50D0*hN         ! clamped end
      ct(Na) = h(Na)

      Do i=1,Na
        i1 = i+1
        i2 = i+2
        rhs(i)  = 3.0D0*( (fp(i2)-fp(i1))/h(i1)
     +                   -(fp(i1)-fp(i) )/h(i)  )
      End Do 

      rhs(1)  = rhs(1)  - 1.5D0*( (fp(2)  -fp(1))/h1 - slope1)  ! clamped end
      rhs(Na) = rhs(Na) + 1.5D0*( (fp(N+1)-fp(N))/hN - slope2)  ! clamped end
 
c-------------------------
c  Solve the N-1 equations
c-------------------------

      call thomas 
     +
     +   (Na
     +   ,at,bt,ct
     +   ,rhs
     +   ,sln
     +   )

      Do i=1,Na
        b(i+1) = sln(i)
      End Do

      b(1)  = -0.5D0*b(2) + 1.5D0*( (fp(2) -fp(1))/h1 - slope1)/h1 ! clamped end
      b(N1) = -0.5D0*b(N) - 1.5D0*( (fp(N1)-fp(N))/hN - slope2)/hN ! clamped end

c---------------------------------
c Compute the coefficients a and c
c---------------------------------

      Do i=1,N
        i1 = i+1
        a(i) = (b(i1)-b(i))/(3.0D0*h(i))
        c(i) = (fp(i1)-fp(i))/h(i) - h(i)*(b(i1)+2.0D0*b(i))/3.0D0
      End Do

c-----
c Done
c-----
      

      Return
      End
