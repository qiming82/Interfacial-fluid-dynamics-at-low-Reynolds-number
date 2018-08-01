      subroutine thomas 
     +
     +   (N
     +   ,A,B,C
     +   ,S
     +   ,X)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------
c Thomas algorithm for tridiagonal systems
c
c LEGEND:
c ------
c
c N: Number of equations
c X: Solution
c 
c------------------------

      Implicit Double Precision (a-h,o-z)
      parameter(neq=1024)
c      Dimension a(512),b(512),c(512),s(512),x(512)
c      Dimension d(512),y(512)
      Dimension a(neq),b(neq),c(neq),s(neq),x(neq)
      Dimension d(neq),y(neq)

      Parameter (tol=0.0000001D0)

c--------
c Prepare
c--------

      Na = N-1
      
	C(1)=0.0D0
	B(N)=0.0D0
c----------
c reduction
c----------

      D(1) = B(1)/A(1)
      Y(1) = S(1)/A(1)

      Do I=1,Na
       I1 = I+1
       Den   = A(I1)-C(I1)*D(I)
       D(I1) = B(I1)/Den
       Y(I1) =(S(I1)-C(I1)*Y(I))/Den
      End Do

c------------------
c Back Substitution
c------------------

      X(N) = Y(N)
      Do I=Na,1,-1
        X(I)= Y(I)-D(I)*X(I+1)
      End Do

c-------------
c Verification
c-------------

      Res = s(1)-a(1)*X(1)-b(1)*X(2)
      If(abs(Res).gt.tol) write (6,*) " thomas: alarm"

      Do I=2,N-1
        Res = s(I)-c(I)*X(I-1)-a(I)*X(I)-b(I)*X(I+1)
        If(abs(Res).gt.tol) write (6,*) " thomas: alarm"
      End Do

      Res = s(N)-c(N)*X(N-1)-a(N)*X(N)
      If(abs(Res).gt.tol) write (6,*) " thomas: alarm"

c-----
c Done
c-----

 100  Format (1x,f15.10)

      Return
      End
