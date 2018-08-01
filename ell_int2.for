      subroutine ell_int2 (RK2,F,E)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c===========================================
c  EVALUATION OF COMPLETE ELLIPTIC INTEGRALS
c  OF THE FIRST AND SECOND KIND
c
c  F:   first kind
c  E:   second kind
c===========================================

      Implicit Double Precision (A-H,O-Z)

      parameter (acc=0.000000001)
      dimension a(10),b(10),c(10),d(10)
c----------
c constants
c----------

      pi  = 3.14159 265358 D0
c      pih = 0.5D0*pi
      F = 0.0D0
	E = 0.0D0

c---------------------
c compute by polynomial approximation
c---------------------
      a(1)=0.09666344259 D0
      a(2)=0.03590092383D0
      a(3)=0.03742563713D0
      a(4)=0.01451196212D0

      b(1)=0.12498593597D0
      b(2)=0.06880248576D0
      b(3)=0.03328355346D0
      b(4)=0.00441787012D0

      c(1)=0.44325141463D0
      c(2)=0.06260601220D0
      c(3)=0.04757383546D0
      c(4)=0.01736506451D0

      d(1)=0.24998368310D0
      d(2)=0.09200180037D0
      d(3)=0.04069697526D0
      d(4)=0.00526449639D0

      s1=0.0D0
      s2=0.0D0
      do i=1,4
	  s1=s1+a(i)*(1.0D0-Rk2)**i
        s2=s2+b(i)*(1.0D0-Rk2)**i
      end do
      F = dlog(4.0D0) + s1 - dlog(1.0D0 - Rk2)*(0.5D0 + s2) 

      s1=0.0D0
      s2=0.0D0
      do i=1,4
	   s1=s1+c(i)*(1.0D0-Rk2)**i
          s2=s2+d(i)*(1.0D0-Rk2)**i
      end do
      E = 1.0D0 + s1 - dlog(1.0D0 - Rk2)*s2

c-----
c Done
c-----

      Return
      End
