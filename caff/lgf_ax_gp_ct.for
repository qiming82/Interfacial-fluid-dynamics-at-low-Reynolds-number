      subroutine lgf_ax_gp_ct
     +
     +   (Iopt,iop
     +   ,X,Y
     +   ,X0,Y0
     +   ,RL,d
     +   ,Nsum,Np
     +   ,G
     +   )

c======================================================

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------
      Iopt = 1
      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c--------
c prepare
c--------
      ssss = y+y0-2.0d0*d
	ssss2 = ssss*ssss
c      subtr = pi4*Y
      RLH   = 0.5D0*RL
	dx = x-x0
	dx0 = dx

c---------------------------
c prepare for fast summation
c using Aitken extrapolation
c---------------------------

      N0 = Nsum
      N1 = N0*Np
      N2 = N1*Np

c-------------
c primary ring
c-------------

      call lgf_ax_fs
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,G
     +   ,GX,GR
     +   )
       G = G - 0.5D0/dsqrt(dx*dx+ssss2)

       G =0.0D0


c------------------------
c sum over periodic rings
c------------------------

      Xsave = X

      Do 1 Ir=1,N2

      X = Xsave+Ir*RL     ! forward ring

      call lgf_ax_fs
     +
     +   (Iopt
     +   ,X,Y,X0,Y0
     +   ,S
     +   ,SX,SR
     +   )
      dx = x-x0
      S = S - 0.5D0/dsqrt(dx*dx+ssss2)
c     +   - 0.25D0/abs(X+RLH)
c      S = S - 0.5D0/abs(X+RLH)


      If   (Ir.eq.N0         ! hold for extrapolation
     +  .or.Ir.eq.N1
     +  .or.Ir.eq.N2) then   

       f = s


      End If

      G = G + S


      X = Xsave-Ir*RL      ! backward ring

      call lgf_ax_fs
     +
     +   (Iopt
     +   ,X,Y,X0,Y0
     +   ,S
     +   ,SX,SR
     +   )
      dx = x-x0
      S = S - 0.5D0/dsqrt(dx*dx+ssss2)
c     +   - 0.25D0/abs(X+RLH)
c      S = S - 0.5D0/abs(X+RLH)

      If   (Ir.eq.N0           ! hold for extrapolation
     +  .or.Ir.eq.N1
     +  .or.Ir.eq.N2) then   

       f  = f  + s


      End If

      G = G + S

      If(Ir.eq.N0) then
        G0  = G - 0.5D0 * f
      Else If(Ir.eq.N1) then
        G1  = G - 0.5D0 * f
      Else If(Ir.eq.N2) then
        G2  = G - 0.5D0 * f
      End If

c---

   1  Continue

      X = Xsave  ! restore

c---
c Aitken extrapolation for the velocity
c---
c	if (dabs(G2-G1).lt.0.000001D0)then
c	G = G2
c      else
      G = (G0*G2-G1**2)/(G2-2.0D0*G1+G0)
c      endif
c	if (dabs(G2-2.0D0*G1+G0).lt.0.00000000000001D0)then
c	    write(*,*) 'G2=', G2
c	    write(*,*) 'G1=', G1
c	    write(*,*) 'G0=', G0
c	    write(*,*) 'num=', G0*G2-G1**2
c	    write(*,*) 'den', G2-2.0D0*G1+G0
c	    write(*,*)'den is zero in lgf_ax_gp_ct'
c	    stop
c	  end if
c      if (X0.lt.0.000000000) then
      if(iop.eq.0)then
      call lgf_ax_fs
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0+RL,Y0
     +   ,Go
     +   ,GXo,GRo
     +   )
	 dx = dx0-RL
c       go = go - 0.5D0/abs(X+RLH)
       Go = Go - 0.5D0/dsqrt(dx*dx+ssss2)
c     +   - 0.25D0/abs(X+RLH)
	   G = G - Go

      end if

c-----
c Done
c-----

 100  format (1x,i3,10(1x,f10.5))
 101  format (10(1x,f12.8))

      Return
      End
