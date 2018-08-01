      subroutine lgf_ax_gp
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,RL
     +   ,Nsum,Np
     +   ,G
     +   )

c--- complementary part of periodic Green's function
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

c      subtr = pi4*Y
      RLH   = 0.5D0*RL

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

c      G = G - subtr/Dabs(X+RLH-RL)
       G =0.0D0

      call lgf_ax_fs
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0+RL,Y0
     +   ,Go
     +   ,GXo,GRo
     +   )
       Go = Go - 0.5D0/abs(X-RL+RLH)

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

      S = S - 0.5D0/abs(X+RLH)


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

      S = S - 0.5D0/abs(X+RLH)

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

      G = (G0*G2-G1**2)/(G2-2.0D0*G1+G0)

	if ((G2-2.0D0*G1+G0).lt.0.000000001D0)then
	    write(*,*)'den is zero in lgf_ax_gp'
	  end if
      if (X0.lt.0.0000000001) then
	   G = G - Go
	end if

c-----
c Done
c-----

 100  format (1x,i3,10(1x,f10.5))
 101  format (10(1x,f12.8))

      Return
      End
