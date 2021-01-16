      subroutine sgf_ax_g
     +
     + (Iopt
     + ,X,Y
     + ,X0,Y0
     + ,SXX,SXY
     + ,SYX,SYY
     + ,Iaxis
     + )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c================================================
c  Axisymmetric Green's function of Stokes flow
c  in free space
c
c  Let b be the strength of the point-force ring located at x;
c  then the induced velocity field is:
c
c  ux(x0) = Sxx(x,x0) * bx + Sxs(x,x0) * bs
c  us(x0) = Ssx(x,x0) * bx + Sss(x,x0) * bs
c
c  The kernel of the axisymmetric double-layer potential is:
c
c   Idlpx(x0) = ux * ( Qxxx * vnx + Qxxs * vns)
c             + us * ( Qxsx * vnx + Qxss * vns)
c
c   Idlps(x0) = ux * ( Qsxx * vnx + Qsxs * vns)
c             + us * ( Qssx * vnx + Qsss * vns)
c
c  arguments of Qxxx are (x,x0)
c
c  This is the flow due to a ring distribution of stresslets
c
c
c------------
c
c  Iopt = 1 produces only the Green's function
c  Iopt = 2 produces the Green's function and the stress tensor
c
c================================================

      Implicit Double Precision (a-h,o-z)

      Parameter (eps=0.00000001)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi

c------------------------
c field point on x axis ?
c------------------------

      Iaxis = 0
      If(y0.lt.eps) Iaxis = 1

c-------
c prepare
c-------

      Y2  = Y*Y
      Y02 = Y0*Y0
      YY2 = Y2+Y02
      YYP = Y*Y0
      DY  = Y-Y0
      DY2 = DY*DY

 
      if(Iaxis.eq.0) then   ! off the axis
        Y4  = Y2*Y2
        Y04 = Y02*Y02
c	  if (YYP.lt.eps) then
c	    write(*,*) 'arguement is negative'
c	    write(*,*) 'YYP is ', YYP
c	    write(*,*) 'Y0 is', Y0
c	    write(*,*) 'Y is',  Y
c	  end if
        YYR = dsqrt(YYP)
        YY3 = YYR**3
        YY5 = YYR**5
        SY  = Y+Y0
        SY2 = SY*SY
      end if

      DX   = X-X0
      DX2  = DX*DX
      DX4  = DX2*DX2
      DR2  = DX2+DY2
      DR   = dsqrt(DR2)
      DXYY = DX2+YY2


c------------------------
      If(Iaxis.eq.0) then       !  point x0 off the axis
c------------------------

        FC1  = 4.0D0*YYP
        FC2  = DX2+SY2
        RK2  = FC1/FC2
        RK   = dsqrt(RK2)
        RK3  = RK2*RK
        RK4  = RK2**2
        RK5  = RK4*RK
        RK2P = 1.0D0-RK2

        call ell_int (RK2,F,E)

        RJ10 = 2.0D0*RK*F/YYR
        RJ11 = RK*(DXYY*F-(DX2+SY2)*E)/YY3
        RJ30 = 2.0D0*RK*E/(YYR*DR2)
        RJ31 = RK*(-F+DXYY*E/DR2)/YY3
        RJ32 = RK*(-DXYY*F+(DX4+2.0*DX2*YY2+Y4+Y04)*E/DR2)/YY5


c-----------
        Else               !  point x0 on the axis
c-----------

        DR3  = DR2*DR
        DR5  = DR3*DR2
        RJ10 = PI2/DR
        RJ11 = 0.0D0
        RJ30 = PI2/DR3
        RJ31 = 0.0D0
        RJ32 = PI/DR3


c-------------
        End If
c-------------

c---
c Build the Green's function
c---

      SXX = Y *    (  RJ10 + DX2*RJ30)
      SXY = Y * DX*(Y*RJ30 - Y0 *RJ31)
      SYX = Y * DX*(Y*RJ31 - Y0 *RJ30)
      SYY = Y *    (  RJ11 + YY2*RJ31-YYP*(RJ30+RJ32))
 
c------

c-----
c Done
c-----

      Return
      End
