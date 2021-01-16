      subroutine df(r,dz,ddz,dr,ddr,ds,curv)

      Implicit Double Precision (a-h,o-z)
      
	if (r.lt.0.0000000001) then
	  curv = (- ddr*dz + dr*ddz)/ds**3 + ddz/(dr*ds)
      else
	  curv = (- ddr*dz + dr*ddz)/ds**3 + dz/(r*ds)
      end if

	return
	end