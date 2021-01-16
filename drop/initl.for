      subroutine initl(z,r,Gama,N1)
      Implicit Double Precision (a-h,o-z)
      parameter(nk=1024*2)
	dimension z(nk),r(nk),Gama(nk)
      open (1,file='xy0.dat',status='old')
	open (2,file='gma0.dat',status='old')
	do i=N1,1,-1
	  read(1,*) z(i), r(i)
	  read(2,*) z(i), Gama(i)
	enddo

	return
	end 
