
	subroutine rdnum(ntime,n,rint)
	implicit none
	integer :: ntime,n,rint(n)
	real :: a(n)
	a=0.0
	rint=0
	call random_seed()
	call random_number(a)
	rint = int(a*(ntime-1))

	return 
	end subroutine rdnum


