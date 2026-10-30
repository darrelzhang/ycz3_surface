	real*8     :: wf_real8(1:231)
	complex*16 :: wf_cmplx(1:231)
	integer    :: i


	open(100,file='1414_real')
	do i=1,231
	   read(100,*)wf_real8(i)
	enddo
	close(100)

	open(100,file='1414_cmplx')
	do i=1,231
	   read(100,*)wf_cmplx(i)
	enddo
	close(100)

	do i=1, 231
	   write(1415,*) real(wf_cmplx(i)/wf_real8(i)), 
     &                   imag(wf_cmplx(i)/wf_real8(i))
	   write(1416,*) ABS(wf_cmplx(i)*CONJG(wf_cmplx(i))), 
     &                   wf_real8(i)**2
	   write(1417,*) wf_cmplx(i)/
     &                   sqrt(ABS(wf_cmplx(i)*CONJG(wf_cmplx(i))))
	enddo

	end
