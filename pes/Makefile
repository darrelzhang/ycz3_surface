#Macro definition

FC = mpif90

#MKL_PATH = /hptc_cluster3/application/compiler/intel/mkl/10.2.5.035/lib/em64t
#MKL_PATH = /hptc_cluster3/application/compiler/intel/mkl/10.1.1.019/lib/em64t
#MKL_PATH = /soft/common/intel/composer_xe_2015.0.090/mkl/lib/intel64
MKL_PATH = /soft/common/intel/composerxe/mkl/lib/intel64
FFLAGS =  -extend_source -O2  -traceback -132 
LDFLAGS = -L$(MKL_PATH)  -lmkl_rt
#-lmkl_lapack64
LIB = libarpack_intel64.a

#end of Macro definition
pes_bm_mpi.exe:    modules.o main.o mass.o readinp.o small.o sub.o bound2d.o \
umbrella.o bf.o cg.o refpotx.o potdriv.o potential.o init.o \
wvfun14.o arpack0.o rottrans0.o schprog.o lib.o trans0.o trans_blas.o \
ch4_ni.o

	 $(FC) $(LDFLAGS) modules.o main.o mass.o readinp.o small.o sub.o bound2d.o \
umbrella.o bf.o cg.o refpotx.o potdriv.o potential.o init.o \
wvfun14.o arpack0.o rottrans0.o schprog.o lib.o trans0.o trans_blas.o \
ch4_ni.o  -openmp  $(LIB) -o pes_bm_mpi.exe

lib.o:          lib.f                                                            
		$(FC) -c -O2 lib.f                                               
schprog.o:      schprog.f                                                        
		$(FC) $(FFLAGS) -c  -openmp schprog.f                                
trans0.o:       trans0.f                                                         
		$(FC) $(FFLAGS) -c -openmp trans0.f                                 
trans_blas.o:   trans_blas.f                                                     
		$(FC) $(FFLAGS) -c -openmp trans_blas.f    
ch4_ni.o:	ch4_ni.F
		${FC} ${FFLAGS} -c ch4_ni.F
arpack0.o:	arpack0.f
		${FC} -c -extend_source arpack0.f


clean:
	rm *.o *.mod *.exe -f
