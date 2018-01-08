FC = gfortran
LIBRARIES = -I/usr/include -lfftw3 -lm
DEFINE = -D$(sum)

ifeq ($(sum), )
  sum = fft
endif

ifeq ($(strip $(mode)),debug)
	FLAGS = -O0 -g -fbounds-check
else
	FLAGS = -O3 
endif
FLAGS += -Jmod -fopenmp

all : mhs_finite mhs_infinite

mhs_finite : harmonics.F90 mhs.F90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $(DEFINE) -Dfinite $^ -o $@

mhs_infinite : harmonics.F90 mhs.F90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $(DEFINE) -Dinfinite $^ -o $@

null_check : harmonics.F90 null_check.f90
	$(FC) $(FLAGS) $(MODULES) $(LIBRARIES) $^ -o $@

clean :
	@rm -r mhs mod/*.mod

datatidy :
	@rm hmi*/synmap*.dat

python : harmpy.f90
	f2py -c -m harmpy harmpy.f90 --f90flags='-fopenmp -g' -lgomp

	# gfortran fftwtest.f90 -o fftwtest -I/usr/include -lfftw3
