CC=g++ -O3 -fopenmp 

CFLAGS=-lblas -llapack
LAPACK= -L /software/cmshared/apps/lapack/gcc/64/3.5.0 -L /usr/lib -L /software/cmshared/apps/blas/gcc/current/lib64 -L /software/cmshared/intel/mkl/lib/intel64 -llapack -lblas -lm -lgomp -fopenmp -ffast-math -lpthread
#LAPACK= -L /cm/shared/apps/lapack/open64/64/3.5.0/ -L /usr/lib -L /software/centos7/lib64/ -L /software/apps/lapack/3.8.0/intel/18.0/lib64/ -llapack -lblas -lm -lgomp -fopenmp -ffast-math -lpthread
BOOST= -I /software/apps/boost/1.68.0/gcc/5.5/openmpi/3.1/include/boost/ -I /software/apps/boost/1.68.0/gcc/5.5/openmpi/3.1/include/


OBJS=	./obj/pentropies_of_states.o \
	./obj/read_inputfile.o \
	./obj/unwrap_pstate.o \
	./obj/Simple_Math_Func.o \
	./obj/Matrix_Functions.o \
	./obj/idsClass.o \
	./obj/Basis_States.o \
	./obj/Momentum_Basis.o \
	./obj/global.o \
	./obj/Entropy.o \
	./obj/pClasses.o
	

EXECUTABLE=pentropy

headers=./src/*.h

MES: $(OBJS) $(headers)
	$(CC) -o ${EXECUTABLE} $(OBJS) ${LAPACK} $(CFLAGS) 

$(OBJS): ./obj/%.o : ./src/%.cpp
	$(CC) -o $@ -c $<

#WallClock.o: ./src/WallClock.h

clean:
	rm -f $(OBJS) $(PROG)
	rm -f $(EXECUTABLE)
