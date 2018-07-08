CC=CC -O3 -fopenmp #-DMPI -fopenmp 
#CC=CC  -O3 -

#CFLAGS= -lgfortran -fopenmp -L -larpack -lparpack

BOOST= -I /sw/xe/boost/1.53.0/sles11.3_gnu4.8.2/boost_1_53_0/boost/ -I /sw/xe/boost/1.53.0/sles11.3_gnu4.8.2/boost_1_53_0/
CFLAGS= -lgfortran -L /u/sciteam/hiteshjc/ARPACK -larpack
#CFLAGS= -L /u/sciteam/hiteshjc/ARPACK -larpack

OBJS=	./obj/ED.o \
	./obj/read_input_file.o \
	./obj/Simple_Math_Func.o \
	./obj/Matrix_Functions.o \
	./obj/idsClass.o \
	./obj/Basis_States.o \
	./obj/sort.o \
	./obj/Heisenberg_Hamiltonian.o \
	./obj/Momentum_Basis.o \
	./obj/ARPACK_sector.o \
	./obj/ARPACK_full_eigs.o \
	./obj/global.o \
	./obj/H_v.o \
	./obj/pClasses.o \
	./obj/lanczos.o \
	./obj/Lanczos_sector.o \
	./obj/Lin_Tables.o
	

EXECUTABLE=ED2

headers=./src/*.h

ED: $(OBJS) $(headers)
	$(CC) -o ED2 $(OBJS) $(CFLAGS) 

#ED.o: ./src/ED.cpp ./src/WallClock.h 
#	$(CC) -c ./src/ED.cpp

#WallClock.o : ./src/WallClock.cpp ./src/WallClock.h
#	$(CC) -c ./src/WallClock.cpp



#$(EXECUTABLE): $(OBJS) $(headers)
#	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJS)


#$(EXECUTABLE): $(OBJS) 
#	$(CC) -o $@ $< $(CFLAGS)

$(OBJS): ./obj/%.o : ./src/%.cpp
	$(CC) $(BOOST) -o $@ -c $<

#WallClock.o: ./src/WallClock.h

clean:
	rm -f $(OBJS) $(PROG)
	rm -f $(EXECUTABLE)
