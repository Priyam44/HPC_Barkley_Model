default: cwprog

libs = -lboost_program_options -fopenmp
prelib = -fopenmp

CC = g++

cwprog: main.o
	$(CC) -o $@ $^ $(libs)

%.o: %.cpp ReactionDiffusion.h
	$(CC) $(prelib) -O1 -c -o $@ $<



#Default Serial Targets for all the Test cases
# np       -> Number of Threads
# filename -> Name of the file with U and V values
.PHONY: test1
test1:
	./cwprog --Nx 101 --Ny 101 --a 0.75 --b 0.06 --eps 50.0 --mu1 5.0 --mu2 0.0 --np 1 --filename "test1.txt"

.PHONY: test2
test2:
	./cwprog --Nx 251 --Ny 251 --a 0.75 --b 0.06 --eps 13.0 --mu1 5.0 --mu2 0.0 --np 1 --filename "test2.txt"

.PHONY: test3
test3:
	./cwprog --Nx 101 --Ny 101 --a 0.5 --b 0.1 --eps 50.0 --mu1 5.0 --mu2 0.0 --np 1 --filename "test3.txt"

.PHONY: test4
test4:
	./cwprog --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --eps 12.5 --mu1 1.0 --mu2 0.01 --np 1 --filename "test4.txt"


#Cleans the executable and compiled .o files
.PHONY: clean
clean:
	rm -f *.o cwprog
