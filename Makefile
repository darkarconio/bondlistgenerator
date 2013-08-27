CC=g++
VPATH=trunk
CFLAGS=-Wall -ansi -I$(VPATH)# -std=c++0x
LDFLAGS=

all: generator

generator: matrix.o point.o matrix3.o atom.o bond.o angle.o fxn.o main.cpp
		$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

fxn.o: fxn.cpp
		$(CC) $(CFLAGS) $(LDFLAGS) -c $< -o $@

atom.o: atom.cpp
		$(CC) $(CFLAGS) -c $< -o $@

matrix.o: matrix.cpp
		$(CC) $(CFLAGS) -c $< -o $@

point.o: point.cpp
		$(CC) $(CFLAGS) -c $< -o $@

matrix3.o: matrix3.cpp
		$(CC) $(CFLAGS) -c $< -o $@

molecule.o: molecule.hpp
		$(CC) $(CFLAGS) -c $< -o $@

bond.o: bond.cpp
		$(CC) $(CFLAGS) -c $< -o $@

angle.o: angle.cpp
		$(CC) $(CFLAGS) -c $< -o $@

fxn.cpp: fxn.hpp point.hpp matrix.hpp 
atom.cpp: atom.hpp point.hpp
point.cpp: point.hpp matrix.hpp
matrix.cpp: matrix.hpp minexcept.hpp
matrix3.cpp: matrix.hpp minexcept.hpp
bond.cpp: bond.hpp atom.hpp
angle.cpp: angle.hpp bond.hpp atom.hpp

clean:
		rm -rf *o generator
