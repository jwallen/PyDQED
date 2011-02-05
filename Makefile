################################################################################
#
#   Makefile for PyDQED
#
################################################################################

F90=gfortran

CFLAGS=-fPIC -O3

all: DQED

DQED: libdqed.a

libdqed.a: build/dqed.o
	ar rcs libdqed.a build/dqed.o

build/dqed.o: dqed.f90
	mkdir -p build
	$(F90) $(CFLAGS) -c dqed.f90 -o build/dqed.o

clean: clean-DQED
	rm -rf build

clean-DQED:
	rm build/dqed.o libdqed.a
