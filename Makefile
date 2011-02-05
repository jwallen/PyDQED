################################################################################
#
#   Makefile for PyDQED
#
################################################################################

F90=gfortran

CFLAGS=-fPIC -O3

all: DQED cython

cython:
	python setup.py build_ext --build-lib . --build-temp build --pyrex-c-in-temp

install:
	python setup.py install

DQED: libdqed.a

libdqed.a: build/dqed.o
	ar rcs libdqed.a build/dqed.o

build/dqed.o: dqed.f90
	mkdir -p build
	$(F90) $(CFLAGS) -c dqed.f90 -o build/dqed.o

clean: clean-DQED clean-cython
	rm -rf build

clean-DQED:
	rm build/dqed.o libdqed.a

clean-cython:
	python setup.py clean --build-temp build
	rm -f *.so *.pyc

help:
	@echo ""
	@echo "This makefile can be used to build PyDQED and its dependencies."
	@echo ""
	@echo "Typing \`make\` with no arguments will compile DQED to a static library"
	@echo "and compile the PyDQED Python module that provides the Python interface"
	@echo "to DQED."
	@echo ""
	@echo "Typing \`make clean\` will delete all of the intermediate build files,"
	@echo "compiled libraries, and compiled Python modules for DQED and the PyDQED"
	@echo "modules."
	@echo ""
	@echo "Individual dependencies can be specified using \`make <target>\`, where"
	@echo "<target> is one of:"
	@echo ""
	@echo "    DQED     to compile the DQED code"
	@echo "    cython   to compile the PyDQED Python wrapper modules"
	@echo ""

