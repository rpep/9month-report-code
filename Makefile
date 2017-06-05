old: clean
	${CC} -g -c fmm.c -fopenmp
	python setup.py build_ext --inplace
	rm -f fmmpy.c fmmpydecl.c
	rm -rf build
	rm -rf __pycache__
	rm -f *.o

clean:
	rm -rf build
	rm -rf __pycache__
	rm -f fmmpy.c
	rm -rf *.so

.PHONY: test

test:
	python -c 'import fmmpy'

