all:
	g++ -std=c++11 main.cpp BinaryTree.cpp Datum.cpp functions.cpp -o Poly -O2  -DARMA_USE_HDF5 -larmadillo `pkg-config --cflags --libs hdf5`
	
	
#-DARMA_HDF5_INCLUDE_DIR=/home/connor/anaconda3/include/   -DARMA_DONT_USE_WRAPPER -DARMA_USE_BLAS -DARMA_USE_LAPACK -lblas -llapack
