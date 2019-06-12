gfortran -c pacotes/K_M_matrix.f90 pacotes/lapack_tools.f pacotes/lapack.f pacotes/lapack_parcer.f08 beam6gdl.f90
gfortran K_M_matrix.o lapack_tools.o lapack.o lapack_parcer.o  beam6gdl.o -o compilabeam
