gfortran -c  K_M_matrix.f90 lapack_tools.f lapack.f lapack_parcer.f08 beam6gdl.f90
gfortran  K_M_matrix.o lapack_tools.o lapack.o lapack_parcer.o  beam6gdl.o -o compilabeam
