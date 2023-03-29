
INC_DIR := -I/usr/local/opt/petsc/include
LIB_DIR := -L/usr/local/opt/petsc/lib

LIB := -lpetsc.3.18

CXX_FLAGS := -O0 -Wall -Werror #-g

#Compilation
all :
	gcc -o project main.c src/poisson.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)
	./project -ksp_type fgmres -pc_type lu
	make clean

#Delete of executable files
clean :
	rm project

#Delete of results
clean_txt :
	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*
