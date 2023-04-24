
INC_DIR := -I/usr/local/opt/petsc/include
LIB_DIR := -L/usr/local/opt/petsc/lib

LIB := -lpetsc.3.18

FLAGS := -Wall -Werror -O3
DEBUG_FLAGS := -Wall -Werror -g -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls

SOURCES = main.c src/poisson.c src/functions.c

exe : 
	clang -o project $(SOURCES) -lm -lmpi $(FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)

#Compilation
all :
	clang -o project $(SOURCES) -lm -lmpi $(FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)
	./project -ksp_type fgmres -pc_type lu
	make clean

debug :
	clang -o project $(SOURCES) -lm -lmpi $(DEBUG_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)
	./project -ksp_type fgmres -pc_type lu
	make clean

#Delete of executable files
clean :
	rm project

#Delete of results
clean_files :
	rm -vf output/*.txt output/P-* output/U-* output/V-* output/Reh-* output/Vtx-* output/Rehw-* output/Div-* output/w-*
