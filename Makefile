# Directories
S = src
L = lib

# ZARCH := $(shell ./get_arch.sh)
ZARCH = native
ZTUNE = native

# Compiler flags.
CFLAGS_S = -Wall -O3 -march=native -mtune=native -lm

#
LPCK_FLAGS = -llapacke -llapack -lblas 
CFLAGS_LPK = -Wall -O3 -march=native -mtune=native ${LPCK_FLAGS} -lm

#
CFLAGS_Z = -O3 -march=${ZARCH} -mtune=${ZTUNE} -Wall -pthread 
PFLAGS_Z = ${CFLAGS_Z} -pthread
#
CFLAGS_N = -Wall -O3 -march=native -mtune=native -lm
PFLAGS_N = ${CFLAGS_N} -pthread 

# Valgrind flags.
VFLAGS_LPK = -Wall -O0 -g -march=native -mtune=native ${LPCK_FLAGS} -lm
VFLAGS = -Wall -O0 -g -march=native -mtune=native -lm

# Default flags.
CFLAGS = ${CFLAGS_S}

ALL: mc_search

mc_search:               \
${L}/main.o              \
${L}/help.o              \
${L}/matrix_operations.o \
${L}/parameters_read.o   \
${L}/data_read.o         \
${L}/random_walk.o       \
${L}/positions_print.o   \
${L}/positions_calc.o
	gcc -o $@ $^ -lm


${L}/main.o:             \
main.c                   \
pcg_random.h             \
prm_def.h                \
${L}/parameters_read.o   \
${L}/data_read.o 		 \
${L}/positions_calc.o    \
${L}/positions_print.o   \
${L}/random_walk.o       \
${L}/help.o
	gcc -o $@ $< ${CFLAGS} -c


${L}/parameters_read.o:  \
parameters_read.c        \
prm_def.h                \
${L}/help.o
	gcc -o $@ $< ${CFLAGS} -c

${L}/data_read.o:        \
data_read.c              \
prm_def.h
	gcc -o $@ $< ${CFLAGS}  -c

${L}/matrix_operations.o: \
matrix_operations.c       \
prm_def.h
	gcc -o $@ $< ${CFLAGS} -c

${L}/positions_calc.o:   \
positions_calc.c         \
prm_def.h
	gcc -o $@ $< ${CFLAGS} -c

${L}/positions_print.o:  \
positions_print.c
	gcc -o $@ $< ${CFLAGS} -c

${L}/random_walk.o:      \
random_walk.c
	gcc -o $@ $< ${CFLAGS} -c

${L}/help.o:             \
help.c
	gcc -o $@ $< ${CFLAGS} -c

clean:
	\rm -rf *~ *~ ${L}/*.o

veryclean: clean
	\rm -rf *mc_search*

strip:
	for f in ${ALL} ; do strip -s $$f ; done
