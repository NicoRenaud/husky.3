# ==== ====================================== ====
# ==== 	        HUSKY.3  Makefile 	      ====
# ==== ====================================== ====
# Automatically generated the  10-01-2017 at 17h25


CC   = gcc

OBJ = OBJ/algebra_cplx.o OBJ/algebra.o OBJ/extract_hamiltonian.o OBJ/spec.o OBJ/huckel.o OBJ/read_input.o OBJ/main.o OBJ/export_file.o OBJ/export_gnuplot.o OBJ/print_screen.o OBJ/computeDOS.o OBJ/computeMO.o OBJ/moCube.o OBJ/current.o OBJ/negf_te.o 

LIBS = -lm -lgslcblas -llapack -lgsl
BIN  =  husky
CFLAGS = -Wall
RM = rm -f

.PHONY: all all-before all-after clean clean-custom
all: all-before $(BIN) all-after

clean: clean-custom
		${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $(BIN) $(LIBS)

OBJ/algebra_cplx.o : ./SRC/algebra_cplx/algebra_cplx.c
	$(CC) -c ./SRC/algebra_cplx/algebra_cplx.c -o OBJ/algebra_cplx.o $(CFLAGS)

OBJ/algebra.o : ./SRC/algebra_real/algebra.c
	$(CC) -c ./SRC/algebra_real/algebra.c -o OBJ/algebra.o $(CFLAGS)

OBJ/extract_hamiltonian.o : ./SRC/algebra_real/extract_hamiltonian.c
	$(CC) -c ./SRC/algebra_real/extract_hamiltonian.c -o OBJ/extract_hamiltonian.o $(CFLAGS)

OBJ/spec.o : ./SRC/algebra_real/spec.c
	$(CC) -c ./SRC/algebra_real/spec.c -o OBJ/spec.o $(CFLAGS)

OBJ/huckel.o : ./SRC/huckel/huckel.c
	$(CC) -c ./SRC/huckel/huckel.c -o OBJ/huckel.o $(CFLAGS)

OBJ/read_input.o : ./SRC/huckel/read_input.c
	$(CC) -c ./SRC/huckel/read_input.c -o OBJ/read_input.o $(CFLAGS)

OBJ/main.o : ./SRC/main.c
	$(CC) -c ./SRC/main.c -o OBJ/main.o $(CFLAGS)

OBJ/export_file.o : ./SRC/manag_output/export_file.c
	$(CC) -c ./SRC/manag_output/export_file.c -o OBJ/export_file.o $(CFLAGS)

OBJ/export_gnuplot.o : ./SRC/manag_output/export_gnuplot.c
	$(CC) -c ./SRC/manag_output/export_gnuplot.c -o OBJ/export_gnuplot.o $(CFLAGS)

OBJ/print_screen.o : ./SRC/manag_output/print_screen.c
	$(CC) -c ./SRC/manag_output/print_screen.c -o OBJ/print_screen.o $(CFLAGS)

OBJ/computeDOS.o : ./SRC/manipMol/computeDOS.c
	$(CC) -c ./SRC/manipMol/computeDOS.c -o OBJ/computeDOS.o $(CFLAGS)

OBJ/computeMO.o : ./SRC/manipMol/computeMO.c
	$(CC) -c ./SRC/manipMol/computeMO.c -o OBJ/computeMO.o $(CFLAGS)

OBJ/moCube.o : ./SRC/manipMol/moCube.c
	$(CC) -c ./SRC/manipMol/moCube.c -o OBJ/moCube.o $(CFLAGS)

OBJ/current.o : ./SRC/negf_te/current.c
	$(CC) -c ./SRC/negf_te/current.c -o OBJ/current.o $(CFLAGS)

OBJ/negf_te.o : ./SRC/negf_te/negf_te.c
	$(CC) -c ./SRC/negf_te/negf_te.c -o OBJ/negf_te.o $(CFLAGS)

