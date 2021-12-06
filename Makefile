# User defines
DEFS?= # e.g. -DN_HIVES=2
UFLAGS?= # e.g. -pg -g

CFL=-Wall -O2 -I src
NVCCFL=-O2 -I src
LIBS=-lm
CUDA_PRELIBS="-L/usr/local/cuda/lib64"
CUDA_LIBS=-lcuda -lcudart

MPI_LIBS=-pthread -Wl,-rpath -Wl,/usr/lib/openmpi/lib -Wl,--enable-new-dtags -L/usr/lib/openmpi/lib -lmpi
MPI_CFL=-I/usr/lib/openmpi/include/openmpi/opal/mca/event/libevent2021/libevent -I/usr/lib/openmpi/include/openmpi/opal/mca/event/libevent2021/libevent/include \
		-I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi

# We take as a rule that if any API changes, everything should be rebuilt.
# Same goes for the makefile itself
HARD_DEPS=migrch.h fitness/gyration.h fitness/CUDA_header.h fitness/fitness_private.h fitness/fitness.h \
          twirmt/twirmt.h abc_alg/hive.h abc_alg/abc_alg.h elf_tree_comm/elf_tree_comm.h numtrd.h config.h \
          solution/solution.h solution/solution_mpi.h solution/solution_structure_private.h \
          shiftmel.h random.h chaininghp.h Makefile

# This is a variable used by Makefile itself
VPATH=src/

all:
	make mpi seq

mpi:
	make milin mquard mptrd milin_threads mcuda

seq:
	make sqline squad seq_threads sqline_threads seq_cuda

milin: main.o numtrd.o measures_linear.o chaininghp.o migrch.o shiftmel.o twirmt.o acaglpal.o elf_tree_comm.o config.o hive.o gyration.o fitness.o random.o solution.o solution_mpi.o
	gcc $(CFL) $(MPI_CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS) $(MPI_LIBS)

mquard: main.o numtrd.o measures_quadratic.o chaininghp.o migrch.o shiftmel.o twirmt.o acaglpal.o elf_tree_comm.o config.o hive.o gyration.o fitness.o random.o solution.o solution_mpi.o
	gcc $(CFL) $(MPI_CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS) $(MPI_LIBS)

mptrd: main.o numtrd.o measures_threads.o chaininghp.o migrch.o shiftmel.o twirmt.o acaglpal.o elf_tree_comm.o config.o hive.o gyration.o fitness.o random.o solution.o solution_mpi.o
	gcc -fopenmp $(CFL) $(MPI_CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS) $(MPI_LIBS)

milin_threads: main.o numtrd.o measures_linear_threads.o chaininghp.o migrch.o shiftmel.o twirmt.o acaglpal.o elf_tree_comm.o config.o hive.o gyration.o fitness.o random.o solution.o solution_mpi.o
	gcc -fopenmp $(CFL) $(MPI_CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS) $(MPI_LIBS)

mcuda: main.o numtrd.o measures_cuda.o CUDA_collision_count.o CUDA_contact_count.o chaininghp.o migrch.o shiftmel.o twirmt.o acaglpal.o elf_tree_comm.o config.o hive.o gyration.o fitness.o random.o solution.o solution_mpi.o
	gcc $(CUDA_PRELIBS) $(CFL) $(MPI_CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS) $(MPI_LIBS) $(CUDA_LIBS)

sqline: main.o numtrd.o measures_linear.o chaininghp.o migrch.o shiftmel.o twirmt.o acalg_seq.o config.o hive.o gyration.o fitness.o random.o solution.o
	gcc $(CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS)

squad: main.o numtrd.o measures_quadratic.o chaininghp.o migrch.o shiftmel.o twirmt.o acalg_seq.o config.o hive.o gyration.o fitness.o random.o solution.o
	gcc $(CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS)

seq_threads: main.o numtrd.o measures_threads.o chaininghp.o migrch.o shiftmel.o twirmt.o acalg_seq.o config.o hive.o gyration.o fitness.o random.o solution.o
	gcc -fopenmp $(CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS)

sqline_threads: main.o numtrd.o measures_linear_threads.o chaininghp.o migrch.o shiftmel.o twirmt.o acalg_seq.o config.o hive.o gyration.o fitness.o random.o solution.o
	gcc -fopenmp $(CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS)

seq_cuda: main.o numtrd.o measures_cuda.o CUDA_collision_count.o CUDA_contact_count.o chaininghp.o migrch.o shiftmel.o twirmt.o acalg_seq.o config.o hive.o gyration.o fitness.o random.o solution.o
	gcc $(CUDA_PRELIBS) $(CFL) $(UFLAGS) $(DEFS) $^ -o $@ $(LIBS) $(CUDA_LIBS)

clean:
	find -name "*~" -type f -exec rm -vf '{}' \;
	find -name "*.o" -type f -exec rm -vf '{}' \;
	rm -vf *~ gmon.out

clean_all: clean
	rm -vf milin mquard mptrd milin_threads mcuda sqline squad seq_threads sqline_threads seq_cuda

dox:
	doxygen Doxyfile


# Explicit non-MPI object rules (WHEN ADDING NEW, MUST ADD TO $(BIN) TARGET TOO)
main.o:               main.c $(HARD_DEPS)
numtrd.o:              numtrd.c $(HARD_DEPS)
measures_quadratic.o: fitness/measures_quadratic.c $(HARD_DEPS)
measures_linear.o:    fitness/measures_linear.c $(HARD_DEPS)
chaininghp.o:            chaininghp.c $(HARD_DEPS)
migrch.o:           migrch.c $(HARD_DEPS)
shiftmel.o:            shiftmel.c $(HARD_DEPS)
twirmt.o:             twirmt/twirmt.c $(HARD_DEPS)
acalg_seq.o: abc_alg/acalg_seq.c $(HARD_DEPS)
config.o:             config.c $(HARD_DEPS)
hive.o:               abc_alg/hive.c $(HARD_DEPS)
gyration.o:           fitness/gyration.c $(HARD_DEPS)
fitness.o:            fitness/fitness.c $(HARD_DEPS)
random.o:             random.c $(HARD_DEPS)
solution.o:           solution/solution.c $(HARD_DEPS)


# Explicit CUDA object rules
CUDA_collision_count.o: fitness/CUDA_collision_count.cu $(HARD_DEPS)
	nvcc $(NVCCFL) -c $(DEFS) -o "$@" "$<" $(LIBS)

CUDA_contact_count.o: fitness/CUDA_contact_count.cu $(HARD_DEPS)
	nvcc $(NVCCFL) -c $(DEFS) -o "$@" "$<" $(LIBS)

measures_cuda.o: fitness/measures_cuda.c $(HARD_DEPS)
	gcc $(CUDA_PRELIBS) $(CFL) -c $(DEFS) -o "$@" "$<" $(LIBS) $(CUDA_LIBS)


# Explicit OpenMP object rules
measures_threads.o: fitness/measures_threads.c $(HARD_DEPS)
	gcc -c $(DEFS) $(CFL) -fopenmp $(UFLAGS) -o "$@" "$<" $(LIBS)

measures_linear_threads.o: fitness/measures_linear_threads.c $(HARD_DEPS)
	gcc -c $(DEFS) $(CFL) -fopenmp $(UFLAGS) -o "$@" "$<" $(LIBS)

# Explicit MPI object rules
acaglpal.o: abc_alg/acaglpal.c $(HARD_DEPS)
	gcc -c $(DEFS) $(CFL) $(MPI_CFL) $(UFLAGS) -o "$@" "$<" $(LIBS) $(MPI_LIBS)

elf_tree_comm.o: elf_tree_comm/elf_tree_comm.c $(HARD_DEPS)
	gcc -c $(DEFS) $(CFL) $(MPI_CFL) $(UFLAGS) -o "$@" "$<" $(LIBS) $(MPI_LIBS)

solution_mpi.o: solution/solution_mpi.c $(HARD_DEPS)
	gcc -c $(DEFS) $(CFL) $(MPI_CFL) $(UFLAGS) -o "$@" "$<" $(LIBS) $(MPI_LIBS)

# Implicit rule for building objects
%.o:
	gcc -c $(DEFS) $(CFL) $(UFLAGS) -o "$@" "$<" $(LIBS)
