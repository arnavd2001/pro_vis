#ifndef SOLUTION_PARALLEL_H
#define SOLUTION_PARALLEL_H

/** \file solution_mpi.h Routines for transmitting Solution objects to and from other nodes within an MPI environment. */

#include <mpi/mpi.h>
#include <elf_tree_comm/elf_tree_comm.h>
#include "solution.h"

#ifndef SOLUTION_PARALLEL_SOURCE_CODE
	#define SOLUTION_PARALLEL_INLINE inline
#else
	#define SOLUTION_PARALLEL_INLINE extern inline
#endif

/** Packs a Solution in the given buffer. */
SOLUTION_PARALLEL_INLINE
void Solution_pack(Solution sol, int hpSize, void *buf, int maxSize, int *position, MPI_Comm comm){
	MPI_Pack(&sol.fitness, 1, MPI_DOUBLE, buf, maxSize, position, comm);
	MPI_Pack(sol.chain, hpSize-1, MPI_CHAR, buf, maxSize, position, comm);
}

/** Unpacks a Solution and returns it. */
SOLUTION_PARALLEL_INLINE
Solution Solution_unpack(int hpSize, void *buf, int maxSize, int *position, MPI_Comm comm){
	Solution sol = Solution_blank(hpSize);
	MPI_Unpack(buf, maxSize, position, &sol.fitness, 1, MPI_DOUBLE, comm);
	MPI_Unpack(buf, maxSize, position, sol.chain, hpSize-1, MPI_CHAR, comm);
	return sol;
}

/** Calculates the fitness for all solutions in the given vector, using all nodes
 *   in the MPI communicator registered in the HIVE (HIVE_COMM.comm).
 */
SOLUTION_PARALLEL_INLINE
void Solution_calculate_fitness_master(Solution *sols, int nSols, int hpSize, MPI_Comm comm){
	int i, j;

	int commSize;
	MPI_Comm_size(comm, &commSize);

	// Allocate buffer for MPI_Scatter / Gather
	int buffSize = commSize * (hpSize - 1);
	shiftmel *buff = malloc(buffSize); // We send mov chains
	double recvBuff[commSize];        // And receive fitnesses

	for(i = 0; i < nSols; i += commSize){
		// Build scatter buffer content
		for(j = 0; j < commSize; j++){
			if((i+j) < nSols){
				memcpy(buff + j*(hpSize-1), sols[i+j].chain, hpSize - 1);
			} else {
				memset(buff + j*(hpSize-1), 0xFEFEFEFE, hpSize - 1);
			}
		}

		// Scatter buffer
		ElfTreeComm_scatter(buff, hpSize - 1, MPI_CHAR, comm);

		// Calculate own fitness
		double fit = FitnessCalc_run2(buff);
		sols[i].fitness = fit;

		// Gather fitnesses
		ElfTreeComm_gather(recvBuff, 1, MPI_DOUBLE, comm);

		// Place fitnesses into the due solutions
		for(j = 1; j < commSize && (i+j) < nSols; j++){
			sols[i+j].fitness = recvBuff[j];

			// For verifying correctness of fitness
			// int good = sols[i+j].fitness == FitnessCalc_run2(buff + j * (hpSize - 1));
		}
	}

	free(buff);
}

/** Tells slaves to return
 * Allocate buffer for MPI_Scatter / Gather */
SOLUTION_PARALLEL_INLINE
void Solution_calculate_fitness_master_kill_slaves(int hpSize, MPI_Comm comm){
	int commSize;
	MPI_Comm_size(comm, &commSize);

	int buffSize = commSize * (hpSize - 1);
	void *buff = malloc(buffSize);
	memset(buff, 0xFFFFFFFF, buffSize);
	ElfTreeComm_scatter(buff, hpSize - 1, MPI_CHAR, comm);
	free(buff);
}

/** Procedure that the slave nodes should execute.
 * Consists of waiting for migrchs, calculating its fitness, and sending the fitness back to node 0.
 * The fitness is sent back with the same MPI_TAG that was received with the migrch.
 * The slave will return once the first element of the migrch received is equal 0xFF.
 */
SOLUTION_PARALLEL_INLINE
void Solution_calculate_fitness_slave(const HPElem *chaininghp, int hpSize, MPI_Comm comm){
	int commSize;
	MPI_Comm_size(comm, &commSize);

	// Create scatter/gather buffers
	int buffSize = commSize * (hpSize - 1);
	shiftmel *buff = malloc(buffSize);
	double sendBuff[commSize];

	while(true){
		ElfTreeComm_scatter(buff, hpSize-1, MPI_CHAR, comm);
		if(0xFF == buff[0]){ // Detect end of work
			free(buff);
			return;
		}

		if(0xFE == buff[0]){ // Detect no-op
			sendBuff[0] = 0;
		} else {
			sendBuff[0] = FitnessCalc_run2(buff);
		}

		ElfTreeComm_gather(sendBuff, 1, MPI_DOUBLE, comm);
	}
}



#endif
