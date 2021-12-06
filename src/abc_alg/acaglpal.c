#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi/mpi.h>

#include <migrch.h>
#include <chaininghp.h>
#include <fitness/fitness.h>
#include <random.h>
#include <solution/solution_mpi.h>

#include "abc_alg.h"
#include "hive.h"

struct {
	MPI_Comm comm;
	int      size;
} HIVE_COMM;

/* Performs the forager phase of the searching cycle
 * Procedure idea:
 *   For each solution, generate a new one in the neighborhood
 *   replace the varied solution if it was improved
 */
static
void parallel_forager_phase(int hpSize){
	int i;
	Solution sols[HIVE_nSols()];

	// Generate new random solutions
	for(i = 0; i < HIVE_nSols(); i++)
		sols[i] = HIVE_perturb_solution(i, hpSize);

	// Calculate fitnesses
	Solution_calculate_fitness_master(sols, HIVE_nSols(), hpSize, HIVE_COMM.comm);

	// Replace solutions in the HIVE
	for(i = 0; i < HIVE_nSols(); i++)
		HIVE_try_replace_solution(sols[i], i, hpSize);
}

/* Performs the onlooker phase of the searching cycle
 * Procedure idea;
 *   Calculate the SUM of fitnesses for all solutions
 *   Fitness can be negative, so we add a BASE that is the lowest fitness found
 *   For each solution SOL, (SOL.fitness/SUM) is its probability PROB of being perturbed
 *   (PROB * nOnlookers) is the number of perturbations that should be generated
 */
static
void parallel_onlooker_phase(int hpSize){
	int i, j;
	int nOnlookers = COLONY_SIZE - (COLONY_SIZE * FORAGER_RATIO);

	Solution sols[nOnlookers + HIVE_nSols()]; // Overestimate due to possible rounding errors.
	int indexes[nOnlookers + HIVE_nSols()];   // Stores indexes where each solution belong
	int nSols;

	// Find the minimum (If no negative numbers, min should be 0)
	double min = 0;
	for(i = 0; i < HIVE_nSols(); i++){
		double fit = Solution_fitness(HIVE_solution(i));
		if(fit < min)
			min = fit;
	}

	// Sum the 'normalized' fitnesses
	double sum = 0;
	for(i = 0; i < HIVE_nSols(); i++){
		sum += Solution_fitness(HIVE_solution(i)) - min;
	}

	// For each solution, count the number of onlooker bees that should perturb it
	//   then add perturbed solutions into the sols vector
	nSols = 0;
	for(i = 0; i < HIVE_nSols(); i++){
		double norm = Solution_fitness(HIVE_solution(i)) - min;
		double prob = norm / sum; // The probability of perturbing such solution

		// Count number of onlookers that should perturb such solution
		int nIter = round(prob * nOnlookers);

		// Generate perturbations
		for(j = 0; j < nIter; j++){
			sols[nSols] = HIVE_perturb_solution(i, hpSize);
			indexes[nSols] = i;
			nSols++;
		}
	}

	// Calculate fitness
	Solution_calculate_fitness_master(sols, nSols, hpSize, HIVE_COMM.comm);

	// Replace solutions where due
	for(i = 0; i < nSols; i++)
		HIVE_try_replace_solution(sols[i], indexes[i], hpSize);
}

/* Performs the scout phase of the searching cycle
 * Procedure idea:
 *   Find all the solutions whose idle_iterations exceeded the limit
 *   Generate replacement solutions
 *   Calculate fitness
 *   Replace solutions
 */
static
void parallel_scout_phase(int hpSize){
	int i;

	Solution sols[HIVE_nSols()];
	int indexes[HIVE_nSols()];
	int nSols = 0;

	// Find idle solutions
	for(i = 0; i < HIVE_nSols(); i++){
		int idle = Solution_idle_iterations(HIVE_solution(i));
		if(idle > IDLE_LIMIT)
			indexes[nSols++] = i;
	}

	// Generate random solutions
	for(i = 0; i < nSols; i++)
		sols[i] = Solution_random(hpSize);

	// Calculate fitness
	Solution_calculate_fitness_master(sols, nSols, hpSize, HIVE_COMM.comm);

	// Replace solutions
	for(i = 0; i < nSols; i++)
		HIVE_force_replace_solution(sols[i], indexes[i]);
}

/* Exchanges solutions among the hives.
 * 'ringComm' should be the communicator containing the masters of each hive.
 */
static
void ring_exchange(MPI_Comm ringComm, int hpSize){
	int commSize, myRank;
	MPI_Comm_size(ringComm, &commSize);
	MPI_Comm_rank(ringComm, &myRank);

	// If there is only 1 process, it's not a ring.
	if(commSize == 1) return;

	// Get solutions to send
	Solution randSol = HIVE_solutions()[urandom_max(HIVE_nSols())];
	Solution bestSol = HIVE_best_sol();

	// Create input/output buffers
	int maxSize = 2 * hpSize + sizeof(double) * 2 + 32; // We overestimate a bit
	char inBuf[maxSize];
	char outBuf[maxSize];

	// Pack data
	int position = 0;
	Solution_pack(bestSol, hpSize, outBuf, maxSize, &position, ringComm);
	Solution_pack(randSol, hpSize, outBuf, maxSize, &position, ringComm);

	// Send data
	int src, dest;
	if(myRank % 2 == 0){
		dest = (myRank+1) % commSize;
		MPI_Send(outBuf, position, MPI_PACKED, dest, 0, ringComm);

		src = myRank == 0 ? commSize-1 : myRank-1;
		MPI_Recv(inBuf, position, MPI_PACKED, src, 0, ringComm, MPI_STATUS_IGNORE);
	} else {
		src = myRank-1;
		MPI_Recv(inBuf, position, MPI_PACKED, src, 0, ringComm, MPI_STATUS_IGNORE);

		dest = (myRank+1) % commSize;
		MPI_Send(outBuf, position, MPI_PACKED, dest, 0, ringComm);
	}

	// Unpack data
	position = 0;
	Solution sol1 = Solution_unpack(hpSize, inBuf, maxSize, &position, ringComm);
	Solution sol2 = Solution_unpack(hpSize, inBuf, maxSize, &position, ringComm);

	int ridx1 = urandom_max(HIVE_nSols());
	HIVE_force_replace_solution(sol1, ridx1);

	int ridx2 = urandom_max(HIVE_nSols());
	HIVE_force_replace_solution(sol2, ridx2);
}

/* Gathers the best solutions among the hives in node 0.
 * 'ringComm' should be the communicator containing the masters of each hive.
 * The HIVE in node 0 is altered so that HIVE.best is the best among all best solutions
 *   of all hives.
 */
static
void ring_gather(MPI_Comm ringComm, int hpSize){
	int i, commSize, myRank;
	MPI_Comm_size(ringComm, &commSize);
	MPI_Comm_rank(ringComm, &myRank);

	// If there is only one process, there is nothing to be done.
	if(commSize == 1) return;

	// Get my solution
	Solution sol = HIVE_best_sol();

	// Create gather buffer
	int maxSize = commSize * (hpSize + sizeof(double) + 32); // We overestimate a bit
	char *gatBuf = malloc(maxSize);

	// Pack my solution
	int position = 0;
	Solution_pack(sol, hpSize, gatBuf, maxSize, &position, ringComm);

	// Gather solutions
	int byteCount = position;
	ElfTreeComm_gather(gatBuf, byteCount, MPI_PACKED, ringComm);

	// Find best solution
	if(myRank == 0){
		position = 0;
		for(i = 0; i < commSize; i++){
			sol = Solution_unpack(hpSize, gatBuf, maxSize, &position, ringComm);

			if(Solution_fitness(sol) > Solution_fitness(HIVE_best_sol())){
				HIVE_replace_best(sol);
			} else {
				Solution_free(sol);
			}
		}
	}

	free(gatBuf);
}

Solution ABC_predict_structure(const HPElem * chaininghp, int hpSize, int nCycles, PredResults *results){
	MPI_Init(NULL, NULL);
	int commSize, myRank;
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	/* We will divide COMM_WORLD into N_HIVES groups with same number of nodes each.
	 * If COMM_WORLD has the nodes:
	 * 0 1 2 3 4 5
	 * And N_HIVES == 2, then we will only place adjacent nodes in each subgroup:
	 * (0 1 2) and (3 4 5)    (in contrast with, for example [0 2 4] and [1 3 5])
	 *
	 * So upon program call, keep this in mind in order keep processes of a single HIVE
	 *   nearby, such as in a single multi-core processor.
	 */
	int nodesPerHive = commSize / N_HIVES;

	if(N_HIVES > commSize){
		if(myRank == 0)
			fprintf(stderr, "Number of Hives cannot be greater than number of launched processes!\n");
		MPI_Finalize();
		exit(0);
	}

	if(commSize % N_HIVES != 0){
		if(myRank == 0)
			fprintf(stderr, "Number of Hives must divide the number of nodes!\n");
		MPI_Finalize();
		exit(0);
	}

	MPI_Comm hiveComm;
	int myColor = myRank / nodesPerHive;
	MPI_Comm_split(MPI_COMM_WORLD, myColor, myRank, &hiveComm);

	HIVE_initialize();
	HIVE_COMM.comm = hiveComm;
	HIVE_COMM.size = nodesPerHive;
	FitnessCalc_initialize(chaininghp, hpSize);

	int myHiveRank, myWorldRank;
	MPI_Comm_rank(hiveComm, &myHiveRank);
	MPI_Comm_rank(MPI_COMM_WORLD, &myWorldRank);

	// We build the communicator for the ring topology
	int color = myHiveRank == 0 ? 0 : MPI_UNDEFINED;
	MPI_Comm ringComm;
	MPI_Comm_split(MPI_COMM_WORLD, color, myWorldRank, &ringComm);

	Solution retval;
	if(myHiveRank != 0){
		Solution_calculate_fitness_slave(chaininghp, hpSize, HIVE_COMM.comm);
		results->fitness = -1;
		results->contactsH = -1;
		results->collisions = -1;
		results->bbGyration = -1;
	} else {
		int i;

		for(i = 0; i < nCycles; i++){

			parallel_forager_phase(hpSize);

			parallel_onlooker_phase(hpSize);

			/* For each solution, check its idle_iterations
			 * If it exceeded the limit, replace it with a new random solution
			 */
			parallel_scout_phase(hpSize);

			const int migCycle = nCycles * 0.1; // Migration cycle
			if( i != 0 && (i % migCycle == 0) ){
				ring_exchange(ringComm, hpSize);
			}

			HIVE_increment_cycle();
		}

		ring_gather(ringComm, hpSize);

		retval = HIVE_best_sol();
		double fit = Solution_fitness(retval);

		if(results && myWorldRank == 0){
			results->fitness = fit;
			FitnessCalc_measures(Solution_chain(retval), &results->contactsH, &results->collisions, &results->bbGyration);
		}

		// Tell slaves to stop
		Solution_calculate_fitness_master_kill_slaves(hpSize, HIVE_COMM.comm);
	}

	MPI_Barrier(hiveComm);
	FitnessCalc_cleanup();
	HIVE_destroy();
	MPI_Comm_free(&hiveComm);
	MPI_Finalize();

	return retval;
}
