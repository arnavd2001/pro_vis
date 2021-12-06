
#include <mpi/mpi.h>
#include <stdio.h>

void ElfTreeComm_scatter(void *buf, int sendCount, MPI_Datatype type, MPI_Comm comm){
	int myRank, commSize;
	MPI_Comm_rank(comm, &myRank);
	MPI_Comm_size(comm, &commSize);

	if(commSize == 1) return;

	int typeSize;
	MPI_Type_size(type, &typeSize);

	// Find lowest power of 2 higher than or equal to commSize
	int hipow2 = 1;
	while(hipow2 < commSize) hipow2 <<= 1;

	int control;
	for(control = hipow2 / 2; control >= 1; control /= 2){
		if( myRank % control != 0 ) continue;

		int myDiv = myRank / control;
		if(myDiv % 2 == 0){
			// Send
			int dest = myRank + control;
			if(dest >= commSize) continue;

			int farDest = dest + control;
			if(farDest > commSize) farDest = commSize;

			int dataCount = (farDest - dest) * sendCount;
			MPI_Send( ((char*) buf) + control * sendCount * typeSize, dataCount, type, dest, 0, comm);
		} else {
			// Receive
			int src = myRank - control;
			
			int farDest = myRank + control;
			if(farDest > commSize) farDest = commSize;

			int dataCount = (farDest - myRank) * sendCount;
			MPI_Recv(buf, dataCount, type, src, 0, comm, MPI_STATUS_IGNORE);
		}
	}
}

void ElfTreeComm_gather(void *buf, int sendCount, MPI_Datatype type, MPI_Comm comm){
	int myRank, commSize;
	MPI_Comm_rank(comm, &myRank);
	MPI_Comm_size(comm, &commSize);

	if(commSize == 1) return;

	int typeSize;
	MPI_Type_size(type, &typeSize);

	// Find lowest power of 2 higher than or equal to commSize
	int hipow2 = 1;
	while(hipow2 < commSize) hipow2 <<= 1;

	int control;
	for(control = 1; control < hipow2; control *= 2){
		if( myRank % control != 0 ) continue;

		int myDiv = myRank / control;
		if(myDiv % 2 == 0){
			// Receive
			int src = myRank + control;
			if(src >= commSize) continue;

			int farDest = src + control;
			if(farDest > commSize) farDest = commSize;

			int dataCount = (farDest - src) * sendCount;
			MPI_Recv( ((char*) buf) + control * sendCount * typeSize, dataCount, type, src, 0, comm, MPI_STATUS_IGNORE);
		} else {
			// Send
			int dest = myRank - control;
			
			int farDest = myRank + control;
			if(farDest > commSize) farDest = commSize;

			int dataCount = (farDest - myRank) * sendCount;
			MPI_Send(buf, dataCount, type, dest, 0, comm);
		}
	}
}


/* DEBUG PROCEDURES

#include <unistd.h>
#include <string.h>

#define STR_SIZE 15

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);

	int myRank, commSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);

	// We will be scattering strings of size STR_SIZE
	// And then gathering them back

	// Create scatter/gather buffer
	int bufSize = commSize * STR_SIZE + 1;
	char buf[bufSize];
	char recv[bufSize];
	buf[bufSize - 1] = '\0';
	recv[bufSize - 1] = '\0';


	int i, j;
	if(myRank == 0){
		for(i = 0; i < commSize; i++){
			for(j = 0; j < STR_SIZE; j++){
				buf[i * STR_SIZE + j] = 'a' + i;
			}
		}

		// We also have to fill recv for the section that belongs to node 0
		// Because it needs to be there when we give 'recv' to the gather procedure
		for(i = 0; i < STR_SIZE; i++)
			recv[i] = 'a';

		printf("%s\n", buf);
	}

	// Scatter and then each process prints their buffer content
	ElfTreeComm_scatter(buf, STR_SIZE, MPI_CHAR, MPI_COMM_WORLD);
	
	// Temporarily add \0 to the buffer
	char aux = buf[STR_SIZE];
	buf[STR_SIZE] = '\0';

	// Print buf
	usleep(1000 * myRank);
	printf("%d: %s\n", myRank, buf);

	// Restore buffer state (only matters for node 0)
	buf[STR_SIZE] = aux;

	MPI_Barrier(MPI_COMM_WORLD);

	// Gather and then process 0 prints the buffer content
	if(myRank != 0){
		ElfTreeComm_gather(buf, STR_SIZE, MPI_CHAR, MPI_COMM_WORLD);
	} else {
		ElfTreeComm_gather(recv, STR_SIZE, MPI_CHAR, MPI_COMM_WORLD);
	}

	if(myRank == 0){
		printf("\n");
		printf("Sizes: %lu %lu\n", strlen(buf), strlen(recv));
		
		for(i = 0; i < bufSize; i++){
			if(buf[i] != recv[i]){
				printf("Difference at: %d\n", i);
			}
		}


		printf("%s\n", recv);
	}

	MPI_Finalize();
	return 0;
}


*/
