#include <cuda.h>
#include <stdlib.h>
#include <stdio.h>

extern "C" {
	#include "CUDA_header.h"
}

/* Divides 'dividend' by 'divisor', rounding up.
 */
static inline
int divisionCeil(int dividend, int divisor){
	return (dividend + divisor - 1) / divisor;
}

// Returns the first power of 2 that is >= 'base'.
static inline
int higherEqualPow2(int base){
	int result = 1;
	while(result < base) result <<= 1;
	return result;
}

/* Multi-block reduce.
 * Accepts only vectors that are power of 2.
 */
__global__ static
void reduce(int *vec, int *result){
	extern __shared__ int sdata[];

	int tid = threadIdx.x;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	sdata[tid] = vec[idx];
	__syncthreads();

	// Reduce
	for(int stride = blockDim.x >> 1; stride > 0; stride >>= 1){
		if(threadIdx.x < stride)
			sdata[threadIdx.x] += sdata[threadIdx.x+stride];

		__syncthreads();
	}

	result[blockIdx.x] = sdata[0];
}


/*
 * Collision Count procedure implemented in CUDA.
 *
 * This procedure parallelizes the sequential algorithm:
 * for i in 0:N-2
 *   for j in i+1:N-1
 *     collisions += (bead[i] == bead[j])
 * by performing just the outer 'for' in parallel.
 */
__global__
void count_collisions_cu(float3 *coords, int *result, int nCoords, int star){
	int baseIdx = blockIdx.x * 1024;
	int horizontalId = threadIdx.x + blockIdx.x * blockDim.x;

	// We read our element in a register (surplus threads will read anything)
	float3 buf = coords[horizontalId % nCoords];

	// Read first 2 blocks into shared memory
	extern __shared__ float3 sCoords[];
	sCoords[threadIdx.x] = coords[ (baseIdx + threadIdx.x) % nCoords ];
	sCoords[threadIdx.x + 1024] = coords[ (baseIdx + threadIdx.x + 1024) % nCoords ];
	__syncthreads();
	
	// Move our base index
	baseIdx = baseIdx + 2048; // We could use modulus here, but doesn't seem necessary

	// Count collisions
	int iterations = 0;
	int collisions = 0;
	int offset = 1;
	while(iterations < star){
		// Do 1024 iterations, or maybe less
		int limit = min(iterations + 1024, star);
		for(; iterations < limit; iterations++){
			// We want to check if sqrt( (Vx - Vy)(Vx - Vy) ) <= 1    Vx and Vy are float3 vectors
			// Which is the same as      (Vx - Vy)(Vx - Vy)   <= 1    and the product is an inner product
			// So we begin by taking the difference
			float3 diff = make_float3(
					buf.x - sCoords[threadIdx.x + offset].x,
					buf.y - sCoords[threadIdx.x + offset].y,
					buf.z - sCoords[threadIdx.x + offset].z
				);

			// horizontalId + iterations + 1 is the element we are comparing to
			if(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z <= 1){
				collisions += 1;
			}

			offset++;
		}
		
		// If offset == 1025, this means beads in shared memory need to be replaced
		if(offset == 1025){
			// Change blocks in shared memory when needed
			// Unfortunately we need to synchronize threads here
			__syncthreads();
			
			// Rewrite older block with earlier block
			sCoords[threadIdx.x] = sCoords[threadIdx.x + 1024];
			// Read new block
			sCoords[threadIdx.x + 1024] = coords[ (baseIdx + threadIdx.x) % nCoords ];

			// We also have to sync here
			__syncthreads();
			
			// Move base index
			baseIdx += 1024;

			offset = 1;
		}
	}

	// If the vector has an even number of elements
	// Because of this, half of the elements must execute one more iteration
	// Notice that the way the 'for...loop' above was implemented, when the
	//   code reach this point, the shared memory has valid elements for one
	//   more iteration, so we don't need to verify it again.
	// Do one more iteration:
	if(horizontalId < nCoords/2 && nCoords%2 == 0){
		float3 diff = make_float3(
				buf.x - sCoords[threadIdx.x + offset].x,
				buf.y - sCoords[threadIdx.x + offset].y,
				buf.z - sCoords[threadIdx.x + offset].z
			);

		// horizontalId + iterations + 1 is the element we are comparing to
		if(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z <= 1){
			collisions += 1;
		}
		offset++;
		iterations++;
	}

	// Sync before reducing collisions on shared memory
	__syncthreads();

	// Fill shared memory with collisions
	// We ignore collision from surplus threads
	extern __shared__ int sdata[];
	sdata[threadIdx.x] = collisions * (horizontalId < nCoords);
	__syncthreads();

	// Reduce 1024 elements
	for(int stride = 512; stride > 0; stride >>= 1){
		if(threadIdx.x < stride)
			sdata[threadIdx.x] += sdata[threadIdx.x+stride];

		__syncthreads();
	}

	// Export result
	if(threadIdx.x == 0){
		result[blockIdx.x] = sdata[0];
	}
}

/* Gets the next cuda stream in the circular list of streams.
 */
static
cudaStream_t get_next_stream(){
	const int nStreams = 8;
	static cudaStream_t streams[nStreams];
	static unsigned int launches = 0;

	// Allocate cuda streams in the first execution
	static int streamInit = 0;
	if(streamInit == 0){
		streamInit = 1;
		for(int i = 0; i < nStreams; i++){
			cudaStreamCreate(&streams[i]);
		}
	}

	launches++;
	return streams[launches%nStreams];
}

/* Given a vector with 3D coordinates of points in the space,
 *   this function calculates the number of collisions among
 *   points, using CUDA-enable GPU.
 *
 * This functions just launches the kernel, returning a
 *   structure that can later be used to fetch the result
 *   back from the device memory.
 */
extern "C" struct CollisionCountPromise
count_collisions_launch(ElfFloat3d *vector, int size){
	if(size == 0){
		CollisionCountPromise retval = { NULL, NULL };
		return retval;
	}

	float3 *d_vector;
	int *d_result;
	cudaStream_t stream = get_next_stream();

	// Allocate cuda vector for the 3D coordinates
	cudaMalloc(&d_vector, sizeof(float3) * size);
	cudaMemcpyAsync(d_vector, vector, sizeof(float3) * size, cudaMemcpyHostToDevice, stream);

	// Prepare kernel launch parameters
	const int elemInShmem = 2048; // 2048 because we need 2 blocks of 1024 elements in shmem.
	int nThreads = 1024;          // We allocate maximum number of threads per block.
	int nBlocks = divisionCeil(size, nThreads);
	int nShMem = elemInShmem * sizeof(float3); // Shared memory required

	// Calculate the number of iterations S* (S star)
	// It is the number of iterations where --all-- threads execute work
	int star;
	if(size%2 == 0)
		star = (size - 2)/2;
	else star = (size - 1)/2;

	// Allocate cuda memory for the number of collisions
	// This will also be used as a working vector for reducing among blocks
	int resultSize = higherEqualPow2(nBlocks);
	cudaMalloc(&d_result, sizeof(int) * resultSize);
	cudaMemsetAsync(d_result, 0, sizeof(int) * resultSize, stream); // Reset is needed due to size overestimation

	// Finally launch kernels
	count_collisions_cu<<<nBlocks, nThreads, nShMem, stream>>>(d_vector, d_result, size, star);

	// Reduce the result vector
	nBlocks = resultSize/1024;
	int workSize = resultSize;
	int *d_toReduce = d_result;
	int *d_reduced  = (int *) d_vector;
	while(true){
		if(nBlocks == 0){
			reduce<<<1, workSize, sizeof(int) * workSize, stream>>>(d_toReduce, d_reduced);
			break;
		}

		reduce<<<nBlocks, 1024, sizeof(int) * 1024, stream>>>(d_toReduce, d_reduced);

		// For the next run, vectors should be swapped
		int *aux = d_reduced;
		d_reduced = d_toReduce;
		d_toReduce = aux;

		// For the next run, the workSize and nBlocks are lower
		workSize = nBlocks;
		nBlocks = workSize/1024;
	}

	const struct CollisionCountPromise ret = { d_toReduce, d_reduced };
	return ret;
}

/* This procedure fetches the result from the call to the
 *   _launch correspondent.
 * The pointers within the promise structure are freed, so
 *   it shouldn't be used anywhere after a call to this function.
 */
extern "C" int count_collisions_fetch(struct CollisionCountPromise promise){
	if(promise.d_toReduce == NULL && promise.d_reduced == NULL)
		return 0;

	const int n = 1;
	int result[n];
	cudaMemcpy(&result, promise.d_reduced, sizeof(int) * n, cudaMemcpyDeviceToHost);

	cudaFree(promise.d_toReduce);
	cudaFree(promise.d_reduced);

	return result[0];
}
