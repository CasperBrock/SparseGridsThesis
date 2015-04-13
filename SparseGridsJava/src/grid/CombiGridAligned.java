package grid;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.RecursiveAction;

import gridFunctions.GridFunction;
import gridFunctions.GridFunctions;

public class CombiGridAligned {

	public int[] levels;
	public int dimensions;
	public double[] grid;
	public int gridSize;
	public int[] pointsPerDimension;
	public int alignment;
	int noAligned;
	int arraySize;
	int[] strides; //used for the recursive code.
	public int recTile=5; //Some values are hard coded in the original.
	public double recTallPar=0.3;
	float [] rf; 
	public int recMaxSpawn=14; //Public, for varying within the test.
	public int recMinSpawn=13; //Public, for varying within the test.



	public static void main(String[] args) {		
		int[] levels = {3, 3, 3};
		CombiGridAligned grid = new CombiGridAligned(levels, 32);
		CombiGridAligned grid2 = new CombiGridAligned(levels, 32);
		//System.out.println("Gridsize: " + grid.gridSize);
		//System.out.println("Arraysize: " + grid.arraySize);
		/*for(int i = 0; i < 10000; i++) {
			//grid.setValues(GridFunctions.ALLONES);
			//grid.hierarchizeOptimized(4);
		}*/
		grid2.setValues(GridFunctions.ALLONES);
		//grid2.hierarchizeOptimized(4);
		//grid2.hierarchizeRecursiveThreads(2);
		grid2.hierarchizeOptimizedThreads(32);
		grid.setValues(GridFunctions.ALLONES);
		grid.hierarchizeRecursive();
		if (grid.compare(grid2)){
			System.out.println("The grids are equal.");
		} else System.out.println("not equal grids. check code.");
		//grid2.printValues();
		 
		grid.printValues();
		//grid2.setValues(GridFunctions.ALLONES);
		//grid2.hierarchizeOptimizedParallelStream(16, 1000);
		//grid2.printValues();
		//System.out.println(grid2.compare(grid));
	}

	public CombiGridAligned(int[] levels, int alignment) {
		// alignment in bytes
		// alignment multiple of 32 bytes for AVX
		// grid needs to be aligned for the blocked (optimized) version of the
		// code= new float[dimensions];
		dimensions = levels.length;
		this.levels = new int[dimensions];
		pointsPerDimension = new int[dimensions];
		int i;
		// lengthen 1st dimensions
		pointsPerDimension[0] = myPow2(levels[0]) - 1;
		this.levels[0] = levels[0];
		noAligned = (int) (Math.ceil((double) (myPow2(levels[0]) - 1)
				/ alignment * 8.0)
				* alignment / 8.0);
		gridSize = pointsPerDimension[0];
		arraySize = noAligned;
		strides = new int[dimensions+1];
		strides[0]=1;
		for (i = 1; i < dimensions; i++) {
			this.levels[i] = levels[i];
			pointsPerDimension[i] = myPow2(levels[i]) - 1;
			gridSize *= pointsPerDimension[i];
			strides[i] = arraySize;
			arraySize *= pointsPerDimension[i];
			//TODO fill Strides[]. See line 454 and forward in c++ code.
		}
		strides[dimensions]=arraySize;
		grid = new double[arraySize];


		rf =  new float[dimensions];
		for(int j=0;j<dimensions;j++) rf[j]=(float) 1.0;


	}

	private void printValues2DArr(int size, int offset, int n0) {
		for (int i = 1; i <= size; i++) {
			System.out.print(grid[offset + i - 1]);
			System.out.print('\t');
			if (i % n0 == 0)
				System.out.print('\n');
		}
		System.out.print('\n');
	}

	public void printValues() {
		System.out.println();
		if (dimensions <= 2)
			printValues2DArr(arraySize, 0, noAligned);
		else {
			int[] currentLevels = new int[dimensions];
			int chunkSize = noAligned * pointsPerDimension[1];
			int numberOfChunks = arraySize / chunkSize;
			for (int ctr = 0; ctr < numberOfChunks; ctr++) {
				printValues2DArr(chunkSize, ctr * chunkSize, noAligned);
				currentLevels[2]++;
				for (int dd = 2; dd < dimensions; dd++) {
					if (currentLevels[dd] == pointsPerDimension[dd]) {
						currentLevels[dd] = 0;
						if (dd + 1 < dimensions)
							currentLevels[dd + 1]++;
					}
				}
			}
		}
	}

	public boolean compare(CombiGridAligned cga) {
		if (cga.arraySize != arraySize)
			return false;
		for(int i = 0; i < arraySize; i++)
			if(Math.abs(cga.grid[i] - grid[i]) > 0.000001) {
				System.out.println(cga.grid[i]);
				System.out.println(grid[i]);
				return false;
			}
		return true;
	}

	public void setValues(GridFunctions func) {
		alignment = alignment / 8;
		double[] stepsize = new double[dimensions];
		double[] x = new double[dimensions];
		int[] levelSets = new int[dimensions];
		levelSets[0] = 1;
		int[] dimensionCounter = new int[dimensions];

		for (int d = 1; d < dimensions; d++) {
			levelSets[d] = levelSets[d - 1] * pointsPerDimension[d - 1];
		}

		for (int d = 0; d < dimensions; d++) {
			stepsize[d] = Math.pow(2, -levels[d]);
			dimensionCounter[d] = 1;
		}

		int pos = 0;
		for (int counter = 0; counter < gridSize; counter++) {
			for (int d = 0; d < dimensions; d++) {
				if (dimensionCounter[d] > pointsPerDimension[d]
						&& d == dimensions - 1) {
					System.out.println("Something went wrong");
					return;
				}

				if (d == 0 && dimensionCounter[d] > pointsPerDimension[d]) { // we
					// are
					// in
					// first
					// dim
					// (padded!)
					// and
					// need
					// to
					// increment
					// pos
					grid[pos] = Double.POSITIVE_INFINITY; // pos points to
					// padded point
					pos++;
					dimensionCounter[d] = 1;
					dimensionCounter[d + 1]++;
				}

				if (dimensionCounter[d] > pointsPerDimension[d]
						&& d < dimensions - 1) {
					dimensionCounter[d] = 1;
					dimensionCounter[d + 1]++;
				}
			}

			for (int d = 0; d < dimensions; d++) {
				x[d] = dimensionCounter[d] * stepsize[d];
			}

			grid[pos] = GridFunction.call(x, func);
			dimensionCounter[0]++;
			pos++;
		}

		grid[arraySize - 1] = Double.POSITIVE_INFINITY;
		return;
	}

	public void hierarchize1DUnoptimized(int start, int stride, int size, int dimension) {
		int level, steps, ctr, offset, parentOffset, stepsize, parOffsetStrided;
		double val1, val2, val3, left, right;

		level = levels[dimension];
		steps = myPow2(level - 1);
		offset = 0;
		stepsize = 2;
		parentOffset = 1;
		left = 0;
		right = 0;

		for(level--; level > 1; level--) {
			parOffsetStrided = parentOffset*stride;
			grid[start + offset * stride] -= 0.5 * grid[start + offset * stride + parOffsetStrided];
			offset += stepsize;
			left = 0.5 * grid[start + offset * stride - parOffsetStrided];
			for (ctr = 1; ctr < steps - 1; ctr++) {				
				val1 = grid[start + offset * stride];
				right = 0.5 * grid[start + offset * stride + parOffsetStrided];
				val2 = val1 - left;
				val3 = val2 - right;
				grid[start+offset*stride] = (int) val3;
				left = right;
				offset += stepsize;
			} 

			grid[start + offset * stride] -= right;
			//steps = steps >> 1;
			steps = steps / 2;
			offset = myPow2(levels[dimension] - level) - 1;
			parentOffset =  stepsize;
			//stepSize = stepSize << 1;
			stepsize = stepsize * 2;
		}

		right = 0.5 * grid[start + (offset + parentOffset) * stride];
		grid[start + offset * stride] -= right;
		offset += stepsize;
		grid[start + offset * stride] -= right;
	}

	public void hierarchize1DOptimized(int start, int stride, int size,
			int dimension, int unroll) {
		// optimized with vector operations. lvl for any dim must be larger
		// than 1.
		int level;
		int steps;
		int offset, parentOffset;
		int stepsize;
		int parOffsetStrided;
		double val1, val2, val3, val4, right1, right2, right3, right4, left1, left2, left3, left4;
		level = levels[dimension];
		steps = myPow2(level - 1);
		offset = 0;
		stepsize = 2;
		parentOffset = 1;

		for (level--; level > 1; level--) {
			parOffsetStrided = parentOffset * stride;

			for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) {
				right1 = grid[start + offset * stride + poleLoop] * 0.5;
				right2 = grid[start + offset * stride + poleLoop + 1] * 0.5;
				right3 = grid[start + offset * stride + poleLoop + 2] * 0.5;
				right4 = grid[start + offset * stride + poleLoop + 3] * 0.5;

				val1 = grid[start + offset * stride + poleLoop];
				val2 = grid[start + offset * stride + poleLoop + 1];
				val3 = grid[start + offset * stride + poleLoop + 2];
				val4 = grid[start + offset * stride + poleLoop + 3];

				grid[start + offset * stride + poleLoop] = val1 - right1;
				grid[start + offset * stride + poleLoop + 1] = val2 - right2;
				grid[start + offset * stride + poleLoop + 2] = val3 - right3;
				grid[start + offset * stride + poleLoop + 3] = val4 - right4;
			} // end poleLoop
			offset += stepsize;

			for (int counter = 1; counter < steps - 1; counter++) {
				for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) {
					left1 = grid[start + offset * stride - parOffsetStrided + poleLoop] * 0.5;
					left2 = grid[start + offset * stride - parOffsetStrided + poleLoop + 1] * 0.5;
					left3 = grid[start + offset * stride - parOffsetStrided + poleLoop + 2] * 0.5;
					left4 = grid[start + offset * stride - parOffsetStrided + poleLoop + 3] * 0.5;

					val1 = grid[start + offset * stride + poleLoop];
					val2 = grid[start + offset * stride + poleLoop + 1];
					val3 = grid[start + offset * stride + poleLoop + 2];
					val4 = grid[start + offset * stride + poleLoop + 3];

					right1 = grid[start + offset * stride + poleLoop] * 0.5;
					right2 = grid[start + offset * stride + poleLoop + 1] * 0.5;
					right3 = grid[start + offset * stride + poleLoop + 2] * 0.5;
					right4 = grid[start + offset * stride + poleLoop + 3] * 0.5;

					grid[start + offset * stride + poleLoop] = val1	- left1 - right1;
					grid[start + offset * stride + poleLoop + 1] = val2	- left2 - right2;
					grid[start + offset * stride + poleLoop + 2] = val3 - left3 - right3;
					grid[start + offset * stride + poleLoop + 3] = val4 - left4 - right4;
				}

				offset += stepsize;
			}

			for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) {
				left1 = grid[start + offset * stride - parOffsetStrided + poleLoop] * 0.5;
				left2 = grid[start + offset * stride - parOffsetStrided + poleLoop + 1] * 0.5;
				left3 = grid[start + offset * stride - parOffsetStrided + poleLoop + 2] * 0.5;
				left4 = grid[start + offset * stride - parOffsetStrided + poleLoop + 3] * 0.5;

				val1 = grid[start + offset * stride + poleLoop];
				val2 = grid[start + offset * stride + poleLoop + 1];
				val3 = grid[start + offset * stride + poleLoop + 2];
				val4 = grid[start + offset * stride + poleLoop + 3];

				grid[start + offset * stride + poleLoop] = val1 - left1;
				grid[start + offset * stride + poleLoop + 1] = val2 - left2;
				grid[start + offset * stride + poleLoop + 2] = val3 - left3;
				grid[start + offset * stride + poleLoop + 3] = val4 - left4;
			}

			steps = steps >> 1;
			offset = myPow2(levels[dimension] - level) - 1;
			parentOffset = stepsize;
			stepsize = stepsize << 1;
		} // end loop over levels

		// level = 2 seperate
		for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) {
			right1 = grid[start + (offset + parentOffset) * stride + poleLoop] * 0.5;
			right2 = grid[start + (offset + parentOffset) * stride + poleLoop + 1] * 0.5;
			right3 = grid[start + (offset + parentOffset) * stride + poleLoop + 2] * 0.5;
			right4 = grid[start + (offset + parentOffset) * stride + poleLoop + 3] * 0.5;

			val1 = grid[start + offset * stride + poleLoop];
			val2 = grid[start + offset * stride + poleLoop + 1];
			val3 = grid[start + offset * stride + poleLoop + 2];
			val4 = grid[start + offset * stride + poleLoop + 3];

			grid[start + offset * stride + poleLoop] = val1 - right1;
			grid[start + offset * stride + poleLoop + 1] = val2	- right2;
			grid[start + offset * stride + poleLoop + 2] = val3	- right3;
			grid[start + offset * stride + poleLoop + 3] = val4	- right4;
		}
		offset += stepsize;

		for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) { 
			left1 = grid[start + (offset - parentOffset) * stride + poleLoop] * 0.5;
			left2 = grid[start + (offset - parentOffset) * stride + poleLoop + 1] * 0.5;
			left3 = grid[start + (offset - parentOffset) * stride + poleLoop + 2] * 0.5;
			left4 = grid[start + (offset - parentOffset) * stride + poleLoop + 3] * 0.5;

			val1 = grid[start + offset * stride + poleLoop];
			val2 = grid[start + offset * stride + poleLoop + 1];
			val3 = grid[start + offset * stride + poleLoop + 2];
			val4 = grid[start + offset * stride + poleLoop + 3];

			grid[start + offset * stride + poleLoop] = val1 - left1;
			grid[start + offset * stride + poleLoop + 1] = val2 - left2;
			grid[start + offset * stride + poleLoop + 2] = val3 - left3;
			grid[start + offset * stride + poleLoop + 3] = val4 - left4;
		} // end PoleLoop for level 2
	}

	public void hierarchizeOptimized() {
		int start;
		int stride = 1;
		int pointsInDimension;
		int numberOfPoles;
		int jump;
		pointsInDimension = noAligned;
		numberOfPoles = arraySize / pointsInDimension;

		for (int i = 0; i < numberOfPoles; i++){
			start = i * pointsInDimension;
			hierarchize1DUnoptimized(start, 1, pointsInDimension, 0);
		}

		for (int dimension = 1; dimension < dimensions; dimension++) { // hierarchize d >=1
			stride *= pointsInDimension;
			pointsInDimension = pointsPerDimension[dimension];
			jump = stride * pointsInDimension;
			int blockSize = stride;
			numberOfPoles = arraySize / pointsInDimension;// do loop over first dim in 1d Parts

			for (int i = 0; i < numberOfPoles; i += blockSize){ // integer operations form bottleneck here -- nested loops are twice as slow
				int div = i / stride;
				start = div * jump + (i % stride);
				hierarchize1DOptimized(start, stride, pointsInDimension, dimension, blockSize);
			}
		}
	}

	public void hierarchizeOptimizedThreads(final int numberOfThreads) {
		int dimension;
		int stride = 1;
		int pointsInDimension;
		int polesPerThread;
		int numberOfPoles;
		int numberOfBlocks;
		int blocksPerThread;
		Thread[] threads = new Thread[numberOfThreads];


		//dimension 1 separate as start of each pole is easier to calculate
		pointsInDimension = noAligned;
		numberOfPoles = arraySize / pointsInDimension;
		polesPerThread = numberOfPoles / numberOfThreads;
		for(int i = 0; i < numberOfThreads; i++) {
			final int n = pointsInDimension;
			final int from = i * polesPerThread;
			final int to = (i + 1 == numberOfThreads) ? numberOfPoles : polesPerThread * (i + 1); 
			threads[i] = new Thread(new Runnable() { public void run() {
				for (int j = from; j < to; j++) {
					int start = j * n;
					hierarchize1DUnoptimized(start, 1, n, 0);
				}// end dimension 1
			}});
		}
		for(Thread t : threads)
			t.start();
		for(Thread t : threads)
			try { t.join(); } catch (InterruptedException e) {}


		for(dimension = 1; dimension < dimensions; dimension++) { // hierarchize all dimensions
			stride *= pointsInDimension;
			pointsInDimension = pointsPerDimension[dimension];
			final int jump = stride * pointsInDimension;
			numberOfPoles = arraySize / pointsInDimension;
			final int blockSize = stride;
			numberOfBlocks = numberOfPoles / blockSize;
			blocksPerThread = numberOfBlocks / numberOfThreads;
			polesPerThread = numberOfPoles / numberOfThreads;
			for(int i = 0; i < numberOfThreads; i++) {
				final int s = stride;
				final int d = dimension;
				final int n = pointsInDimension;
				//final int from = i * polesPerThread;
				final int from = i * blocksPerThread;
				//final int to = (i + 1 == numberOfThreads) ? numberOfPoles : polesPerThread * (i + 1);
				final int to = (i + 1 == numberOfThreads) ? numberOfBlocks : blocksPerThread * (i + 1);
				threads[i] = new Thread(new Runnable() { public void run() {
					for (int j = from; j < to; j++) { //loops over blocks, inside loop we multiply loop variable with blockSize
						int k = j * blockSize;
						int div = k / s;
						int start = div * jump + (k % s);
						hierarchize1DOptimized(start, s, n, d, blockSize);
					}
				}});
			} // end loop over dimension 2 to d
			for(Thread t : threads)
				t.start();
			for(Thread t : threads)
				try { t.join(); } catch (InterruptedException e) {}
		}
	}

	public void hierarchizeOptimizedThreadsOnce(final int blockSize, final int numberOfThreads) {
		final CyclicBarrier barrier = new CyclicBarrier(numberOfThreads);
		Thread[] threads = new Thread[numberOfThreads];

		for(int k = 0; k < numberOfThreads; k++) {
			final int i = k;
			threads[i] = new Thread(new Runnable() { public void run() {
				int dimension;
				int stride = 1;
				int pointsInDimension;
				int polesPerThread;
				int numberOfBlocks;
				int blocksPerThread;
				int numberOfPoles;
				int jump;
				int from, to;
				//dimension 1 separate as start of each pole is easier to calculate
				pointsInDimension = noAligned;
				numberOfPoles = arraySize / pointsInDimension;
				polesPerThread = numberOfPoles / numberOfThreads;
				from = i * polesPerThread;
				to = (i + 1 == numberOfThreads) ? numberOfPoles : polesPerThread * (i + 1);
				for (int j = from; j < to; j++) {
					int start = j * pointsInDimension;
					hierarchize1DUnoptimized(start, 1, pointsInDimension, 0);
				}// end dimension 1

				//Wait for all threads to be done with the first dimension
				try { barrier.await(); } catch (InterruptedException | BrokenBarrierException e) {}


				for(dimension = 1; dimension < dimensions; dimension++) { // hierarchize all dimensions
					stride *= pointsInDimension;
					pointsInDimension = pointsPerDimension[dimension];
					jump = stride * pointsInDimension;
					numberOfPoles = arraySize / pointsInDimension;
					//polesPerThread = numberOfPoles / numberOfThreads;
					numberOfBlocks = numberOfPoles / blockSize;
					blocksPerThread = numberOfBlocks / numberOfThreads;
					from = i * blocksPerThread;
					to = (i + 1 == numberOfThreads) ? numberOfBlocks : blocksPerThread * (i + 1);
					for (int j = from; j < to; j++) { // integer operations form bottleneck here -- nested loops are twice as slow
						int m = j * blockSize;
						int div = m / stride;
						int start = div * jump + (m % stride);
						hierarchize1DOptimized(start, stride, pointsInDimension, dimension, blockSize);
					}

					//Wait for all threads to be done with this dimension
					try { barrier.await(); } catch (InterruptedException | BrokenBarrierException e) {}
				} // end loop over dimension 2 to d
			}});
		}
		for(Thread t : threads)
			t.start();
		for(Thread t : threads)
			try { t.join(); } catch (InterruptedException e) {}
	}

	public void hierarchizeOptimizedTasks(final int blockSize, int numberOfTasks) {
		int dimension;
		int stride = 1;
		int pointsInDimension;
		int polesPerTask;
		int numberOfPoles;
		int numberOfBlocks;
		int blocksPerTask;
		ExecutorService executor = Executors.newFixedThreadPool(numberOfTasks);
		//ExecutorService executor = Executors.newWorkStealingPool();
		List<Future<?>> futures = new ArrayList<Future<?>>();

		//dimension 1 separate as start of each pole is easier to calculate
		pointsInDimension = noAligned;
		numberOfPoles = arraySize / pointsInDimension;
		polesPerTask = numberOfPoles / numberOfTasks;
		for(int i = 0; i < numberOfTasks; i++) {
			final int n = pointsInDimension;
			final int from = i * polesPerTask;
			final int to = (i + 1 == numberOfTasks) ? numberOfPoles : polesPerTask * (i + 1); 
			futures.add(executor.submit(new Runnable() { public void run() {
				for (int j = from; j < to; j++) {
					int start = j * n;
					hierarchize1DUnoptimized(start, 1, n, 0);
				}// end dimension 1
			}}));
		}

		try { for (Future<?> fut : futures) fut.get(); } catch (Exception e) {}
		futures.clear();

		for(dimension = 1; dimension < dimensions; dimension++) { // hierarchize all dimensions
			stride *= pointsInDimension;
			pointsInDimension = pointsPerDimension[dimension];
			final int jump = stride * pointsInDimension;
			numberOfPoles = arraySize / pointsInDimension;
			numberOfBlocks = numberOfPoles / blockSize;
			blocksPerTask = numberOfBlocks / numberOfTasks;
			//polesPerTask = numberOfPoles / numberOfTasks;
			for(int i = 0; i < numberOfTasks; i++) {
				final int s = stride;
				final int d = dimension;
				final int n = pointsInDimension;
				final int from = i * blocksPerTask;
				final int to = (i + 1 == numberOfTasks) ? numberOfBlocks : blocksPerTask * (i + 1);
				futures.add(executor.submit(new Runnable() { public void run() {
					for (int j = from; j < to; j++) { // integer operations form bottleneck here -- nested loops are twice as slow
						int k = j * blockSize;
						int div = k / s;
						int start = div * jump + (k % s);
						hierarchize1DOptimized(start, s, n, d, blockSize);
					}
				}}));
			} // end loop over dimension 2 to d
			try { for (Future<?> fut : futures) fut.get(); } catch (Exception e) {}
			futures.clear();
		}
	}



	//	public void hierarchizeOptimizedParallelStream(int blockSize, int numberOfChunks) {
	//		int dimension;
	//		int stride = 1;
	//		int pointsInDimension;
	//		int numberOfPoles;
	//		int jump;
	//		int polesPerChunk;
	//		int numberOfBlocks;
	//		int blocksPerChunk;
	//		List<PoleBlockAligned> blocks = new ArrayList<PoleBlockAligned>();
	//
	//
	//		//dimension 1 separate as start of each pole is easier to calculate
	//		pointsInDimension = noAligned;
	//		numberOfPoles = arraySize / pointsInDimension;
	//		polesPerChunk = numberOfPoles / numberOfChunks;
	//
	//		for(int i = 0; i < numberOfChunks; i++) {
	//			int from = i * polesPerChunk;
	//			int to = (i + 1 == numberOfChunks) ? numberOfPoles : polesPerChunk * (i + 1);
	//			blocks.add(new PoleBlockAligned(this, 1, pointsInDimension, 0, pointsInDimension, blockSize, from, to));
	//		}
	//
	//		blocks
	//		.parallelStream()
	//		.forEach(block -> block.hierarchize());
	//
	//		for(dimension = 1; dimension < dimensions; dimension++) { // hierarchize all dimensions
	//			blocks.clear();
	//			stride *= pointsInDimension;
	//			pointsInDimension = pointsPerDimension[dimension];
	//			jump = stride * pointsInDimension;
	//			numberOfPoles = arraySize / pointsInDimension;
	//			numberOfBlocks = numberOfPoles / blockSize;
	//			blocksPerChunk = numberOfBlocks / numberOfChunks;
	//			for(int i = 0; i < numberOfChunks; i++) {
	//				int from = i * blocksPerChunk;
	//				int to = (i + 1 == numberOfChunks) ? numberOfBlocks : blocksPerChunk * (i + 1);
	//				blocks.add(new PoleBlockAligned(this, stride, pointsInDimension, dimension, jump, blockSize, from, to));
	//			}
	//
	//			blocks
	//			.parallelStream()
	//			.forEach(block -> block.hierarchize());
	//		}
	//		// end loop over dimension 2 to d
	//	}

	//	public void hierarchizeOptimizedParallelStream(int blockSize, int numberOfChunks) {
	//		int dimension;
	//		int stride = 1;
	//		int pointsInDimension;
	//		int numberOfPoles;
	//		int jump;
	//		int polesPerChunk;
	//		int numberOfBlocks;
	//		int blocksPerChunk;
	//		List<PoleBlockAligned> blocks = new ArrayList<PoleBlockAligned>();
	//
	//
	//		//dimension 1 separate as start of each pole is easier to calculate
	//		pointsInDimension = noAligned;
	//		numberOfPoles = arraySize / pointsInDimension;
	//		polesPerChunk = numberOfPoles / numberOfChunks;
	//
	//		for(int i = 0; i < numberOfChunks; i++) {
	//			int from = i * polesPerChunk;
	//			int to = (i + 1 == numberOfChunks) ? numberOfPoles : polesPerChunk * (i + 1);
	//			blocks.add(new PoleBlockAligned(this, 1, pointsInDimension, 0, pointsInDimension, blockSize, from, to));
	//		}
	//
	//		blocks
	//		.parallelStream()
	//		.forEach(block -> block.hierarchize());
	//
	//		for(dimension = 1; dimension < dimensions; dimension++) { // hierarchize all dimensions
	//			blocks.clear();
	//			stride *= pointsInDimension;
	//			pointsInDimension = pointsPerDimension[dimension];
	//			jump = stride * pointsInDimension;
	//			numberOfPoles = arraySize / pointsInDimension;
	//			numberOfBlocks = numberOfPoles / blockSize;
	//			blocksPerChunk = numberOfBlocks / numberOfChunks;
	//			for(int i = 0; i < numberOfChunks; i++) {
	//				int from = i * blocksPerChunk;
	//				int to = (i + 1 == numberOfChunks) ? numberOfBlocks : blocksPerChunk * (i + 1);
	//				blocks.add(new PoleBlockAligned(this, stride, pointsInDimension, dimension, jump, blockSize, from, to));
	//			}
	//
	//			blocks
	//			.parallelStream()
	//			.forEach(block -> block.hierarchize());
	//		}
	//		// end loop over dimension 2 to d
	//	}

	private int myPow2(int i) {
		return 1 << i;
	}

	private int pos( int index[]){
		int retPos = 0;
		for(int i =0; i< dimensions; i++) {
			assert( index[i] < pointsPerDimension[i] );
			assert( 0 <= index[i] );
			retPos += index[i]*strides[i];
		}
		return retPos;
	};


	//METHODS RELATED TO RECURSION.
	double stencil( double center, double left, double right) { return center - .5*left -.5* right;};

	public void hierarchizeRecursive(){ //Overall recursive call
		// this method starts the recursion, using the hierarchizeRec-call.

		int centerInd[] = new int[dimensions];
		Content fullInterval = new Content();
		for(int i = 0; i < dimensions; i++) {
			fullInterval.l[i] = levels[i] - 1; // boundary need not be split away
			centerInd[i] = myPow2(levels[i] - 1) - 1;
		}
		
		fullInterval.l[6] = 0 ; // no predecessors to the left
		fullInterval.l[7] = 0 ; // no predecessors to the right
		int center = pos(centerInd);
		hierarchizeRec(0, dimensions, center, fullInterval, 1);

	}

	private Content copyContent(Content inputContent){
		Content outputContent = new Content(inputContent.asInt, inputContent.l); 
		return outputContent;
	}

	public void hierarchizeRec(int s, int t, int center, Content inputContent, int level){
		//This is the recursive code. This method calls itself, and the hierarchizeApplyStencil4v4, when divided completely.

		Content ic = copyContent(inputContent);

		// chosedim
		int localSize = 0; // sum of levels
		for (int i = 1; i < dimensions; i++) {
			if(ic.l[i] > 0) {
				localSize += ic.l[i];
			}
		}

		if(localSize == 0) { // singleton cache line
			if(ic.l[0] <= 0) { // real singleton
				for(int i = s; i < t; i++) { 
					int rmask = myPow2(i);
					int dist = myPow2(-ic.l[i]);
					double lVal, rVal;
					int posLeft, posRight;
					if((ic.l[6] & rmask) !=0) { //Checks if the bitwise combination equals to 1.
						posLeft = center - dist*strides[i];
						lVal = grid[posLeft];
					}

					else {
						posLeft = -1;
						lVal = 0.0;
					}
					
					if((ic.l[7] & rmask) !=0) { //centerInd[i] + dist < n[i] ) { // it should be == n[i], but hey
						posRight = center + dist*strides[i];
						rVal = grid[posRight];
					}

					else {
						posRight = -1;
						rVal = 0.0;
					}
					
					grid[center] = stencil(grid[center], lVal, rVal);
				}
			} 
			
			else {
				if( s == 0 ) { // actually hierarchize in dir 0
					int rmask = (1 << 0); // replace by iterative?
					int dist = myPow2(ic.l[0]);
					double leftBdVal, rightBdVal;
					int leftBdPos, rightBdPos;
					// hierarchize1DUnoptimized(CGIndex start, CGIndex stride, CGIndex size, int dim) does not fit because it never uses boundary
					// if we don't split in dim0, we know we are at the boundary...
					if((ic.l[6] & rmask) !=0) { //centerInd[i] - dist >= 0 ) { // it should be == -1, but hey
						leftBdPos = center - dist;
						leftBdVal = grid[leftBdPos];
					}
					
					else {
						leftBdPos = -1;
						leftBdVal = 0.0;
					}
					
					if((ic.l[7] & rmask) !=0) { //centerInd[i] + dist < n[i] ) { // it should be == n[i], but hey
						rightBdPos = center + dist;
						rightBdVal = grid[rightBdPos];
					}
					
					else {
						rightBdPos = -1;
						rightBdVal = 0.0;
					}
					
					int step = 1;
					while(step < dist) {
						int start = center - dist + step;

						grid[(start)] = stencil(grid[start], leftBdVal, grid[(start + step)]);

						start += 2*step;
						while( start < center + dist - step ) {

							grid[start] = stencil(grid[start],grid[start-step],grid[start+step]);
							start += 2*step;
						}
						
						assert( start == center+dist-step );
						grid[start] = stencil(grid[start], grid[start-step], rightBdVal);
						step *= 2;
					}
					// while of levels
					grid[(center)] = stencil(grid[center], leftBdVal, rightBdVal);
					s = 1; // hierarchized in dim 0
				}
				
				int d0dist = myPow2(ic.l[0]);
				int first = - d0dist +1;
				int last = + d0dist -1;
				for(int dim=s; dim<t; dim++) {
					int rmask = (1 << dim); // replace by iterative?
					int dist = myPow2(-ic.l[dim]);
					assert(0== (center+first) %4 );
					assert(2== (center+last) %4 );
					if(((ic.l[6] & rmask) !=0) && ((ic.l[7] & rmask) !=0)) {
						for(int i = first; i <= last - 3; i += 4) {
							//System.out.println("Center: " + center + '\t' + "i: " + i + '\t' + "Right + Left");
							hierarchizeApplyStencil4v4(center+i, dist*strides[dim],true,true,dim);
						}
						 //System.out.println("Center: " + center + '\t' + "i: " + (last-2) + '\t' + "Right + Left");
						hierarchizeApplyStencil3v4(center+last-2, dist*strides[dim],true,true,dim);
					}

					if(((ic.l[6] & rmask) !=0) && ((ic.l[7] & rmask) ==0)) {
						for(int i=first; i<= last-3; i+= 4) {
							//System.out.println("Center: " + center + '\t' + "i: " + i + '\t' + "Right");
							hierarchizeApplyStencil4v4(center+i, dist*strides[dim],true,false,dim);
						}
						//System.out.println("Center: " + center + '\t' + "i: " + (last-2) + '\t' + "Right");
						hierarchizeApplyStencil3v4(center+last-2, dist*strides[dim],true,false,dim);
					}
					if(((ic.l[6] & rmask) == 0) && ((ic.l[7] & rmask) !=0)) {
						for(int i=first; i<= last-3; i+= 4) {
							//System.out.println("Center: " + center + '\t' + "i: " + i + '\t' + "Left");
							hierarchizeApplyStencil4v4(center+i, dist*strides[dim],false,true,dim);
						}
						//System.out.println("Center: " + center + '\t' + "i: " + (last-2) + '\t' + "Left");
						hierarchizeApplyStencil3v4(center+last-2, dist*strides[dim],false,true,dim);
					} 
				}
			}
		}
		
		else { 
			int r = 0;
			//We added the cast to int in the following line.
			int maxl = ic.l[0] - recTile - ((int) (recTallPar * localSize)); // block size of pseudo singletons, tall cache assumption. 
			for(int i = 1; i < dimensions; i++){
				if(rf[i] * ic.l[i] > maxl) {
					r = i;
					maxl = (int) (rf[i] * ic.l[i]);
				}
			}
			//System.out.println();
			//System.out.println("r: " + r);
			//ic.printInterval();
			//System.out.println();
			// ic used for right
			Content midI, leftI;
			midI = copyContent(ic);
			midI.l[r] = -midI.l[r];
			ic.l[r]--;
			int dist = myPow2(ic.l[r]); // already reduced!
			//leftI.asInt = ic.asInt;
			leftI = copyContent(ic);
			int rmask = myPow2(r);
			ic.l[6] |= rmask;
			leftI.l[7] |= rmask;
			dist *= strides[r]; // already reduced!
			if(r < s) r = s; // avoid calls!
			if(r > t) r = t;
			if((localSize >= recMinSpawn) && (localSize <= recMaxSpawn) && (r != 0))
			{
				if(r > s) {System.out.println("s: " + s + '\t' + "r: " + r + '\t' + "center: " + center + '\t' + "Mid"); hierarchizeRec(s, r, center, midI, level + 1); }
				System.out.println("s: " + s + '\t' + "r: " + r + '\t' + "center: " + (center-dist) + '\t' + "Left");
				hierarchizeRec(s, t, center - dist, leftI, level + 1);
				System.out.println("s: " + s + '\t' + "r: " + r + '\t' + "center: " + (center+dist) + '\t' + "Right");
				hierarchizeRec(s, t, center + dist, ic, level + 1);
				if(t > r) {System.out.println("s: " + s + '\t' + "r: " + r + '\t' + "center: " + center + '\t' + "Mid"); hierarchizeRec(r, t, center, midI, level + 1); }
			}
			
			else {
				if(r > s) {/*System.out.println("s: " + s + '\t' + "r: " + r + '\t' + "center: " + center + '\t' + "Mid" + '\t' + "Level: " + level);*/ hierarchizeRec(s, r, center, midI, level + 1); }
				//System.out.println("s: " + s + '\t' + "r: " + r + '\t' + "center: " + (center-dist) + '\t' + "Left" + '\t' + "Level: " + level);
				hierarchizeRec(s, t, center - dist, leftI, level + 1);
				//System.out.println("s: " + s + '\t' + "r: " + r + '\t' + "center: " + (center+dist) + '\t' + "Right" + '\t' + "Level: " + level);
				hierarchizeRec(s, t, center + dist, ic, level + 1);
				if(t > r) {/*System.out.println("r: " + r + '\t' + "t: " + t + '\t' + "center: " + center + '\t' + "Mid" + '\t' + "Level: " + level);*/ hierarchizeRec(r, t, center, midI, level + 1); }
			}
		}
	}
	
	public void hierarchizeRecursiveThreads(int NumberOfMaxThreads){ //Overall threaded recursive call
		// this method starts the recursion, using the hierarchizeRec-call.
		//It uses Javas inbuilt threaded recursion-system, which maintains a threadpool.
		ForkJoinPool pool = new ForkJoinPool(NumberOfMaxThreads); //For starting and distributing the tasks.
		
		
		
		int centerInd[] = new int[dimensions];
		Content fullInterval = new Content();
		for(int i =0; i < dimensions; i++) {
			fullInterval.l[i] = levels[i] - 1; // boundary need not be split away
			centerInd[i] = myPow2(levels[i] - 1) - 1;
		}
		
		fullInterval.l[6] = 0 ; // no predecessors to the left
		fullInterval.l[7] = 0 ; // no predecessors to the right
		int center = pos(centerInd);
		
		hierarchizeRecThreads HT = new hierarchizeRecThreads(0, dimensions, center, fullInterval);
		pool.invoke(HT);
	}
	
	
		
	class hierarchizeRecThreads extends RecursiveAction{
		private static final long serialVersionUID = 1L; //Needed for serialization internally in the Java forkjoin-framwork.
		//This is the recursive code. This method calls itself, and the hierarchizeApplyStencil4v4, when divided completely.
		Content ic;
		int s;
		int t;
		int center;
		public hierarchizeRecThreads(int sInput, int tInput, int centerInput, Content inputContent){
			ic = copyContent(inputContent);
			t=tInput;
			s=sInput;
			center=centerInput;
		};
		
		@Override
		protected void compute(){ //This is the content of the recursive call.

			// chosedim
			int localSize = 0; // sum of levels
			for (int i = 1; i < dimensions; i++) {
				if(ic.l[i] > 0) {
					localSize += ic.l[i];
				}
			}

			if(localSize == 0) { // singleton cache line
				if(ic.l[0] <= 0) { // real singletons, t, center, inputContent
					for(int i = s; i < t; i++) { 
						int rmask = myPow2(i);
						int dist = myPow2(-ic.l[i]);
						double lVal, rVal;
						int posLeft, posRight;
						if((ic.l[6] & rmask) !=0) { //Checks if the bitwise combination equals to 1.
							posLeft = center - dist*strides[i];
							lVal = grid[posLeft];
						}

						else {
							posLeft = -1;
							lVal = 0.0;
						}

						if((ic.l[7] & rmask) !=0) { //centerInd[i] + dist < n[i] ) { // it should be == n[i], but hey
							posRight = center + dist*strides[i];
							rVal = grid[posRight];
						}

						else {
							posRight = -1;
							rVal = 0.0;
						}

						grid[center] = stencil(grid[center], lVal, rVal);
					}
				} 

				else {
					if( s == 0 ) { // actually hierarchize in dir 0
						int rmask = (1 << 0); // replace by iterative?
						int dist = myPow2(ic.l[0]);
						double leftBdVal, rightBdVal;
						int leftBdPos, rightBdPos;
						// hierarchize1DUnoptimized(CGIndex start, CGIndex stride, CGIndex size, int dim) does not fit because it never uses boundary
						// if we don't split in dim0, we know we are at the boundary...
						if((ic.l[6] & rmask) !=0) { //centerInd[i] - dist >= 0 ) { // it should be == -1, but hey
							leftBdPos = center - dist;
							leftBdVal = grid[leftBdPos];
						}

						else {
							leftBdPos = -1;
							leftBdVal = 0.0;
						}

						if((ic.l[7] & rmask) !=0) { //centerInd[i] + dist < n[i] ) { // it should be == n[i], but hey
							rightBdPos = center + dist;
							rightBdVal = grid[rightBdPos];
						}

						else {
							rightBdPos = -1;
							rightBdVal = 0.0;
						}

						int step = 1;
						while(step < dist) {
							int start = center - dist + step;

							grid[(start)] = stencil(grid[start], leftBdVal, grid[(start + step)]);

							start += 2*step;
							while( start < center + dist - step ) {

								grid[start] = stencil(grid[start],grid[start-step],grid[start+step]);
								start += 2*step;
							}

							assert( start == center+dist-step );
							grid[start] = stencil(grid[start], grid[start-step], rightBdVal);
							step *= 2;
						}
						// while of levels
						grid[(center)] = stencil(grid[center], leftBdVal, rightBdVal);
						s = 1; // hierarchized in dim 0
					}

					int d0dist = myPow2(ic.l[0]);
					int first = - d0dist +1;
					int last = + d0dist -1;
					for(int dim=s; dim<t; dim++) {
						int rmask = (1 << dim); // replace by iterative?
						int dist = myPow2(-ic.l[dim]);
						assert(0== (center+first) %4 );
						assert(2== (center+last) %4 );
						if(((ic.l[6] & rmask)) !=0 && (ic.l[7] & rmask)!=0) {
							for(int i = first; i <= last - 3; i += 4) {
								hierarchizeApplyStencil4v4(center+i, dist*strides[dim],true,true,dim);
							}
							hierarchizeApplyStencil3v4(center+last-2, dist*strides[dim],true,true,dim);
						}

						if( (ic.l[6] & rmask)!=0 && (ic.l[7] & rmask)==0 ) {
							for(int i=first; i<= last-3; i+= 4) {
								hierarchizeApplyStencil4v4(center+i, dist*strides[dim],true,false,dim);
							}
							hierarchizeApplyStencil3v4(center+last-2, dist*strides[dim],true,false,dim);
						}
						if( (ic.l[6] & rmask)==0 && (ic.l[7] & rmask)!=0 ) {
							for(int i=first; i<= last-3; i+= 4) {
								hierarchizeApplyStencil4v4(center+i, dist*strides[dim],false,true,dim);
							}

							hierarchizeApplyStencil3v4(center+last-2, dist*strides[dim],false,true,dim);
						} 
					}
				}
			}

			else { 
				int r=0;
				//We added the cast to int in the following line.
				int maxl=ic.l[0] - recTile - ((int) (recTallPar * localSize)); // block size of pseudo singletons, tall cache assumption. 
				for(int i=1;i<dimensions;i++){
					if( rf[i]*ic.l[i] > maxl ) {
						r = i;
						maxl = (int) (rf[i]*ic.l[i]);
					}
				}
				// ic used for right
				Content midI, leftI;
				midI = copyContent(ic);
				midI.l[r] = -midI.l[r];
				ic.l[r]--;
				int dist = myPow2(ic.l[r]); // already reduced!
				//leftI.asInt = ic.asInt;
				leftI = copyContent(ic);
				int rmask = myPow2(r);
				ic.l[6] |= rmask;
				leftI.l[7] |= rmask;
				dist *= strides[r]; // already reduced!
				if(r < s) r = s; // avoid calls!
				if(r > t) r = t;
				if((localSize >= recMinSpawn) && (localSize <= recMaxSpawn) && (r != 0))
				{ //It doesn't seem necessary to have this if/else?

					final int sFin = s;
					final int rFin = r;
					final int tFin = t;
					final int centerFin = center;
					final int distFin=dist;
					final Content midFin = midI;
					final Content leftFin = leftI;
					final Content icFin = ic;
					hierarchizeRecThreads HT1=null;
					hierarchizeRecThreads HT4=null;
					if(r > s) {
						HT1= new hierarchizeRecThreads(sFin, rFin, centerFin, midFin); 
					}

					hierarchizeRecThreads HT2 = new hierarchizeRecThreads(sFin, tFin, centerFin - distFin, leftFin);
					hierarchizeRecThreads HT3= new hierarchizeRecThreads(sFin, tFin, centerFin + distFin, icFin);

					if(t > r) {
						HT4= new hierarchizeRecThreads(rFin, tFin, centerFin, midFin);
					}

					if (r>s) HT1.fork(); //Put at the end. Note that Java does not currently seem to support optimization of tailrecursion.
					HT2.fork();
					HT3.fork();
					if (t>r) HT4.fork();
					if (r>s) HT1.join();
					HT2.join();
					HT3.join();
					if (t>r) HT4.join();


				} //Recursive threading stops here. The following lines are for the last recursions.

				else {
					if(r > s) hierarchizeRec(s, r, center, midI, 0);
					hierarchizeRec(s, t, center - dist, leftI, 0);
					hierarchizeRec(s, t, center + dist, ic, 0);
					if(t > r) hierarchizeRec(r, t, center, midI, 0);
				}
			}
		}
	}
	
	public void hierarchizeApplyStencil4v4(int center, int offset, boolean left, boolean right, int r ) {
		double val1, val2, val3, val4, temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, temp4 = 0.0,
				right1 = 0.0, right2 = 0.0, right3 = 0.0, right4 = 0.0, left1 = 0.0, left2 = 0.0, left3 = 0.0, left4 = 0.0;
		int pLeft = center - offset;
		int pRight = center + offset;

		if(!left)  pLeft  = -9;
		if(!right) pRight = -9;

		if(left || right) {
			val1 = grid[center];
			val2 = grid[center + 1];
			val3 = grid[center + 2];
			val4 = grid[center + 3];

			if(left) {
				left1 = grid[pLeft];
				left2 = grid[pLeft + 1];
				left3 = grid[pLeft + 2];
				left4 = grid[pLeft + 3];

				if(right) {
					left1 = left1 * 0.5;
					left2 = left2 * 0.5;
					left3 = left3 * 0.5;
					left4 = left4 * 0.5;
				}

				else {
					temp1 = left1 * 0.5;
					temp2 = left2 * 0.5;
					temp3 = left3 * 0.5;
					temp4 = left4 * 0.5;
				}
			}

			if(right) {
				right1 = grid[pRight];
				right2 = grid[pRight + 1];
				right3 = grid[pRight + 2];
				right4 = grid[pRight + 3];

				if(left) {
					right1 = right1 * 0.5;
					right2 = right2 * 0.5;
					right3 = right3 * 0.5;
					right4 = right4 * 0.5;
				}

				else {
					temp1 = right1 * 0.5;
					temp2 = right2 * 0.5;
					temp3 = right3 * 0.5;
					temp4 = right4 * 0.5;
				}
			}

			if(left &&  right) {
				temp1 = right1 + left1;
				temp2 = right2 + left2;
				temp3 = right3 + left3;
				temp4 = right4 + left4;
			}

			grid[center] = val1 - temp1;
			grid[center + 1] = val2 - temp2;
			grid[center + 2] = val3 - temp3;
			grid[center + 3] = val4 - temp4;
		}
	}

	public void hierarchizeApplyStencil3v4(int center, int offset, boolean left, boolean right, int r ) {
		double val1, val2, val3, temp1 = 0.0, temp2 = 0.0, temp3 = 0.0,
				right1 = 0.0, right2 = 0.0, right3 = 0.0, left1 = 0.0, left2 = 0.0, left3 = 0.0;
		int pLeft = center - offset;
		int pRight = center + offset;

		//Unsure about this, if left == '-' pLeft is never used...
		if(!left)  pLeft  = -9;
		if(!right) pRight = -9;

		if(left || right) {
			val1 = grid[center];
			val2 = grid[center + 1];
			val3 = grid[center + 2];

			if(left) {
				left1 = grid[pLeft];
				left2 = grid[pLeft + 1];
				left3 = grid[pLeft + 2];

				if(right) {
					left1 = left1 * 0.5;
					left2 = left2 * 0.5;
					left3 = left3 * 0.5;
				}

				else {
					temp1 = left1 * 0.5;
					temp2 = left2 * 0.5;
					temp3 = left3 * 0.5;
				}
			}

			if(right) {
				right1 = grid[pRight];
				right2 = grid[pRight + 1];
				right3 = grid[pRight + 2];

				if(left) {
					right1 = right1 * 0.5;
					right2 = right2 * 0.5;
					right3 = right3 * 0.5;
				}

				else {
					temp1 = right1 * 0.5;
					temp2 = right2 * 0.5;
					temp3 = right3 * 0.5;
				}
			}

			if(left &&  right) {
				temp1 = right1 + left1;
				temp2 = right2 + left2;
				temp3 = right3 + left3;
			}

			grid[center] = val1 - temp1;
			grid[center + 1] = val2 - temp2;
			grid[center + 2] = val3 - temp3;
		}
	}
}

class PoleBlockAligned {
	CombiGridAligned grid;
	int stride;
	int pointsInDimension;
	int dimension;
	int blockSize;
	int jump;
	int from, to;

	public PoleBlockAligned(CombiGridAligned grid, int stride, int pointsInDimension, int dimension, int jump, int blockSize, int from, int to) {
		this.grid = grid;
		this.stride = stride;
		this.pointsInDimension = pointsInDimension;
		this.dimension = dimension;
		this.blockSize = blockSize;
		this.jump = jump;
		this.from = from;
		this.to = to;
	}

	public void hierarchize() {
		if(dimension == 0) {
			for (int i = from; i < to; i++) { 
				int div = i / stride;
				int start = div * jump + (i % stride);
				grid.hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
			}
		}
		else {
			for (int i = from; i < to; i++) {
				int k = i * blockSize;
				int div = k / stride;
				int start = div * jump + (k % stride);
				grid.hierarchize1DOptimized(start, stride, pointsInDimension, dimension, blockSize);
			}
		}
	}
}
