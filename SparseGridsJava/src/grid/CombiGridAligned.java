package grid;

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

	public static void main(String[] args) {
		int[] levels = {10, 5, 4, 3, 2};
		CombiGridAligned grid = new CombiGridAligned(levels, 32);
		CombiGridAligned grid2 = new CombiGridAligned(levels, 32);
		System.out.println("Gridsize: " + grid.gridSize);
		System.out.println("Arraysize: " + grid.arraySize);
		grid.setValues(GridFunctions.ALLONES);
		grid.hierarchizeOptimized(16);
		//grid.printValues();
		grid2.setValues(GridFunctions.ALLONES);
		grid2.hierarchizeOptimizedThreads(16, 400);
		//grid.printValues();
		System.out.println(grid2.compare(grid));
	}

	CombiGridAligned(int[] levels, int alignment) {
		// alignment in bytes
		// alignment multiple of 32 bytes for AVX
		// grid needs to be aligned for the blocked (optimized) version of the
		// code
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
		for (i = 1; i < dimensions; i++) {
			this.levels[i] = levels[i];
			pointsPerDimension[i] = myPow2(levels[i]) - 1;
			gridSize *= pointsPerDimension[i];
			arraySize *= pointsPerDimension[i];
		}

		grid = new double[arraySize];
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

	private void hierarchize1DUnoptimized(int start, int stride, int size, int dimension) {
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

	private void hierarchize1DOptimized(int start, int stride, int size,
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
				right2 = grid[start + offset * stride + poleLoop	+ 1] * 0.5;
				right3 = grid[start + offset * stride + poleLoop + 2] * 0.5;
				right4 = grid[start + offset * stride + poleLoop	+ 3] * 0.5;

				val1 = grid[start + offset * stride + poleLoop];
				val2 = grid[start + offset * stride + poleLoop+1];
				val3 = grid[start + offset * stride + poleLoop+2];
				val4 = grid[start + offset * stride + poleLoop+3];

				grid[start + offset * stride + poleLoop] = val1 - right1;
				grid[start + offset * stride + poleLoop + 1] = val2 - right2;
				grid[start + offset * stride + poleLoop + 2] = val3 - right3;
				grid[start + offset * stride + poleLoop + 3] = val4 - right4;
			} // end poleLoop
			offset += stepsize;

			for (int counter = 1; counter < steps - 1; counter++) {
				for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) {
					left1 = grid[start + offset * stride	- parOffsetStrided + poleLoop] * 0.5;
					left2 = grid[start + offset * stride	- parOffsetStrided + poleLoop + 1] * 0.5;
					left3 = grid[start + offset * stride	- parOffsetStrided + poleLoop + 2] * 0.5;
					left4 = grid[start + offset * stride	- parOffsetStrided + poleLoop + 3] * 0.5;

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
				val2 = grid[start + offset * stride + poleLoop+1];
				val3 = grid[start + offset * stride + poleLoop+2];
				val4 = grid[start + offset * stride + poleLoop+3];

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

	public void hierarchizeOptimized(int blockSize) {
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
			numberOfPoles = arraySize / pointsInDimension;// do loop over first dim in 1d Parts

			for (int i = 0; i < numberOfPoles; i += blockSize){ // integer operations form bottleneck here -- nested loops are twice as slow
				int div = i / stride;
				start = div * jump + (i % stride);
				hierarchize1DOptimized(start, stride, pointsInDimension, dimension, blockSize);
			}
		}
	}

	public void hierarchizeOptimizedThreads(int blockSize, int numberOfThreads) {
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

	private int myPow2(int i) {
		return 1 << i;
	}
}
