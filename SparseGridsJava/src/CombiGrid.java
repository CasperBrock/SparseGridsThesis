import java.util.Arrays;


public class CombiGrid {

	int[] levels;
	int dimensions;
	double[] grid;
	int gridSize;
	int[] pointsPerDimension;

	public static void main(String[] args) {
		int[] levels = {2, 2, 2, 2, 2};
		CombiGrid grid = new CombiGrid(5, levels);
		Arrays.fill(grid.grid, 1.0);
		//grid.hierarchizeUnoptimized();
		grid.hierarchizeUnoptimizedThreads(8);
		grid.printValues();
	}

	public CombiGrid(int dimensions, int[] levels) {
		this.dimensions = dimensions;
		pointsPerDimension = new int[dimensions];
		this.levels = new int[dimensions];
		int size = 1;

		for(int i = 0; i < dimensions; i++) {
			this.levels[i] = levels[i];
			pointsPerDimension[i] = myPow2(levels[i]) - 1;
			size *= pointsPerDimension[i];
		}
		gridSize = size;
		grid = new double[gridSize];
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
		if (dimensions <= 2) printValues2DArr(gridSize, 0, pointsPerDimension[0]);
		else {
			int[] currentLevels = new int[dimensions];
			int chunkSize = pointsPerDimension[0] * pointsPerDimension[1];
			int numberOfChunks = gridSize / chunkSize;
			for (int ctr = 0; ctr < numberOfChunks; ctr++){
				printValues2DArr(chunkSize, ctr*chunkSize, pointsPerDimension[0]);
				currentLevels[2]++;
				for (int dd = 2; dd < dimensions; dd++){
					if (currentLevels[dd] == pointsPerDimension[dd]) {
						currentLevels[dd] = 0;
						if(dd + 1 < dimensions)
							currentLevels[dd + 1]++;
					}
				}
			}
		}
	}

	/***
	 * Hierarchizes a single dimension with an unoptimized algorithm
	 * 
	 * @param start Where in the grid to start
	 * @param stride The distance between points to work on
	 * @param size How many points to hierachize
	 * @param dimension The dimension to work in
	 */
	public void hierarchize1DUnoptimized(int start, int stride, int size, int dimension) {
		int level, steps, ctr, offset, parentOffset, stepSize, parOffsetStrided;
		double val1, val2, val3, left, right;

		level = levels[dimension];
		steps = myPow2(level - 1);
		offset = 0;
		stepSize = 2;
		parentOffset = 1;
		left = 0;
		right = 0;

		for(level--; level > 1; level--) {
			parOffsetStrided = parentOffset*stride;
			grid[start + offset * stride] -= 0.5 * grid[start + offset * stride + parOffsetStrided];
			offset += stepSize;
			left = 0.5 * grid[start + offset * stride - parOffsetStrided];
			for (ctr = 1; ctr < steps - 1; ctr++) {
				val1 = grid[start+offset*stride];
				right = 0.5 * grid[start+offset*stride+parOffsetStrided];
				val2 = val1 - left;
				val3 = val2 - right;
				grid[start+offset*stride] = (int) val3;
				left = right;
				offset += stepSize;
			} 

			grid[start + offset * stride] -= right;
			//steps = steps >> 1;
			steps = steps / 2;
			offset = myPow2(levels[dimension] - level) - 1;
			parentOffset =  stepSize;
			//stepSize = stepSize << 1;
			stepSize = stepSize * 2;
		}

		right = 0.5 * grid[start + (offset + parentOffset) * stride];
		grid[start+offset*stride] -= right;
		offset += stepSize;
		grid[start+offset*stride] -= right;
	}

	/***
	 * Hierachizes the grid with a sequential unoptimized algorithm
	 */
	public void hierarchizeUnoptimized() {
		int dimension;
		int start;
		int stride = 1;
		int pointsInDimension;
		int numberOfPoles;
		int jump;


		//dimension 1 separate as start of each pole is easier to calculate
		pointsInDimension = pointsPerDimension[0];
		numberOfPoles = gridSize / pointsInDimension;
		for (int i = 0; i < numberOfPoles; i++){
			start = i * pointsInDimension;
			hierarchize1DUnoptimized(start, 1, pointsInDimension, 0);
		}
		// end dimension 1

		for(dimension = 1; dimension < dimensions; dimension++){ // hierarchize all dimensions
			stride *= pointsInDimension;
			pointsInDimension = pointsPerDimension[dimension];
			jump = stride * pointsInDimension;
			numberOfPoles = gridSize / pointsInDimension;
			for (int i = 0; i < numberOfPoles; i++){ // integer operations form bottleneck here -- nested loops are twice as slow
				int div = i / stride;
				start = div * jump + (i % stride);
				hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
			}
		} // end loop over dimension 2 to d
	}

	public void hierarchizeUnoptimizedThreads(int numberOfThreads) {
		int dimension;
		int stride = 1;
		final int polesPerThread;
		Thread[] threads = new Thread[numberOfThreads];


		//dimension 1 separate as start of each pole is easier to calculate
		final int pointsInDimension = pointsPerDimension[0];
		final int numberOfPoles = gridSize / pointsInDimension;
		polesPerThread = numberOfPoles / numberOfThreads;
		for(int i = 0; i < numberOfThreads; i++) {
			final int k = i;
			threads[i] = new Thread(new Runnable() { public void run() {
				if(k + 1 == numberOfThreads) {
					for (int j = k * polesPerThread; j < numberOfPoles; j++) {
						int start = j * pointsInDimension;
						hierarchize1DUnoptimized(start, 1, pointsInDimension, 0);
					}
				}
				else {
					for (int j = k * polesPerThread; j < (k + 1) * polesPerThread; j++) {
						int start = j * pointsInDimension;
						hierarchize1DUnoptimized(start, 1, pointsInDimension, 0);
					}
				}// end dimension 1
			}});
		}
		for(Thread t : threads)
			t.start();
		for(Thread t : threads)
			try { t.join(); } catch (InterruptedException e) {}


		for(dimension = 1; dimension < dimensions; dimension++) { // hierarchize all dimensions
			stride *= pointsInDimension;
			final int pointsInDimension2 = pointsPerDimension[dimension];
			final int jump = stride * pointsInDimension2;
			final int numberOfPoles2 = gridSize / pointsInDimension2;
			final int polesPerThread2 = numberOfPoles2 / numberOfThreads;
			for(int i = 0; i < numberOfThreads; i++) {
				final int k = i;
				final int strides = stride;
				final int d = dimension;
				threads[i] = new Thread(new Runnable() { public void run() {
					if(k + 1 == numberOfThreads) {
						for(int j = k * polesPerThread2; j < numberOfPoles2; j++) {
							int div = j / strides;
							int start = div * jump + (j % strides);
							hierarchize1DUnoptimized(start, strides, pointsInDimension2, d);
						}
					}
					else {
						for (int j = k * polesPerThread2; j < (k + 1) * polesPerThread2; j++) { // integer operations form bottleneck here -- nested loops are twice as slow
							int div = j / strides;
							int start = div * jump + (j % strides);
							hierarchize1DUnoptimized(start, strides, pointsInDimension2, d);
						}
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
