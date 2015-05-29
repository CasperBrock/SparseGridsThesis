package grid;

import gridFunctions.*;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/***
 * Class containing the simple component grid and all the hierarchization methods
 * run on the standard grid.
 */
public class CombiGrid {

	public int[] levels;
	public int dimensions;
	public double[] grid;
	public int gridSize;
	public int[] pointsPerDimension;

	/**
	 * Create a new CombiGrid based on the given level vector.
	 * @param levels The level vector describing the structure of the grid.
	 */
	public CombiGrid(int[] levels) {
		dimensions = levels.length;
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

	/**
	 * Prints the values of the grid to the console.
	 */
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

	/**
	 * Returns a string of the level vector for the grid as well as the size of the vector
	 * 
	 * @return String describing the level vector of the grid
	 */
	public String getLevels() {
		StringBuilder s = new StringBuilder();
		int sum = 0;
		for(int i : levels) {
			sum += i;
			s.append("" +  i + " ");
		}
		s.append("[" + sum + "]");
		return s.toString();
	}

	/**
	 * Compares this grid to another CombiGrid. Returns whether or not the grids
	 * contain the same values.
	 * @param cg Grid to compare with.
	 * @return True if both grids have the same values for all grid points, false otherwise.
	 */
	public boolean compare(CombiGrid cg) {
		if (cg.gridSize != gridSize)
			return false;
		for(int i = 0; i < gridSize; i++)
			if(Math.abs(cg.grid[i] - grid[i]) > 0.000001) {
				System.out.println(cg.grid[i]);
				System.out.println(grid[i]);
				return false;
			}
		return true;
	}

	/**
	 * Sets the values of the grid according to a given function.
	 * 
	 * @param func The function to set the grid values.
	 */
	public void setValues(GridFunctions func) {
		double[] stepsize = new double[dimensions];
		double[] x = new double[dimensions];
		int[] levelSets = new int[dimensions];
		levelSets[0] = 1;
		int[] dimensionCounter = new int[dimensions];

		for (int d = 1; d < dimensions; d++){
			levelSets[d] = levelSets[d - 1] * pointsPerDimension[d - 1];
		}

		for (int d = 0; d < dimensions; d++){
			stepsize[d] = Math.pow(2, -levels[d]);
			dimensionCounter[d] = 1;
		}

		for (int counter = 0; counter < gridSize; counter++){
			for (int d = 0; d < dimensions; d++){
				if (dimensionCounter[d] > pointsPerDimension[d] && d == dimensions - 1) {
					return;
				}
				if (dimensionCounter[d] > pointsPerDimension[d] && d < dimensions - 1) {
					dimensionCounter[d] = 1;
					dimensionCounter[d + 1]++;
				}
			}
			for (int d = 0; d < dimensions; d ++){
				x[d] = dimensionCounter[d] * stepsize[d];
			}

			grid[counter] = GridFunction.call(x, func);
			dimensionCounter[0]++;
		}
		return ;
	}

	/**
	 * Hierarchizes a single dimension in the grid.
	 * 
	 * @param start Where in the grid to start
	 * @param stride The distance between points to work on
	 * @param size How many points to hierachize
	 * @param dimension The dimension to work in
	 */
	private void hierarchize1DUnoptimized(int start, int stride, int size, int dimension) {
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
				val1 = grid[start + offset * stride];
				right = 0.5 * grid[start + offset * stride + parOffsetStrided];
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
		grid[start + offset * stride] -= right;
		offset += stepSize;
		grid[start + offset * stride] -= right;
	}

	/**
	 * Hierarchizes a single dimension in the grid. This method uses loop-unrolling of 4
	 * in an attempt to reach higher performance.
	 * 
	 * @param start Where in the grid to start
	 * @param stride The distance between points to work on
	 * @param size How many points to hierachize
	 * @param dimension The dimension to work in
	 * @param unroll How far we can work in the grid, must be a multiple of 4
	 */
	private void hierarchize1DOptimized4(int start, int stride, int size,
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

	/**
	 * Hierarchizes a single dimension in the grid. This method uses loop-unrolling by 3
	 * in an attempt to reach higher performance.
	 * 
	 * @param start Where in the grid to start
	 * @param stride The distance between points to work on
	 * @param size How many points to hierachize
	 * @param dimension The dimension to work in
	 * @param unroll How far we can go in the grid. Must be a multiple of 3
	 */
	private void hierarchize1DOptimized3(int start, int stride, int size,
			int dimension, int unroll) {
		// optimized with vector operations. lvl for any dim must be larger
		// than 1.
		int level;
		int steps;
		int offset, parentOffset;
		int stepsize;
		int parOffsetStrided;
		double val1, val2, val3, right1, right2, right3, left1, left2, left3;
		level = levels[dimension];
		steps = myPow2(level - 1);
		offset = 0;
		stepsize = 2;
		parentOffset = 1;

		for (level--; level > 1; level--) {
			parOffsetStrided = parentOffset * stride;

			for (int poleLoop = 0; poleLoop < unroll; poleLoop += 3) {
				right1 = grid[start + offset * stride + poleLoop] * 0.5;
				right2 = grid[start + offset * stride + poleLoop + 1] * 0.5;
				right3 = grid[start + offset * stride + poleLoop + 2] * 0.5;

				val1 = grid[start + offset * stride + poleLoop];
				val2 = grid[start + offset * stride + poleLoop + 1];
				val3 = grid[start + offset * stride + poleLoop + 2];

				grid[start + offset * stride + poleLoop] = val1 - right1;
				grid[start + offset * stride + poleLoop + 1] = val2 - right2;
				grid[start + offset * stride + poleLoop + 2] = val3 - right3;
			} // end poleLoop
			offset += stepsize;

			for (int counter = 1; counter < steps - 1; counter++) {
				for (int poleLoop = 0; poleLoop < unroll; poleLoop += 3) {
					left1 = grid[start + offset * stride - parOffsetStrided + poleLoop] * 0.5;
					left2 = grid[start + offset * stride - parOffsetStrided + poleLoop + 1] * 0.5;
					left3 = grid[start + offset * stride - parOffsetStrided + poleLoop + 2] * 0.5;

					val1 = grid[start + offset * stride + poleLoop];
					val2 = grid[start + offset * stride + poleLoop + 1];
					val3 = grid[start + offset * stride + poleLoop + 2];

					right1 = grid[start + offset * stride + poleLoop] * 0.5;
					right2 = grid[start + offset * stride + poleLoop + 1] * 0.5;
					right3 = grid[start + offset * stride + poleLoop + 2] * 0.5;

					grid[start + offset * stride + poleLoop] = val1	- left1 - right1;
					grid[start + offset * stride + poleLoop + 1] = val2	- left2 - right2;
					grid[start + offset * stride + poleLoop + 2] = val3 - left3 - right3;
				}

				offset += stepsize;
			}

			for (int poleLoop = 0; poleLoop < unroll; poleLoop += 3) {
				left1 = grid[start + offset * stride - parOffsetStrided + poleLoop] * 0.5;
				left2 = grid[start + offset * stride - parOffsetStrided + poleLoop + 1] * 0.5;
				left3 = grid[start + offset * stride - parOffsetStrided + poleLoop + 2] * 0.5;

				val1 = grid[start + offset * stride + poleLoop];
				val2 = grid[start + offset * stride + poleLoop + 1];
				val3 = grid[start + offset * stride + poleLoop + 2];

				grid[start + offset * stride + poleLoop] = val1 - left1;
				grid[start + offset * stride + poleLoop + 1] = val2 - left2;
				grid[start + offset * stride + poleLoop + 2] = val3 - left3;
			}

			steps = steps >> 1;
			offset = myPow2(levels[dimension] - level) - 1;
			parentOffset = stepsize;
			stepsize = stepsize << 1;
		} // end loop over levels

		// level = 2 seperate
		for (int poleLoop = 0; poleLoop < unroll; poleLoop += 3) {
			right1 = grid[start + (offset + parentOffset) * stride + poleLoop] * 0.5;
			right2 = grid[start + (offset + parentOffset) * stride + poleLoop + 1] * 0.5;
			right3 = grid[start + (offset + parentOffset) * stride + poleLoop + 2] * 0.5;

			val1 = grid[start + offset * stride + poleLoop];
			val2 = grid[start + offset * stride + poleLoop + 1];
			val3 = grid[start + offset * stride + poleLoop + 2];

			grid[start + offset * stride + poleLoop] = val1 - right1;
			grid[start + offset * stride + poleLoop + 1] = val2	- right2;
			grid[start + offset * stride + poleLoop + 2] = val3	- right3;
		}
		offset += stepsize;

		for (int poleLoop = 0; poleLoop < unroll; poleLoop += 3) { 
			left1 = grid[start + (offset - parentOffset) * stride + poleLoop] * 0.5;
			left2 = grid[start + (offset - parentOffset) * stride + poleLoop + 1] * 0.5;
			left3 = grid[start + (offset - parentOffset) * stride + poleLoop + 2] * 0.5;

			val1 = grid[start + offset * stride + poleLoop];
			val2 = grid[start + offset * stride + poleLoop + 1];
			val3 = grid[start + offset * stride + poleLoop + 2];

			grid[start + offset * stride + poleLoop] = val1 - left1;
			grid[start + offset * stride + poleLoop + 1] = val2 - left2;
			grid[start + offset * stride + poleLoop + 2] = val3 - left3;
		} // end PoleLoop for level 2
	}

	/**
	 * Hierarchizes a single dimension in the grid. This method does not loop unroll, but handles a
	 * large number of poles in a single call.
	 * 
	 * @param start Where in the grid to start
	 * @param stride The distance between points to work on
	 * @param size How many points to hierachize
	 * @param dimension The dimension to work in
	 * @param unroll How far in the grid we can continue
	 */
	private void hierarchize1DOptimizedNoUnroll(int start, int stride, int size,
			int dimension, int unroll) {
		// optimized with vector operations. lvl for any dim must be larger
		// than 1.
		int level;
		int steps;
		int offset, parentOffset;
		int stepsize;
		int parOffsetStrided;
		double right, left, val;
		level = levels[dimension];
		steps = myPow2(level - 1);
		offset = 0;
		stepsize = 2;
		parentOffset = 1;

		for (level--; level > 1; level--) {
			parOffsetStrided = parentOffset * stride;

			for (int poleLoop = 0; poleLoop < unroll; poleLoop += 1) {
				right = grid[start + offset * stride + poleLoop] * 0.5;

				val = grid[start + offset * stride + poleLoop];

				grid[start + offset * stride + poleLoop] = val - right;
			} // end poleLoop
			offset += stepsize;

			for (int counter = 1; counter < steps - 1; counter++) {
				for (int poleLoop = 0; poleLoop < unroll; poleLoop += 1) {
					left = grid[start + offset * stride - parOffsetStrided + poleLoop] * 0.5;

					val = grid[start + offset * stride + poleLoop];

					right = grid[start + offset * stride + poleLoop] * 0.5;

					grid[start + offset * stride + poleLoop] = val	- left - right;
				}

				offset += stepsize;
			}

			for (int poleLoop = 0; poleLoop < unroll; poleLoop += 1) {
				left = grid[start + offset * stride - parOffsetStrided + poleLoop] * 0.5;

				val = grid[start + offset * stride + poleLoop];

				grid[start + offset * stride + poleLoop] = val - left;
			}

			steps = steps >> 1;
			offset = myPow2(levels[dimension] - level) - 1;
			parentOffset = stepsize;
			stepsize = stepsize << 1;
		} // end loop over levels

		// level = 2 seperate
		for (int poleLoop = 0; poleLoop < unroll; poleLoop += 1) {
			right = grid[start + (offset + parentOffset) * stride + poleLoop] * 0.5;

			val = grid[start + offset * stride + poleLoop];

			grid[start + offset * stride + poleLoop] = val - right;
		}
		offset += stepsize;

		for (int poleLoop = 0; poleLoop < unroll; poleLoop += 1) { 
			left = grid[start + (offset - parentOffset) * stride + poleLoop] * 0.5;

			val = grid[start + offset * stride + poleLoop];

			grid[start + offset * stride + poleLoop] = val - left;
		} // end PoleLoop for level 2
	}

	/***
	 * Hierarchizes the grid with a sequential algorithm.
	 * This method performs the entire hierarchization sequentially.
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
			for (int i = 0; i < numberOfPoles; i++) {
				int div = i / stride;
				start = div * jump + (i % stride);
				hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
			}
		} // end loop over dimension 2 to d
	}

	/***
	 * Hierarchizes the grid with a sequential algorithm using loop-unrolling by 4 and 3
	 * to try and reach higher performance.
	 * This method performs the entire hierarchization sequentially.
	 */
	public void hierarchizeOptimized() {
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

		for (dimension = 1; dimension < dimensions; dimension++) { // hierarchize d >=1
			stride *= pointsInDimension;
			pointsInDimension = pointsPerDimension[dimension];
			jump = stride * pointsInDimension;
			numberOfPoles = gridSize / pointsInDimension;
			int blockSize = (stride / 4 * 4 < 4) ? 4 : stride / 4 * 4;
			int blocks = stride / blockSize;
			int rem = stride % blockSize;
			int i = 0;
			start = 0;

			while(i < numberOfPoles) {
				int div = i / stride;
				start = div * jump + (i % stride);

				for(int j = 0; j < blocks; j++) {
					hierarchize1DOptimized4(start, stride, pointsInDimension, dimension, blockSize);
					start += blockSize;
				}

				int j = rem;				
				int block3 = j / 3;
				hierarchize1DOptimized3(start, stride, pointsInDimension, dimension, 3 * block3);
				j = j % 3;
				start += 3 * block3;
				while(j > 0) {
					hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
					start++;
					j--;
				}

				i += stride;
			}
		}
	}

	/***
	 * Hierarchizes the grid with a parallel algorithm.
	 * This method creates new threads for each dimension, resulting in a overhead
	 * on thread creation.
	 * 
	 * @param numberOfThreads The amount of threads to use for the parallelisation
	 */
	public void hierarchizeUnoptimizedThreads(final int numberOfThreads) {
		int dimension;
		int stride = 1;
		int pointsInDimension;
		int polesPerThread;
		int numberOfPoles;
		Thread[] threads = new Thread[numberOfThreads];


		//dimension 1 separate as start of each pole is easier to calculate
		pointsInDimension = pointsPerDimension[0];
		numberOfPoles = gridSize / pointsInDimension;
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
			numberOfPoles = gridSize / pointsInDimension;
			polesPerThread = numberOfPoles / numberOfThreads;
			for(int i = 0; i < numberOfThreads; i++) {
				final int s = stride;
				final int d = dimension;
				final int n = pointsInDimension;
				final int from = i * polesPerThread;
				final int to = (i + 1 == numberOfThreads) ? numberOfPoles : polesPerThread * (i + 1);
				threads[i] = new Thread(new Runnable() { public void run() {
					for (int j = from; j < to; j++) { // integer operations form bottleneck here -- nested loops are twice as slow
						int div = j / s;
						int start = div * jump + (j % s);
						hierarchize1DUnoptimized(start, s, n, d);
					}
				}});
			} // end loop over dimension 2 to d
			for(Thread t : threads)
				t.start();
			for(Thread t : threads)
				try { t.join(); } catch (InterruptedException e) {}
		}
	}

	/***
	 * Hierarchizes the grid with a parallel algorithm. This method uses loop-unrolling by 4 and 3 
	 * to try and reach higher performance.
	 * This method creates new threads for each dimension, resulting in a overhead
	 * on thread creation.
	 * 
	 * @param numberOfThreads The amount of threads to use for the parallelisation
	 */
	public void hierarchizeOptimizedThreads(final int numberOfThreads) {
		int dimension;
		int stride = 1;
		int pointsInDimension;
		int polesPerThread;
		int numberOfPoles;
		Thread[] threads = new Thread[numberOfThreads];


		//dimension 1 separate as start of each pole is easier to calculate
		pointsInDimension = pointsPerDimension[0];
		numberOfPoles = gridSize / pointsInDimension;
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


		for(dimension = 1; dimension < dimensions - 1; dimension++) { // hierarchize from d = 1 to dimensions - 1
			stride *= pointsInDimension;
			pointsInDimension = pointsPerDimension[dimension];
			final int jump = stride * pointsInDimension;
			numberOfPoles = gridSize / pointsInDimension;
			final int jumps = numberOfPoles / stride;
			final int actualThreads;
			if(jumps < numberOfThreads) {
				actualThreads = jumps;
				threads = new Thread[actualThreads];
			}
			else
				actualThreads = numberOfThreads;
			final int jumpsPerThread = jumps / actualThreads;
			for(int i = 0; i < actualThreads; i++) {
				final int blockSize = (stride / 4 * 4 < 4) ? 4 : stride / 4 * 4;
				final int s = stride;
				final int d = dimension;
				final int n = pointsInDimension;
				final int from = i * jumpsPerThread;
				final int to = (i + 1 == actualThreads) ? jumps : jumpsPerThread * (i + 1);
				threads[i] = new Thread(new Runnable() { public void run() {
					//Loop over jumps
					for(int k = from; k < to; k++) {
						int blocks = s / blockSize;
						int rem = s % blockSize;
						int i = k * s;
						int div = i / s;
						int start = div * jump + (i % s);

						for(int j = 0; j < blocks; j++) {
							hierarchize1DOptimized4(start, s, n, d, blockSize);
							start += blockSize;
						}

						int j = rem;				
						int block3 = j / 3;
						if(block3 > 0)
							hierarchize1DOptimized3(start, s, n, d, 3 * block3);
						j = j % 3;
						start += 3 * block3;
						while(j > 0) {
							hierarchize1DUnoptimized(start, s, n, d);
							start++;
							j--;
						}
					}
				}});
			} // end loop over dimension 2 to d-1

			for(Thread t : threads)
				t.start();
			for(Thread t : threads)
				try { t.join(); } catch (InterruptedException e) {}
		}

		threads = new Thread[numberOfThreads];
		//Final dimension has 0 jumps, handle differently
		stride *= pointsInDimension;
		pointsInDimension = pointsPerDimension[dimension];
		final int jump = stride * pointsInDimension;
		numberOfPoles = gridSize / pointsInDimension;
		polesPerThread = numberOfPoles / numberOfThreads;
		for(int i = 0; i < numberOfThreads; i++) {
			final int s = stride;
			final int d = dimension;
			final int n = pointsInDimension;
			final int from = i * polesPerThread;
			final int blockSize = (polesPerThread / 4 * 4 < 4) ? 4 : polesPerThread / 4 * 4;
			final int to = (i + 1 == numberOfThreads) ? numberOfPoles : polesPerThread * (i + 1);
			threads[i] = new Thread(new Runnable() { public void run() {			
				int blocks = (to - from) / blockSize;
				int rem = (to - from) % blockSize;
				int i = from;
				int start = 0;

				while(i < to) {
					int div = i / s;
					start = div * jump + (i % s);

					for(int j = 0; j < blocks; j++) {
						hierarchize1DOptimized4(start, s, n, d, blockSize);
						start += blockSize;
					}

					int j = rem;				
					int block3 = j / 3;
					if(block3 > 0)
						hierarchize1DOptimized3(start, s, n, d, 3 * block3);
					j = j % 3;
					start += 3 * block3;
					while(j > 0) {
						hierarchize1DUnoptimized(start, s, n, d);
						start++;
						j--;
					}

					i += s;
				}
			}});
		}

		for(Thread t : threads)
			t.start();
		for(Thread t : threads)
			try { t.join(); } catch (InterruptedException e) {}
	}

	/**
	 * Hierarchizes the grid using a parallel algorithm. This method eliminates method calls
	 * by handling a large number of poles per call.
	 * 
	 * @param numberOfThreads Number of threads to use for the parallelisation
	 */
	public void hierarchizeOptimizedThreadsNoUnroll(final int numberOfThreads) {
		int dimension;
		int stride = 1;
		int pointsInDimension;
		int polesPerThread;
		int numberOfPoles;
		Thread[] threads = new Thread[numberOfThreads];


		//dimension 1 separate as start of each pole is easier to calculate
		pointsInDimension = pointsPerDimension[0];
		numberOfPoles = gridSize / pointsInDimension;
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


		for(dimension = 1; dimension < dimensions - 1; dimension++) { // hierarchize from d = 1 to dimensions - 1
			stride *= pointsInDimension;
			pointsInDimension = pointsPerDimension[dimension];
			final int jump = stride * pointsInDimension;
			numberOfPoles = gridSize / pointsInDimension;
			final int jumps = numberOfPoles / stride;
			final int actualThreads;
			if(jumps < numberOfThreads) {
				actualThreads = jumps;
				threads = new Thread[actualThreads];
			}
			else
				actualThreads = numberOfThreads;
			final int jumpsPerThread = jumps / actualThreads;
			for(int i = 0; i < actualThreads; i++) {
				final int s = stride;
				final int d = dimension;
				final int n = pointsInDimension;
				final int from = i * jumpsPerThread;
				final int to = (i + 1 == actualThreads) ? jumps : jumpsPerThread * (i + 1);
				threads[i] = new Thread(new Runnable() { public void run() {
					//Loop over jumps
					for(int k = from; k < to; k++) {
						int i = k * s;
						int div = i / s;
						int start = div * jump + (i % s);
						hierarchize1DOptimizedNoUnroll(start, s, n, d, s);
					}
				}});
			} // end loop over dimension 2 to d-1

			for(Thread t : threads)
				t.start();
			for(Thread t : threads)
				try { t.join(); } catch (InterruptedException e) {}
		}

		threads = new Thread[numberOfThreads];
		//Final dimension has 0 jumps, handle differently
		stride *= pointsInDimension;
		pointsInDimension = pointsPerDimension[dimension];
		final int jump = stride * pointsInDimension;
		numberOfPoles = gridSize / pointsInDimension;
		polesPerThread = numberOfPoles / numberOfThreads;
		for(int i = 0; i < numberOfThreads; i++) {
			final int s = stride;
			final int d = dimension;
			final int n = pointsInDimension;
			final int from = i * polesPerThread;
			final int to = (i + 1 == numberOfThreads) ? numberOfPoles : polesPerThread * (i + 1);
			threads[i] = new Thread(new Runnable() { public void run() {
				int i = from;
				int start = 0;
				int div = i / s;
				start = div * jump + (i % s);
				hierarchize1DOptimizedNoUnroll(start, s, n, d, (to - from));
			}});
		}

		for(Thread t : threads)
			t.start();
		for(Thread t : threads)
			try { t.join(); } catch (InterruptedException e) {}
	}

	/***
	 * Hierarchizes the grid using an parallel algorithm using threads.
	 * This method creates the threads only once and runs the looping calculations inside the threads.
	 * Creates less overhead on thread creation, but all loop calculations are done numberOfThreads times. 
	 * 
	 * @param numberOfThreads The number of threads to use for the parallelisation
	 */
	public void hierarchizeUnoptimizedThreadsOnce(final int numberOfThreads) {
		final CyclicBarrier barrier = new CyclicBarrier(numberOfThreads);
		Thread[] threads = new Thread[numberOfThreads];

		for(int k = 0; k < numberOfThreads; k++) {
			final int i = k;
			threads[i] = new Thread(new Runnable() {public void run() {
				int dimension;
				int stride = 1;
				int pointsInDimension;
				int polesPerThread;
				int numberOfPoles;
				int jump;
				int from, to;
				//dimension 1 separate as start of each pole is easier to calculate
				pointsInDimension = pointsPerDimension[0];
				numberOfPoles = gridSize / pointsInDimension;
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
					numberOfPoles = gridSize / pointsInDimension;
					polesPerThread = numberOfPoles / numberOfThreads;
					from = i * polesPerThread;
					to = (i + 1 == numberOfThreads) ? numberOfPoles : polesPerThread * (i + 1);
					for (int j = from; j < to; j++) { // integer operations form bottleneck here -- nested loops are twice as slow
						int div = j / stride;
						int start = div * jump + (j % stride);
						hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
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

	/***
	 * Hierarchizes the grid with a parallel unoptimized algorithm using the Java Task framework.
	 * Creates an overhead in task creation and some in execution, but avoids the overhead of
	 * creating many threads. All thread creation and maintenance is handled by the Java
	 * ExecutorService.
	 * 
	 * @param numberOfTasks The number of tasks to use
	 */
	public void hierarchizeUnoptimizedTasks(final int numberOfTasks) {
		int dimension;
		int stride = 1;
		int pointsInDimension;
		int polesPerTask;
		int numberOfPoles;
		ExecutorService executor = Executors.newCachedThreadPool();
		List<Future<?>> futures = new ArrayList<Future<?>>();

		//dimension 1 separate as start of each pole is easier to calculate
		pointsInDimension = pointsPerDimension[0];
		numberOfPoles = gridSize / pointsInDimension;
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
			numberOfPoles = gridSize / pointsInDimension;
			polesPerTask = numberOfPoles / numberOfTasks;
			for(int i = 0; i < numberOfTasks; i++) {
				final int s = stride;
				final int d = dimension;
				final int n = pointsInDimension;
				final int from = i * polesPerTask;
				final int to = (i + 1 == numberOfTasks) ? numberOfPoles : polesPerTask * (i + 1);
				futures.add(executor.submit(new Runnable() { public void run() {
					for (int j = from; j < to; j++) { // integer operations form bottleneck here -- nested loops are twice as slow
						int div = j / s;
						int start = div * jump + (j % s);
						hierarchize1DUnoptimized(start, s, n, d);
					}
				}}));
			} // end loop over dimension 2 to d
			try { for (Future<?> fut : futures) fut.get(); } catch (Exception e) {}
			futures.clear();
		}

		executor.shutdown();
	}

	/**
	 * Hierarchizes the grid using a sequential algorithm. This method eliminates method calls
	 * by handling a large number of poles per call to hierarchize1DOptimizedNoUnroll
	 */
	public void hierarchizeOptimizedNoUnroll() {
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
			for(int i = 0; i < numberOfPoles; i += stride) {
				int div = i / stride;
				start = div * jump + (i % stride);
				hierarchize1DOptimizedNoUnroll(start, stride, pointsInDimension, dimension, stride);
			}
		} // end loop over dimension 2 to d
	}

	private int myPow2(int i) {
		return 1 << i;
	}
}