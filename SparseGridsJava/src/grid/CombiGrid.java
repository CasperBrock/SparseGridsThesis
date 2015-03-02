package grid;
import gridFunctions.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.BrokenBarrierException;
import java.util.concurrent.Callable;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class CombiGrid {

	public int[] levels;
	public int dimensions;
	public double[] grid;
	public int gridSize;
	public int[] pointsPerDimension;

	public static void main(String[] args) {
		int[] levels = {3, 2, 3, 2, 3};
		CombiGrid grid = new CombiGrid(levels);
		//int[] levels = {3, 3};
		//CombiGrid grid = new CombiGrid(2, levels);
		//Arrays.fill(grid.grid, 1.0);
		grid.setValues(GridFunctions.ALLONES);
		//System.out.println("Array size is: " + grid.grid.length);
		//grid.hierarchizeUnoptimized();
		//grid.hierarchizeUnoptimizedThreads(8);
		//grid.hierarchizeUnoptimizedThreadsOnce(8);
		//grid.hierarchizeUnoptimizedTasks(100);
		//grid.hierarchizeUnoptimizedParallelStream();
		//grid.hierarchizeUnoptimizedParallelStream(100);
		grid.printValues();
	}

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

	/***
	 * Hierarchizes the grid with a sequential unoptimized algorithm.
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
			for (int i = 0; i < numberOfPoles; i++){ // integer operations form bottleneck here -- nested loops are twice as slow
				int div = i / stride;
				start = div * jump + (i % stride);
				hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
			}
		} // end loop over dimension 2 to d
	}

	/***
	 * Hierarchizes the grid with a parallel unoptimized algorithm.
	 * This method creates new threads for each dimension, resulting in a overhead
	 * on thread creation.
	 * 
	 * @param numberOfThreads The amount of threads to use for the parallelization
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
	 * Hierarchizes the grid using an unoptimized algorithm using threads.
	 * This method creates the threads only once and runs the looping calculations inside the threads.
	 * Creates less overhead on thread creation, but all loop calculations are done numberOfThreads times. 
	 * 
	 * @param numberOfThreads The number of threads to use for the parallelization
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
		ExecutorService executor = Executors.newWorkStealingPool();
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
	}

	/***
	 * Hierarchizes the grid using a parallel unoptimized algorithm using the Java Parallel Stream framework.
	 * This method creates a pole object for each pole, and then uses the parallelStream() method to run all the poles
	 * in parallel. All the parallelism is done by the Java framework.
	 */
	public void hierarchizeUnoptimizedParallelStream() {
		int dimension;
		int start;
		int stride = 1;
		int pointsInDimension;
		int numberOfPoles;
		int jump;
		List<Pole> poles = new ArrayList<Pole>();


		//dimension 1 separate as start of each pole is easier to calculate
		pointsInDimension = pointsPerDimension[0];
		numberOfPoles = gridSize / pointsInDimension;
		for (int i = 0; i < numberOfPoles; i++){
			start = i * pointsInDimension;
			poles.add(new Pole(start, 1, pointsInDimension, 0));
		}
		// end dimension 1
		poles
		.parallelStream()
		.forEach(pole -> hierarchize1DUnoptimized(pole.start, pole.stride, pole.pointsInDimension, pole.dimension));

		for(dimension = 1; dimension < dimensions; dimension++){ // hierarchize all dimensions
			poles.clear();
			stride *= pointsInDimension;
			pointsInDimension = pointsPerDimension[dimension];
			jump = stride * pointsInDimension;
			numberOfPoles = gridSize / pointsInDimension;
			for (int i = 0; i < numberOfPoles; i++){ // integer operations form bottleneck here -- nested loops are twice as slow
				int div = i / stride;
				start = div * jump + (i % stride);
				poles.add(new Pole(start, stride, pointsInDimension, dimension));
			}
			poles
			.parallelStream()
			.forEach(pole -> hierarchize1DUnoptimized(pole.start, pole.stride, pole.pointsInDimension, pole.dimension));
		} // end loop over dimension 2 to d
	}

	/***
	 * Hierarchizes the grid using a parallel unoptimized algorithm using the Java Parallel Stream framework.
	 * This method creates numberOfBlocks poleBlock objects per dimension to hierarchize a subset of the poles.
	 * Then is uses the parallelStream() method to run the hierarchization in parallel.
	 * All parallism is handled by the Java framework.
	 *
	 * @param numberOfChunks The number of block objects to use
	 */
	public void hierarchizeUnoptimizedParallelStream(int numberOfChunks) {
		int dimension;
		int stride = 1;
		int pointsInDimension;
		int numberOfPoles;
		int jump;
		int polesPerBlock;
		List<PoleBlock> blocks = new ArrayList<PoleBlock>();


		//dimension 1 separate as start of each pole is easier to calculate
		pointsInDimension = pointsPerDimension[0];
		numberOfPoles = gridSize / pointsInDimension;
		polesPerBlock = numberOfPoles / numberOfChunks;

		for(int i = 0; i < numberOfChunks; i++) {
			int from = i * polesPerBlock;
			int to = (i + 1 == numberOfChunks) ? numberOfPoles : polesPerBlock * (i + 1);
			blocks.add(new PoleBlock(this, 1, pointsInDimension, 0, pointsInDimension, from, to));
		}

		blocks
		.parallelStream()
		.forEach(block -> block.hierarchize());

		for(dimension = 1; dimension < dimensions; dimension++) { // hierarchize all dimensions
			blocks.clear();
			stride *= pointsInDimension;
			pointsInDimension = pointsPerDimension[dimension];
			jump = stride * pointsInDimension;
			numberOfPoles = gridSize / pointsInDimension;
			polesPerBlock = numberOfPoles / numberOfChunks;
			for(int i = 0; i < numberOfChunks; i++) {
				int from = i * polesPerBlock;
				int to = (i + 1 == numberOfChunks) ? numberOfPoles : polesPerBlock * (i + 1);
				blocks.add(new PoleBlock(this, stride, pointsInDimension, dimension, jump, from, to));
			}

			blocks
			.parallelStream()
			.forEach(block -> block.hierarchize());
		}
		// end loop over dimension 2 to d
	}

	private int myPow2(int i) {
		return 1 << i;
	}
}

/***
 * Object to contain the variables needed for a pole.
 * Used in the hierarchizeUnoptimizedParallelStream method.
 */
class Pole {
	public int start;
	public int stride;
	public int pointsInDimension;
	public int dimension;

	public Pole(int start, int stride, int pointsInDimension, int dimension) {
		this.start = start;
		this.stride = stride;
		this.pointsInDimension = pointsInDimension;
		this.dimension = dimension;
	}
}

class PoleBlock {
	CombiGrid grid;
	int stride;
	int pointsInDimension;
	int dimension;
	int jump;
	int from, to;

	public PoleBlock(CombiGrid grid, int stride, int pointsInDimension, int dimension, int jump, int from, int to) {
		this.grid = grid;
		this.stride = stride;
		this.pointsInDimension = pointsInDimension;
		this.dimension = dimension;
		this.jump = jump;
		this.from = from;
		this.to = to;
	}

	public void hierarchize() {
		for (int i = from; i < to; i++) { 
			int div = i / stride;
			int start = div * jump + (i % stride);
			grid.hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
		}
	}
}
