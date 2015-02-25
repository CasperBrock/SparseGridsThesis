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
		int[] levels = { 2, 2, 2, 2 };
		CombiGridAligned grid = new CombiGridAligned(levels, 32);
		System.out.println("Gridsize: " + grid.gridSize);
		System.out.println("Arraysize: " + grid.arraySize);
		grid.setValues(GridFunctions.ALLONES);
		grid.printValues();
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
	
	private void hierarchizeOptimized(int blockSize){
		int dim;
		int start;
		int stride =1 ;
		int ndim;
		int nbrOfPoles;
		int jump;
		int divresult_qoutient;
		int divresult_remainder;
		ndim = noAligned;
		nbrOfPoles = arraySize /ndim;

		for (int kk = 0; kk < nbrOfPoles; kk++){
			start = kk *ndim;
			hierarchize1DUnoptimized(start, 1, ndim,0);
			}
		
		for (dim = 1; dim < dimensions; dim ++){ // hierarchize d >=1
			stride *=ndim;
			ndim = pointsPerDimension[dim];
			jump = stride*ndim;
			nbrOfPoles = arraySize/ndim;// do loop over first dim in 1d Parts
			
			for (int nn = 0; nn< nbrOfPoles; nn+=blockSize){ // integer operations form bottleneck here -- nested loops are twice as slow
				divresult_qoutient = nn/stride;
				divresult_remainder = nn % stride;
			start = divresult_qoutient*jump +divresult_remainder;
			hierarchize1DOptimized(start, stride, ndim,dim, blockSize);
			}
			}
		
	}

	private void hierarchize1DOptimized(int start, int stride, int size,
			int dim, int unroll) {
		// optimized with vector operations. lvl for any dim must be larger
		// than 1.
		int ll;
		int steps;
		int ctr;
		int offset, parentOffset;
		int stepsize;
		int parOffsetStrided;
		ll = levels[dim];
		steps = myPow2(ll - 1);
		offset = 0;
		stepsize = 2;
		parentOffset = 1;

		for (ll--; ll > 1; ll--) {
			parOffsetStrided = parentOffset * stride;
			for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) { // multiply
																		// right
																		// parents
																		// by
																		// 0.5.
				double parRight1_05 = grid[start + offset * stride + poleLoop] * 0.5;
				double parRight2_05 = grid[start + offset * stride + poleLoop
						+ 1] * 0.5;
				double parRight3_05 = grid[start + offset * stride + poleLoop
						+ 2] * 0.5;
				double parRight4_05 = grid[start + offset * stride + poleLoop
						+ 3] * 0.5;

				grid[start + offset * stride + poleLoop] = parRight1_05;
				grid[start + offset * stride + poleLoop + 1] = parRight2_05;
				grid[start + offset * stride + poleLoop + 2] = parRight3_05;
				grid[start + offset * stride + poleLoop + 3] = parRight4_05;
			} // end poleLoop
			offset += stepsize;
			for (ctr = 1; ctr < steps - 1; ctr++) {
				for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) { // subtract
																			// 0.5*rightpar
																			// and
																			// 0.5*leftpar
																			// from
																			// org.
					double parLeft1_05 = grid[start + offset * stride
							- parOffsetStrided + poleLoop] * 0.5;
					double parLeft2_05 = grid[start + offset * stride
							- parOffsetStrided + poleLoop + 1] * 0.5;
					double parLeft3_05 = grid[start + offset * stride
							- parOffsetStrided + poleLoop + 2] * 0.5;
					double parLeft4_05 = grid[start + offset * stride
							- parOffsetStrided + poleLoop + 3] * 0.5;

					double Org1 = grid[start + offset * stride + poleLoop];
					double Org2 = grid[start + offset * stride + poleLoop + 1];
					double Org3 = grid[start + offset * stride + poleLoop + 2];
					double Org4 = grid[start + offset * stride + poleLoop + 3];

					double parRight1_05 = grid[start + offset * stride
							+ poleLoop] * 0.5;
					double parRight2_05 = grid[start + offset * stride
							+ poleLoop + 1] * 0.5;
					double parRight3_05 = grid[start + offset * stride
							+ poleLoop + 2] * 0.5;
					double parRight4_05 = grid[start + offset * stride
							+ poleLoop + 3] * 0.5;

					grid[start + offset * stride + poleLoop] = Org1
							- parLeft1_05 - parRight1_05;
					grid[start + offset * stride + poleLoop + 1] = Org2
							- parLeft2_05 - parRight2_05;
					;
					grid[start + offset * stride + poleLoop + 2] = Org3
							- parLeft3_05 - parRight3_05;
					;
					grid[start + offset * stride + poleLoop + 3] = Org4
							- parLeft4_05 - parRight4_05;
					;
				}
				offset += stepsize;
			}
			for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) {
				double parLeft1_05 = grid[start + offset * stride
						- parOffsetStrided + poleLoop] * 0.5;
				double parLeft2_05 = grid[start + offset * stride
						- parOffsetStrided + poleLoop + 1] * 0.5;
				double parLeft3_05 = grid[start + offset * stride
						- parOffsetStrided + poleLoop + 2] * 0.5;
				double parLeft4_05 = grid[start + offset * stride
						- parOffsetStrided + poleLoop + 3] * 0.5;

				grid[start + offset * stride + poleLoop] = parLeft1_05;
				grid[start + offset * stride + poleLoop + 1] = parLeft2_05;
				grid[start + offset * stride + poleLoop + 2] = parLeft3_05;
				grid[start + offset * stride + poleLoop + 3] = parLeft4_05;
			}
			steps = steps >> 1;
			offset = myPow2(levels[dim] - ll) - 1;
			parentOffset = stepsize;
			stepsize = stepsize << 1;
		} // end loop over levels
			// level = 2 seperate
		for (int poleLoop = 0; poleLoop < unroll; poleLoop += 4) {
			double parRight1_05 = grid[start + (offset + parentOffset) * stride
					+ poleLoop] * 0.5;
			double parRight2_05 = grid[start + (offset + parentOffset) * stride
					+ poleLoop + 1] * 0.5;
			double parRight3_05 = grid[start + (offset + parentOffset) * stride
					+ poleLoop + 2] * 0.5;
			double parRight4_05 = grid[start + (offset + parentOffset) * stride
					+ poleLoop + 3] * 0.5;

			double org1_R = grid[start + offset * stride + poleLoop];
			double org2_R = grid[start + offset * stride + poleLoop + 1];
			double org3_R = grid[start + offset * stride + poleLoop + 2];
			double org4_R = grid[start + offset * stride + poleLoop + 3];

			grid[start + offset * stride + poleLoop] = org1_R - parRight1_05;
			grid[start + offset * stride + poleLoop + 1] = org2_R
					- parRight2_05;
			grid[start + offset * stride + poleLoop + 2] = org3_R
					- parRight3_05;
			grid[start + offset * stride + poleLoop + 3] = org4_R
					- parRight4_05;

			offset += stepsize;
			for (poleLoop = 0; poleLoop < unroll; poleLoop += 4) { //c++ let's you redeclare int poleloop. We reuse the value instead, since the first loop is done.
				double parLeft1_05 = grid[start + (offset - parentOffset)
						* stride + poleLoop] * 0.5;
				double parLeft2_05 = grid[start + (offset - parentOffset)
						* stride + poleLoop + 1] * 0.5;
				double parLeft3_05 = grid[start + (offset - parentOffset)
						* stride + poleLoop + 2] * 0.5;
				double parLeft4_05 = grid[start + (offset - parentOffset)
						* stride + poleLoop + 3] * 0.5;

				double org1_L = grid[start + offset * stride + poleLoop];
				double org2_L = grid[start + offset * stride + poleLoop + 1];
				double org3_L = grid[start + offset * stride + poleLoop + 2];
				double org4_L = grid[start + offset * stride + poleLoop + 3];

				grid[start + offset * stride + poleLoop] = org1_L - parLeft1_05;
				grid[start + offset * stride + poleLoop + 1] = org2_L
						- parLeft2_05;
				grid[start + offset * stride + poleLoop + 2] = org3_L
						- parLeft3_05;
				grid[start + offset * stride + poleLoop + 3] = org4_L
						- parLeft4_05;
			} // end PoleLoop for level 2
		}
	}
	
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

	private int myPow2(int i) {
		return 1 << i;
	}
}
