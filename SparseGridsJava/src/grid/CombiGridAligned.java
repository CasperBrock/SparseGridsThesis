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
		int[] levels = {2, 2, 2, 2};
		CombiGridAligned grid = new CombiGridAligned(levels, 32);
		System.out.println("Gridsize: " + grid.gridSize);
		System.out.println("Arraysize: " + grid.arraySize);
		grid.setValues(GridFunctions.ALLONES);
		grid.printValues();
	}
	
	CombiGridAligned(int[] levels, int alignment) {
		// alignment in bytes
		// alignment multiple of 32 bytes for AVX
		//grid needs to be aligned for the blocked (optimized) version of the code
		dimensions = levels.length;
		this.levels = new int[dimensions];
		pointsPerDimension = new int[dimensions];
		int i;
		// lengthen 1st dimensions
		pointsPerDimension[0] = myPow2(levels[0]) - 1;
		this.levels[0]= levels[0];
		noAligned = (int) (Math.ceil((double) (myPow2(levels[0]) - 1) / alignment * 8.0) * alignment / 8.0);
		gridSize = pointsPerDimension[0];
		arraySize = noAligned;
		for (i = 1; i < dimensions; i++){
			this.levels[i]= levels[i];
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
		if (dimensions <= 2) printValues2DArr(arraySize, 0, noAligned);
		else {
			int[] currentLevels = new int[dimensions];
			int chunkSize = noAligned * pointsPerDimension[1];
			int numberOfChunks = arraySize / chunkSize;
			for (int ctr = 0; ctr < numberOfChunks; ctr++){
				printValues2DArr(chunkSize, ctr*chunkSize, noAligned);
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
	
	public void setValues(GridFunctions func) {
		alignment = alignment / 8;
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

		int pos = 0;
		for (int counter = 0; counter < gridSize; counter++){
			for (int d = 0; d < dimensions; d++){
				if (dimensionCounter[d] > pointsPerDimension[d] && d == dimensions - 1) {
					System.out.println("Something went wrong");
					return;
				}
				
				if (d == 0 && dimensionCounter[d] > pointsPerDimension[d]) { // we are in first dim (padded!) and need to increment pos
					grid[pos] = Double.POSITIVE_INFINITY; // pos points to padded point
					pos++;
					dimensionCounter[d] = 1;
					dimensionCounter[d + 1]++;
				}
				
				if (dimensionCounter[d] > pointsPerDimension[d] && d < dimensions - 1) {
					dimensionCounter[d] = 1;
					dimensionCounter[d + 1]++;
				}
			}
			
			for (int d = 0; d < dimensions; d ++){
				x[d] = dimensionCounter[d] * stepsize[d];
			}
			
			grid[pos] = GridFunction.call(x, func);
			dimensionCounter[0]++;
			pos++;
		}
		
		grid[arraySize - 1] = Double.POSITIVE_INFINITY;
		return ;
	}
	
	private int myPow2(int i) {
		return 1 << i;
	}
}
