package testing;

import grid.*;
import gridFunctions.*;

public class Verification {

	public static void main(String[] args) {
		//verifyUnoptimized();
		verifyOptimized();
	}

	public static void verifyUnoptimized() {
		boolean test1, test2, test3;
		CombiGrid real, test;

		int[] levels1 = {2, 2};
		int[] levels2 = {5, 5, 5, 5, 5};
		int[] levels3 = {10, 5, 4, 3, 2};

		real = new CombiGrid(levels1);
		real.setValues(GridFunctions.ALLONES);
		real.hierarchizeUnoptimized();
		test = new CombiGrid(levels1);
		test1 = true;
		for(int i = 1; i < 16; i++) {
			test.setValues(GridFunctions.ALLONES);
			test.hierarchizeOptimized();
			test1 = test1 && real.compare(test);
		}

		real = new CombiGrid(levels2);
		real.setValues(GridFunctions.ALLONES);
		real.hierarchizeUnoptimized();
		test = new CombiGrid(levels2);
		test2 = true;
		for(int i = 1; i < 16; i++) {
			test.setValues(GridFunctions.ALLONES);
			test.hierarchizeOptimized();
			test2 = test2 && real.compare(test);
		}

		real = new CombiGrid(levels3);
		real.setValues(GridFunctions.ALLONES);
		real.hierarchizeUnoptimized();
		test = new CombiGrid(levels3);
		test3 = true;
		for(int i = 1; i < 16; i++) {
			test.setValues(GridFunctions.ALLONES);
			test.hierarchizeOptimized();
			test3 = test3 && real.compare(test);
		}

		System.out.println("Result: " + (test1 && test2 && test3));
	}
	
	public static void verifyOptimized() {
		boolean test1, test2, test3;
		CombiGridAligned real, test;

		int[] levels1 = {2, 2};
		int[] levels2 = {5, 5, 5, 5, 5};
		int[] levels3 = {10, 5, 4, 3, 2};

		real = new CombiGridAligned(levels1, 32);
		real.setValues(GridFunctions.ALLONES);
		real.hierarchizeOptimized();
		test = new CombiGridAligned(levels1, 32);
		test1 = true;
		for(int i = 1; i < 16; i++) {
			test.setValues(GridFunctions.ALLONES);
			test.hierarchizeOptimizedThreadsOnce(4, i);
			test1 = test1 && real.compare(test);
		}

		real = new CombiGridAligned(levels2, 32);
		real.setValues(GridFunctions.ALLONES);
		real.hierarchizeOptimized();
		test = new CombiGridAligned(levels2, 32);
		test2 = true;
		for(int i = 1; i < 16; i++) {
			test.setValues(GridFunctions.ALLONES);
			test.hierarchizeOptimizedThreadsOnce(4, i);
			test2 = test2 && real.compare(test);
		}

		real = new CombiGridAligned(levels3, 32);
		real.setValues(GridFunctions.ALLONES);
		real.hierarchizeOptimized();
		test = new CombiGridAligned(levels3, 32);
		test3 = true;
		for(int i = 1; i < 16; i++) {
			test.setValues(GridFunctions.ALLONES);
			test.hierarchizeOptimizedThreadsOnce(4, i);
			test3 = test3 && real.compare(test);
		}

		System.out.println("Result: " + (test1 && test2 && test3));
	}

}
