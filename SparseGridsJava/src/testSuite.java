import java.util.Arrays;

public class testSuite {

	/**
	 * @param args
	 */
	
	static long millisStart;
	static long millisStop;
	
	//This class can run various tests on the CombiGrid.java code.
	//supposed to test amount of cores and threads, threading methods and time. 
	public static void main(String[] args) {
		int cores = Runtime.getRuntime().availableProcessors();
		long mem = Runtime.getRuntime().freeMemory()/(1024*1024);
		System.out.println("Running tests on "+cores+" cores using "+mem+" MB of memory.");
		
		for (int i = 2; i<8; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<5;j++){// this loop will relate to number of levels.
				for (int k = 1; k<cores*2;k++){//TODO this loop will relate to number of threads, when implemented.
				
					int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.
				
					//lines that create the objects we are to work on. Should not be timed, ofc.
					Arrays.fill(levels, j);
					//try {
						CombiGrid CG = new CombiGrid(i, levels);
						
					//} catch (OutOfMemoryError e){
					//	System.out.println("The system ran out of memory running the final test. ");
					//}
					Arrays.fill(CG.grid, 1);
				
					System.out.println("Test with " + i + " dimensions, all with level " + j + ", running in "+k+" threads.");
						
					millisStart = System.currentTimeMillis();
				
					CG.hierarchizeUnoptimizedThreads(k);
		
					//	CG.printValues(); //print only to check that everything is being calculated correctly.
					millisStop = System.currentTimeMillis();
				
								
					System.out.println("Test D" + i + "-L" + j + "-T"+k+" took " + ((millisStop-millisStart)) + " milliseconds to complete." );
				
				}
			}
		}		
	}
}