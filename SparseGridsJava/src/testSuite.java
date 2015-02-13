import java.util.Arrays;

public class testSuite {

	/**
	 * @param args
	 */

	static long nanoStart;
	static long nanoStop;
	static long currentTime;
	static CombiGrid CG;
	static int time;
	static int k; //so we can call it from outside its loop.
	static int best; //this will hold the current best time for each test.
	static int bestThreads;

	//This class can run various tests on the CombiGrid.java code.
	//supposed to test amount of cores and threads, threading methods and time. 
	public static void main(String[] args) {
		int cores = Runtime.getRuntime().availableProcessors();
		long mem = Runtime.getRuntime().maxMemory()/(1024*1024);
		System.out.println("Running tests on "+cores+" cores using up to "+mem+" MB of memory.");

		for (int i = 2; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<6;j++){// this loop will relate to number of levels.
				bestThreads=0;
				int best = Integer.MAX_VALUE;
				for (k = 1; k<cores*3;k++){//this loop iterates over number of threads.
					System.out.println("Test with " + i + " dimensions, all with level " + j + ", running in "+k+" threads.");
					time = 0;
					for (int l = 0;l<10;l++){

						int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

						//lines that create the objects we are to work on. Should not be timed, ofc.
						Arrays.fill(levels, j);
						CG = new CombiGrid(i, levels);
						Arrays.fill(CG.grid, 1);					
						nanoStart = System.nanoTime();
						CG.hierarchizeUnoptimizedThreads(k);
						//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
						nanoStop = System.nanoTime();

						currentTime = nanoStop-nanoStart; 		
						time += currentTime;
					}
					
					if (time/10 < best) {
						best=time/10;
						bestThreads=k;
					}
					System.out.println("Test D" + i + "-L" + j + "-T"+k+" took " + (time/10) + " nanoseconds to complete in average." );
					System.out.println("");
				}
				System.out.println("Best number of threads for D"+i+"-L"+j+" is "+bestThreads+" in " + best + " nanoseconds.");
				System.out.println("");
			}
		}		
	}
}