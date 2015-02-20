package testing;
import grid.CombiGrid;
import gridFunctions.GridFunction;
import gridFunctions.GridFunctions;

import java.util.Arrays;

public class TestSuite {

	/**
	 * @param args
	 */

	static long nanoStart;
	static long nanoStop;
	static long currentTime;
	static long currentTimeStream; //used to hold time for the Streaming method.
	static CombiGrid CG;
	static double time;
	static double timeStream;
	static int k; //so we can call it from outside its loop.
	static double best; //this will hold the current best time for each test.
	static int bestThreads;
	static int cores;
	static long mem;
	static long nanoStartStream;
	static long nanoStopStream;
	static int bestThreadAmount;
	static long bestThreaded;
	static long nanoStartTask;
	static long nanoStopTask;
	static long currentTimeTask;
	static long timeTask;
	static double currentBestTaskTime;
	static int timeStreamBlock;
	static double bestBlockTime;
	static int bestBlocks;

	//This class can run various tests on the CombiGrid.java code.
	//supposed to test amount of cores and threads, threading methods and time. 
	public static void main(String[] args) {
		cores = Runtime.getRuntime().availableProcessors();
		mem = Runtime.getRuntime().maxMemory()/(1024*1024);
		System.out.println("Running tests on "+cores+" cores using up to "+mem+" MB of memory.");
		
		//Which tests to run. The following four will increase grid size gradually, and test the average runtime over 10 iterations.
		//testFullThreads(true); //Note that these are designed to crash, to you can only run one.
		//testFullTasks();
		//testFullPStream();
		//testFullPStreamBlocks(true); //true will print all. False, will only print the result for the best number of blocks.
		//testFullThreadsOnce();
		
		//This test will compare the best results for each method, and only print those.
		//testSideways();
		int[] test = {4,4,4,4,4};
		testAllSimple(test, GridFunctions.ALLONES); //input an array of levels, and a GridFunctions to use.
		
	}
	static void testFullThreads(boolean PrintAll){ //Will test with an different amounts of threads.
		System.out.println("dim\tlev\tthreads/tasks\ttime\ttype");
		
		for (int i = 2; i<8; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<6;j++){// this loop will relate to number of levels.
				bestThreads=0;
				double bestThreadsTime=Integer.MAX_VALUE;
				double best = Integer.MAX_VALUE;
				for (k = 1; k<cores*3;k++){//this loop iterates over number of threads.
					//System.out.println("Test with " + i + " dimensions, all with level " + j + ", running in "+k+" threads.");
					time = 0;
					for (int l = 0;l<10;l++){

						int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

						Arrays.fill(levels, j); //lines that create the objects we are to work on. Should not be timed, ofc.
						CG = new CombiGrid(levels);
						Arrays.fill(CG.grid, 1);					
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedThreads(k);
						//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
						nanoStop = System.currentTimeMillis();

						currentTime = nanoStop-nanoStart; 		
						time += currentTime;
					}
					if (PrintAll) System.out.println(i+"\t"+j+"\t"+k+"\t"+time/10+"\tThreads"); // will print info for Thread
					
					}
				
					if (time/10<bestThreadsTime){
						bestThreadsTime=time/10;
						bestThreads=k;
					
				}
				if (!PrintAll) System.out.println(i+"\t"+j+"\t"+bestThreads+"\t"+bestThreadsTime+"\tThreads");
			}
		}		
	}
	
	//TESTS Threads, where each thread is maintained and handed a new worktask, to minimize overhead.
	static void testFullThreadsOnce(){
		System.out.println("Testing threading with thread re-use.");
		System.out.println("dim\tlev\tthreads/tasks\ttime\ttype");
		
		for (int i = 2; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<7;j++){// this loop will relate to number of levels.
				bestThreads=0;
				double best = Integer.MAX_VALUE;
				for (k = 1; k<cores*3;k++){//this loop iterates over number of threads.
					//System.out.println("Test with " + i + " dimensions, all with level " + j + ", running in "+k+" threads.");
					time = 0;
					for (int l = 0;l<10;l++){

						int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

						Arrays.fill(levels, j); //lines that create the objects we are to work on. Should not be timed, ofc.
						CG = new CombiGrid(levels);
						Arrays.fill(CG.grid, 1);					
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedThreadsOnce(k);
						//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
						nanoStop = System.currentTimeMillis();

						currentTime = nanoStop-nanoStart; 		
						time += currentTime;
					}
					System.out.println(i+"\t"+j+"\t"+k+"\t"+time/10+"\tThreads"); // will print info for Thread
				}
			}
		}		
	}
	
	
	
	static void testFullPStream(){
		System.out.println("dim\tlev\tthreads/tasks\ttime in ms\ttype");
		
		for (int i = 2; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<7;j++){// this loop will relate to number of levels.
				bestThreads=0;
				double best = Integer.MAX_VALUE;
					time = 0;
					for (int l = 0;l<10;l++){

						int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

						Arrays.fill(levels, j); //lines that create the objects we are to work on. Should not be timed, ofc.
						CG = new CombiGrid(levels);
						Arrays.fill(CG.grid, 1);					
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedParallelStream();
						//	CG.printValnew int[10,5,4,3,2]ues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
						nanoStop = System.currentTimeMillis();

						currentTime = nanoStop-nanoStart; 		
						time += currentTime;
					}
					System.out.println(i+"\t"+j+"\t"+"not known"+"\t"+time/10+"\tThreads"); // will print info for Thread
				}
			}
		}		
	
	static void testFullPStreamBlocks(boolean PrintAllResults){ 
		if (!PrintAllResults) System.out.println("Printing only best for dimension/level:");
		
		System.out.println("Dim\tlevel\tblocks\ttime");
		int bestBlocks =0;
		double best = Integer.MAX_VALUE;
		for (int i = 3; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<7;j++){// this loop will relate to number of levels.
				currentTime=0;
				best=Integer.MAX_VALUE;
				bestBlocks=0;
				time = 0;
				
					for (int block = 1; block < 100; block++){
						
						for (int l = 0;l<10;l++){ //10 iterations to take average over, for precision.

							
							int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

							Arrays.fill(levels, j); //lines that create the objects we are to work on. Should not be timed, ofc.
							CG = new CombiGrid(levels);
							Arrays.fill(CG.grid, 1);					
							nanoStart = System.currentTimeMillis();
							CG.hierarchizeUnoptimizedParallelStream(block);
							//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
							nanoStop = System.currentTimeMillis();

							currentTime = nanoStop-nanoStart; 		
							time += currentTime;
							
						}
						if (time/10 < best)  {
							bestBlocks=block;
							best=currentTime;
						}
							
						if (PrintAllResults) System.out.println(i+"\t"+j+"\t"+block+"\t"+currentTime);
					}
					if (!PrintAllResults) System.out.println(i+"\t"+j+"\t"+bestBlocks+"\t"+best); // will print info for Thread
					System.out.println("");
				}
			}
		}
	
	static void testFullTasks(){ 
		System.out.println("dim\tlev\tthreads/tasks\ttime in ms\ttype");
		
		for (int i = 2; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<7;j++){// this loop will relate to number of levels.
				bestThreads=0;
				double best = Integer.MAX_VALUE;
				for (k = 1; k<10;k++){//this loop iterates over number of threads.
					//System.out.println("Test with " + i + " dimensions, all with level " + j + ", running in "+k+" threads.");
					time = 0;
					for (int l = 0;l<10;l++){

						int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

						Arrays.fill(levels, j); //lines that create the objects we are to work on. Should not be timed, ofc.
						CG = new CombiGrid(levels);
						Arrays.fill(CG.grid, 1);					
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedTasks(k);
						nanoStop = System.currentTimeMillis();

						currentTime = nanoStop-nanoStart; 		
						time += currentTime;
					}
					System.out.println(i+"\t"+j+"\t"+k+"\t"+time/10+"\tTasks"); // will print info for Thread
				}
			}
		}		
	}
	
	
	static void testAllSimple(int[] levels, GridFunctions GF){
		//System.out.println("Running test on grid with "+ levels.length " dimensions, and ");
		//for (each  
		//test threads basic
		System.out.println("Threads");
		System.out.println("Threads\tmilliseconds");
		for (k = 1; k<11;k++){//this loop iterates over number of threads.
				CG = new CombiGrid(levels);
				CG.setValues(GF);
				
				nanoStart = System.currentTimeMillis();
				CG.hierarchizeUnoptimizedThreads(k);
				nanoStop = System.currentTimeMillis();

				currentTime = nanoStop-nanoStart; 		
				System.out.println(k+"\t"+currentTime);
			}
		
				System.out.println("ThreadsOnce");
				System.out.println("threads\tmilliseconds");
				for (k = 1; k<11;k++){//this loop iterates over number of threads.
						CG = new CombiGrid(levels);
						CG.setValues(GF);		
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedThreadsOnce(k);
						nanoStop = System.currentTimeMillis();

						currentTime = nanoStop-nanoStart; 		
						System .out.println(k+"\t"+currentTime);
					}
			
				
				System.out.println("P-stream");
				System.out.println("blocks\tmilliseconds");
				for (k = 1; k<100;k++){//this loop iterates over number of blocks.

					CG = new CombiGrid(levels);
					CG.setValues(GF);				
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedParallelStream(k);;
						nanoStop = System.currentTimeMillis();
						currentTime = nanoStop-nanoStart; 		
						System.out.println(k+"\t"+currentTime);
					}

				//test tasks basic
				System.out.println("Tasks");
				System.out.println("Tasks\tmilliseconds");
				for (k = 1; k<100;k++){//this loop iterates over number of threads.
						CG = new CombiGrid(levels);
						CG.setValues(GF);		
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedTasks(k);
						nanoStop = System.currentTimeMillis();
						currentTime = nanoStop-nanoStart; 		
						System.out.println(k+"\t"+currentTime);
					}
				System.out.println("Successfully finished.");
			}

	static void testSideways(){ //This method makes a big run, comparing the best result of each method, for each successive dim/lev.
		System.out.println("dim\tlev\tthreads\ttime in ms\ttype");
		int bestTaskAmount=0;
		float bestTimeTask=0;
		float currentTimeONCE = 0;
		float bestThreadedONCE =0;
		int bestThreadAmountONCE=0;

		for (int i = 3; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<7;j++){// this loop will relate to number of levels.

				time = 0;
				timeStream=0;
				currentTime=0;
				timeTask=0;
				bestThreadAmount=0;

				for (int l = 0;l<10;l++){
					time = 0;
					timeStream=0;
					currentTime=0;
					timeTask=0;
					bestThreadAmount=0;
					bestTaskAmount=0;
					bestTimeTask=0;
					currentTimeONCE = 0;
					bestThreadedONCE =0;
					bestThreadAmountONCE=0;
					
					int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

					//lines that create the objects we are to work on. Should not be timed, ofc.
					Arrays.fill(levels, j);
					CG = new CombiGrid(levels);
					Arrays.fill(CG.grid, 1);		
					bestThreaded=Integer.MAX_VALUE; //for holding the best value, so we only compare the best results.

					for (k = 1; k<cores*3;k++){//this loop iterates over number of threads.
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedThreads(k);
						nanoStop = System.currentTimeMillis();
						currentTime = nanoStop-nanoStart;
						if (currentTime < bestThreaded){
							bestThreaded = currentTime;
							bestThreadAmount=k;
						}
					}
					
					//THREADS ONCE
					Arrays.fill(levels, j);
					CG = new CombiGrid(levels);
					Arrays.fill(CG.grid, 1);		
					bestThreadedONCE=Integer.MAX_VALUE; //for holding the best value, so we only compare the best results.

					for (int o = 1; o<cores*3;o++){//this loop iterates over number of threads.
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedThreadsOnce(o);
						//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
						nanoStop = System.currentTimeMillis();
						currentTimeONCE = nanoStop-nanoStart;
						if (currentTimeONCE < bestThreadedONCE){
							bestThreadedONCE = currentTimeONCE;
							bestThreadAmountONCE=o;
						}
					}
					
					//Runs P-STREAM
					for (int n=0; n<10; n++){
						CG = new CombiGrid(levels); //get a new grid for comparison.
						Arrays.fill(CG.grid, 1);

						nanoStartStream = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedParallelStream();
						//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
						nanoStopStream = System.currentTimeMillis();

						currentTimeStream = nanoStopStream-nanoStartStream; 		
						timeStream += currentTimeStream;
					}
					
					//Runs P-STREAM finding and using the optimal amount of blocks.
					bestBlocks=0;
					bestBlockTime=Integer.MAX_VALUE;
					for (int blocks = 1; blocks<100;blocks++) {
						timeStreamBlock=0;
						for (int n=0; n<10; n++){

							CG = new CombiGrid(levels); //get a new grid for comparison.
							Arrays.fill(CG.grid, 1);

							nanoStartStream = System.currentTimeMillis();
							CG.hierarchizeUnoptimizedParallelStream(blocks);
							//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
							nanoStopStream = System.currentTimeMillis();

							currentTimeStream = nanoStopStream-nanoStartStream; 		
							timeStreamBlock += currentTimeStream;
						}
						if (timeStreamBlock/10<bestBlockTime){
							bestBlocks=blocks;
							bestBlockTime=timeStreamBlock/10;
						}
					}
					
					currentBestTaskTime=Integer.MAX_VALUE;
					
					//Runs tasks
					for (int m=1; m<10; m++){
						CG = new CombiGrid(levels); //get a new grid for comparison.
						Arrays.fill(CG.grid, 1);

						nanoStartTask = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedTasks(m);
						//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
						nanoStopTask = System.currentTimeMillis();

						currentTimeTask = nanoStopTask-nanoStartTask; 		
						timeTask += currentTimeTask;	//Runs tasks
						
						if (currentTimeTask<currentBestTaskTime){
							currentBestTaskTime=currentTimeTask;
							bestTaskAmount=m;
						}
					}
				}
				System.out.println(i+"\t"+j+"\t"+bestThreadAmount+"\t"+bestThreaded+"\tThreads"); // will print info for Thread
				System.out.println(i+"\t"+j+"\t"+bestThreadAmountONCE+"\t"+bestThreadedONCE+"\tThreads ONCE"); // will print info for Thread
				System.out.println(i+"\t"+j+"\tnull\t"+timeStream/10+"\tP-Stream");
				System.out.println(i+"\t"+j+"\t"+bestBlocks+"\t"+bestBlockTime+"\tP-Stream Blocks");
				System.out.println(i+"\t"+j+"\t"+bestTaskAmount+"\t"+currentBestTaskTime+"\tTasks");
				System.out.println("");
			}
		}
	}		
}
