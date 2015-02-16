import java.util.Arrays;

public class testSuite {

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
	

	//This class can run various tests on the CombiGrid.java code.
	//supposed to test amount of cores and threads, threading methods and time. 
	public static void main(String[] args) {
		cores = Runtime.getRuntime().availableProcessors();
		mem = Runtime.getRuntime().maxMemory()/(1024*1024);
		System.out.println("Running tests on "+cores+" cores using up to "+mem+" MB of memory.");
		
		//Which tests to run. The following three will increase grid size gradually, and test the average runtime over 10 iterations.
		//testFullThreads(); //Note that these are designed to crash, to you can only run one.
		//testFullTasks();
		//testFullPStream();
		
		//This test will compare the best results for each method, and only print those.
		testSideways();
		
	}
	static void testFullThreads(){ //TODO make this method compare the different threading methods.
		System.out.println("dim\tlev\tthreads/tasks\ttime\ttype");
		
		for (int i = 2; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<10;j++){// this loop will relate to number of levels.
				bestThreads=0;
				double best = Integer.MAX_VALUE;
				for (k = 1; k<cores*3;k++){//this loop iterates over number of threads.
					//System.out.println("Test with " + i + " dimensions, all with level " + j + ", running in "+k+" threads.");
					time = 0;
					for (int l = 0;l<10;l++){

						int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

						Arrays.fill(levels, j); //lines that create the objects we are to work on. Should not be timed, ofc.
						CG = new CombiGrid(i, levels);
						Arrays.fill(CG.grid, 1);					
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedThreads(k);
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
	
	static void testFullPStream(){ //TODO make this method compare the different threading methods.
		System.out.println("dim\tlev\tthreads/tasks\ttime in ms\ttype");
		
		for (int i = 2; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<10;j++){// this loop will relate to number of levels.
				bestThreads=0;
				double best = Integer.MAX_VALUE;
					time = 0;
					for (int l = 0;l<10;l++){

						int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

						Arrays.fill(levels, j); //lines that create the objects we are to work on. Should not be timed, ofc.
						CG = new CombiGrid(i, levels);
						Arrays.fill(CG.grid, 1);					
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedParallelStream();
						//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
						nanoStop = System.currentTimeMillis();

						currentTime = nanoStop-nanoStart; 		
						time += currentTime;
					}
					System.out.println(i+"\t"+j+"\t"+"not known"+"\t"+time/10+"\tThreads"); // will print info for Thread
				}
			}
		}		
	
	
	static void testFullTasks(){ //TODO make this method compare the different threading methods.
		System.out.println("dim\tlev\tthreads/tasks\ttime in ms\ttype");
		
		for (int i = 2; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<10;j++){// this loop will relate to number of levels.
				bestThreads=0;
				double best = Integer.MAX_VALUE;
				for (k = 1; k<10;k++){//this loop iterates over number of threads.
					//System.out.println("Test with " + i + " dimensions, all with level " + j + ", running in "+k+" threads.");
					time = 0;
					for (int l = 0;l<10;l++){

						int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

						Arrays.fill(levels, j); //lines that create the objects we are to work on. Should not be timed, ofc.
						CG = new CombiGrid(i, levels);
						Arrays.fill(CG.grid, 1);					
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedTasks(k);
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
	
	

	static void testSideways(){ //TODO make this method compare the different threading methods.
		System.out.println("dim\tlev\tthreads\ttime in ms\ttype");

		for (int i = 3; i<7; i++){// this loop will relate to the number of dimensions, as this scales strongest		
			for (int j = 2; j<7;j++){// this loop will relate to number of levels.
				time = 0;
				timeStream=0;
				bestThreadAmount=0;
				int bestTaskAmount=0;
				float bestTimeTask=0;
				for (int l = 0;l<10;l++){

					int[] levels = new int[i]; //array of levels, must be as long as the amount of dimensions.

					//lines that create the objects we are to work on. Should not be timed, ofc.
					Arrays.fill(levels, j);
					CG = new CombiGrid(i, levels);
					Arrays.fill(CG.grid, 1);		
					bestThreaded=Integer.MAX_VALUE; //for holding the best value, so we only compare the best results.

					for (k = 1; k<cores*3;k++){//this loop iterates over number of threads.
						nanoStart = System.currentTimeMillis();
						CG.hierarchizeUnoptimizedThreads(k);
						//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
						nanoStop = System.currentTimeMillis();
						currentTime = nanoStop-nanoStart;
						if (currentTime < bestThreaded){
							bestThreaded = currentTime;
							bestThreadAmount=k;
						}
					}
					

					//Runs P-STREAM
					CG = new CombiGrid(i, levels); //get a new grid for comparison.
					Arrays.fill(CG.grid, 1);

					nanoStartStream = System.currentTimeMillis();
					CG.hierarchizeUnoptimizedParallelStream();
					//	CG.printValues(); //print only to check that everything is being calculated correctly. Takes a lot of extra time.
					nanoStopStream = System.currentTimeMillis();

					currentTimeStream = nanoStopStream-nanoStartStream; 		
					timeStream += currentTimeStream;
					
					float currentBestTaskTime=Integer.MAX_VALUE;
					//Runs tasks
					for (int m =1; m<10; m++){
						CG = new CombiGrid(i, levels); //get a new grid for comparison.
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
						
							bestTimeTask=currentBestTaskTime;
							timeTask += currentTimeTask;	//Runs tasks
					}
					
				}
				System.out.println(i+"\t"+j+"\t"+bestThreadAmount+"\t"+bestThreaded+"\tThreads"); // will print info for Thread
				System.out.println(i+"\t"+j+"\tnull\t"+timeStream/10+"\tP-Stream");
				//System.out.println("D" + i + "-L" + j + "\t" + (timeStream/10) + " millis to complete in average for P-Stream");
				System.out.println(i+"\t"+j+"\t"+bestTaskAmount+"\t"+bestTimeTask+"\tTasks");
				System.out.println("");
			}
		}
	}		
}