package testing;
import grid.*;
import gridFunctions.*;

import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.util.List;

public class ScalingTest {

	public static void main(String[] args) {
		//int[] levels = {8, 6, 5, 3, 2};
		printInfo();
		int[] levels = {5, 5, 5, 5, 5};
		CombiGrid cg = new CombiGrid(levels);
		CombiGridAligned cga = new CombiGridAligned(levels, 32);
		//test(cg, 10);
		//testOptimized(cg, 10);
		//testOptimized(cga, 10);
		//testOptimizedNoUnroll(cg, 10);
		//testThreads(cg, 10);
		//testOptimizedThreads(cg, 10);
		//testOptimizedThreads(cga, 10);
		//testOptimizedThreadsNoUnroll(cg, 10);
		testRecursive(cga, 10);
		testRecursiveThread(cga, 10);
		
//		int[] levels2 = {10, 10};
		//CombiGrid cg = new CombiGrid(levels);
//		RecursiveThreadedFindOptimalValues(cga, 10, 3);
		
		/*for (int i=20;i<31;i++) {
			System.out.println("Now running on a grid of size "+ i);
			CombiGridAligned cga = createAlignedGrid(i);
			//CombiGridAligned cga = new CombiGridAligned(levels2, 32);
			RecursiveThreadedFindOptimalValues(cga, 10, 32);
			cga=null;
			System.gc();
		}*/
		//CombiGridAligned cga = new CombiGridAligned(levels2, 32);
		testRecursiveLinear(cga, 10);
//		cga=null;testRecursiveThreaded
//		System.gc();
//		
//		cga = new CombiGridAligned(levels2, 32);
		testRecursiveThreaded(cga, 10, 32);
//		cga=null;
//		System.gc();
//		
//		cga = new CombiGridAligned(levels2, 32);
//		testOptimizedThreads(cga, 4);
//		cga=null;
//		System.gc();
		
		
		//CombiGridAligned cga = new CombiGridAligned(levels, 32);
		//testAllUnoptimized();
		//testAllOptimized();
		
		/*for (int i = 25; i <32 ; i++){
			CombiGridAligned cga = createAlignedGrid(i);
			System.out.println("-------------------------------------");
			System.out.println("Now working on a grid of size " + i);
			System.out.println("Level array: ");
			int j =0;
			while (j<cga.levels.length){
				System.out.print(cga.levels[j] + " ");
				j++;
			}
		System.out.println("-------------------------------------");
			
			
			testOptimized(cga, 5);
			testOptimizedThreads(cga, 5);
			testOptimizedThreadsOnce(cga, 5);
			testOptimizedTasks(cga, 5);
			
			cga=null;
			CombiGrid cg = createGrid(i);

			test(cg, 5);
			testThreads(cg, 5);
			testThreadsOnce(cg, 5);
			testTasks(cg, 5);

			cg=null;
			//testParallelStream(cg, 10);
			//testParallelStreamChunks(cg, 10);
			//testOptimizedParallelStreamChunks(cga, 10);
		}
		*/
		//testAllOptimized();
	}
	
	private static void testAllNonAligned() {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;
		double avgTimeThreads;
		double avgTimeThreadsOnce;
		double avgTimeTasks;
		//double avgTimeParallelStream;
		//double avgTimeParallelStreamBlocks;
		double avgTimeOptimized;
		double avgTimeOptimizedThreads;
		double avgTimeNoUnroll;
		double avgTimeOptimizedThreadsNoUnroll;
		
		System.out.println("Testing all Unoptimized");
		//System.out.println("Gridsize" + '\t' + "Sequentiel" + '\t' + "Threads" + '\t' + "ThreadsOnce" + '\t' + "Tasks" + '\t' + "ParallelStream" + '\t' + "ParallelStreamChunks");
		System.out.println("Gridsize" + '\t' + "Sequentiel" + '\t' + "Threads" + '\t' + "ThreadsOnce" + '\t' + "Tasks" + '\t' + "Optimized" + '\t' + "OptimizedThreads" + '\t' + "OptimizedNoUnroll" + '\t' + "avgTimeOptimizedThreadsNoUnroll");
		for(int size = 19; size <= 30; size++) {
			CombiGrid cg = createGrid(size);
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimized();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedThreads(16);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeThreads = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedThreadsOnce(16);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeThreadsOnce = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedTasks(16);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeTasks = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimized();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeOptimized = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreads(16);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeOptimizedThreads = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedNoUnroll();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeNoUnroll = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreadsNoUnroll(16);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeOptimizedThreadsNoUnroll = totalTime / 10;
			
//			totalTime = 0;
//			for(int i = 0; i < 10; i++) {
//				cg.setValues(GridFunctions.ALLONES);
//				start = System.currentTimeMillis();
//				cg.hierarchizeUnoptimizedParallelStream();
//				end = System.currentTimeMillis();
//				time = end - start;
//				totalTime += time;
//			}
//			avgTimeParallelStream = totalTime / 10;
			
//			totalTime = 0;
//			for(int i = 0; i < 10; i++) {
//				cg.setValues(GridFunctions.ALLONES);
//				start = System.currentTimeMillis();
//				cg.hierarchizeUnoptimizedParallelStream(4);
//				end = System.currentTimeMillis();
//				time = end - start;
//				totalTime += time;
//			}
//			avgTimeParallelStreamBlocks = totalTime / 10;
//			System.out.println("" + size + '\t' + avgTime + '\t' + avgTimeThreads + '\t' + avgTimeThreadsOnce + '\t' + avgTimeTasks + '\t' + avgTimeParallelStream + '\t' + avgTimeParallelStreamBlocks);
			System.out.println("" + size + '\t' + avgTime + '\t' + avgTimeThreads + '\t' + avgTimeThreadsOnce + '\t' + avgTimeTasks + '\t' + avgTimeOptimized + '\t' + avgTimeOptimizedThreads + '\t' + avgTimeNoUnroll + '\t' + avgTimeOptimizedThreadsNoUnroll);
			cg = null;
			System.gc();
		}
	}
	
	private static void testAllOptimized() {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;
		double avgTimeThreads;
		double avgTimeThreadsOnce;
		double avgTimeTasks;
		//double avgTimeParallelStreamBlocks;
		double avgTimeRecursive;
		double avgTimeRecursiveThreads;
		
		System.out.println("Testing all Optimized");
		//System.out.println("Gridsize" + '\t' + "Sequentiel" + '\t' + "Threads" + '\t' + "ThreadsOnce" + '\t' + "Tasks" + '\t' + "ParallelStreamChunks");
		System.out.println("Gridsize" + '\t' + "Sequentiel" + '\t' + "Threads" + '\t' + "ThreadsOnce" + '\t' + "Tasks"+'\t'+"Recursive Linear"+'\t'+"Recursive Threaded");
		
		for(int size = 19; size <= 30; size++) {
			CombiGridAligned cg = createAlignedGrid(size);
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimized();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreads(16);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeThreads = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreadsOnce(16);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeThreadsOnce = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedTasks(16);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeTasks = totalTime / 10;
			
//			totalTime = 0;
//			for(int i = 0; i < 10; i++) {
//				cg.setValues(GridFunctions.ALLONES);
//				start = System.currentTimeMillis();
//				cg.hierarchizeOptimizedParallelStream(4, 4);
//				end = System.currentTimeMillis();
//				time = end - start;
//				totalTime += time;
//			}
//			avgTimeParallelStreamBlocks = totalTime / 10;
//			System.out.println("" + size + '\t' + avgTime + '\t' + avgTimeThreads + '\t' + avgTimeThreadsOnce + '\t' + avgTimeTasks + '\t' + avgTimeParallelStreamBlocks);
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeRecursive();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeRecursive = totalTime / 10;

			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeRecursiveThreads(2);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeRecursiveThreads = totalTime / 10;
			
			System.out.println("" + size + '\t' + avgTime + '\t' + avgTimeThreads + '\t' + avgTimeThreadsOnce + '\t' + avgTimeTasks+'\t'+avgTimeRecursive+'\t'+avgTimeRecursiveThreads);
		}
	}

	private static void test(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimized");
		System.out.println("Run" + '\t' + "Time in ms");
		
		for(int run = 1; run <= 1; run++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimized();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + run + '\t' + avgTime);
		}
	}
	
	private static void testRecursiveLinear(CombiGridAligned cg, int repititions){
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;
		
		System.out.println("HierarchizeRecursiveLinear");
		System.out.println("Run" + '\t' + "Time in ms");
		
		for(int run = 1; run <= 1; run++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeRecursive();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + run + '\t' + avgTime);
		}
	}
	
	private static void testRecursiveThreaded(CombiGridAligned cg, int repititions, int NumberOfThreads){
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;
		
		System.out.println("HierarchizeRecursiveThreaded");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//The following lines will set the recursion values.
		int LevelMax =0;
		for (int l : cg.levels){
			if (l>LevelMax) LevelMax = l; //For varying the recMin and recMax.
		}
		if (LevelMax>3) cg.recMaxSpawn=LevelMax-2; // See page 20 in SGA2014.
		else cg.recMaxSpawn = LevelMax;
		
		if (LevelMax>3) cg.recMinSpawn = cg.recMaxSpawn - 2; //Also set based on info in paper for now.
		else cg.recMinSpawn = 2; //At least this should not affect the performance poorly.
		
		
		
		cg.recTile=3; //More here doesn't seem to make it faster currently.
		
		System.out.println("Recursion values:");
		System.out.println("Max: " + cg.recMaxSpawn+", Min: " + cg.recMinSpawn);
		
		for(int run = 1; run <= 1; run++) {
			totalTime = 0;
			
			for (int j = 1; j <= NumberOfThreads; j++){
				totalTime = 0;
				for(int i = 0; i < repititions; i++) {
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeRecursiveThreads(j);
					end = System.currentTimeMillis();
					time = end - start;
					totalTime += time;
					//System.out.println(time);
				}
			avgTime = totalTime / repititions;
			System.out.println("" + j + '\t' + avgTime);
		}
	}
	}
	
	private static void RecursiveThreadedFindOptimalValues(CombiGridAligned cg, int repititions, int MaxNumberOfThreads){
		//This method will vary all the input variables of the recursive method, including number of threads up to the input, but also when to stop the recursion.
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;
		
		System.out.println("RecursiveThreadedFindOptimalValues");
		System.out.println("recMaxSpawn" +"\t"+ "RecMinSpawn"+"\t" + "recTile" + '\t' +"Threads"+"\t"+ "Time in ms");
		int LevelMax = 0;
		for (int l : cg.levels){
			if (l>LevelMax) LevelMax = l; //For varying the recMin and recMax.
		}

		double currentBestTime = Integer.MAX_VALUE;
		int bestMin=0;
		int bestMax=0;
		int bestTile=0;
		int bestThreads=0;
		totalTime = 0;
		int NumberOfThreads=1;
		cg.recTile=0;
		cg.recMaxSpawn=1;
		cg.recMinSpawn=0;

		for (cg.recMaxSpawn++; cg.recMaxSpawn<=LevelMax;cg.recMaxSpawn++){
			cg.recMinSpawn=0;

			for (cg.recMinSpawn++;cg.recMinSpawn<cg.recMaxSpawn;cg.recMinSpawn++){
				cg.recTile=0;

				for (cg.recTile++;cg.recTile<=cg.levels[0];cg.recTile++) { //recTile set to max lvl in dim 0, according to paper. 
					NumberOfThreads=0;

					for (NumberOfThreads++;NumberOfThreads<=MaxNumberOfThreads;NumberOfThreads++){	
						totalTime=0;
						
						for(int i = 0; i < repititions; i++) {

							cg.setValues(GridFunctions.ALLONES);
							start = System.currentTimeMillis();
							cg.hierarchizeRecursiveThreads(NumberOfThreads);
							end = System.currentTimeMillis();
							time = end - start;
							totalTime += time;
							System.gc();
							//System.out.println(time); //Un-comment, to get each time output.
						}
						avgTime = totalTime / repititions;
						if (avgTime < currentBestTime) {
							currentBestTime=avgTime;
							bestMin = cg.recMinSpawn;
							bestMax = cg.recMaxSpawn;
							bestThreads=NumberOfThreads;
							bestTile=cg.recTile;
						}
						System.out.println(cg.recMaxSpawn+"\t" +cg.recMinSpawn+"\t"+ cg.recTile+ '\t'+NumberOfThreads+"\t" + avgTime);
					}
				}
			}
		}


		System.out.println("Best values:");
		System.out.println("Time:\t" + currentBestTime+" ms");
		System.out.println("RecMinLevel:\t"+ bestMin);
		System.out.println("RecMaxLevel:\t"+bestMax);
		System.out.println("Threads:\t"+bestThreads);
		System.out.println("RecTile:\t"+bestTile);
	}
	
	private static void testOptimized(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimized (non-aligned)");
		System.out.println("Run" + '\t' + "Time in ms");
		
		for(int run = 1; run <= 1; run++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimized();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + run + '\t' + avgTime);
		}
	}
	
	private static void testOptimizedNoUnroll(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimizedNoUnroll (non-aligned)");
		System.out.println("Run" + '\t' + "Time in ms");
		
		for(int run = 1; run <= 1; run++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedNoUnroll();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + run + '\t' + avgTime);
		}
	}
	
	private static void testOptimized(CombiGridAligned cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimized (aligned)");
		System.out.println("Run" + '\t' + "Time in ms");
		
		for(int run = 1; run <= 1; run++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimized();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + run + '\t' + avgTime);
		}
	}

	private static void testThreads(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedThreads");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-32 threads
		for(int threads = 1; threads <= 32; threads++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedThreads(threads);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}
	
	private static void testOptimizedThreads(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimizedThreads (non-aligned)");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-32 threads
		for(int threads = 1; threads <= 32; threads++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreads(threads);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}
	
	private static void testOptimizedThreadsNoUnroll(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimizedThreadsNoUnroll (non-aligned)");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-32 threads
		for(int threads = 1; threads <= 32; threads++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreadsNoUnroll(threads);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}
	
	private static void testOptimizedThreads(CombiGridAligned cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimizedThreads");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-32 threads
		for(int threads = 1; threads <= 32; threads++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreads(threads);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}

	private static void testThreadsOnce(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedThreadsOnce");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-32 threads
		for(int threads = 1; threads <= 32; threads++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedThreadsOnce(threads);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}
	
	private static void testOptimizedThreadsOnce(CombiGridAligned cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimizedThreadsOnce");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-32 threads
		for(int threads = 1; threads <= 32; threads++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreadsOnce(threads);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}
	
	private static void testTasks(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedTasks");
		System.out.println("Tasks" + '\t' + "Time in ms");
		
		//Run for 1-32 tasks
		for(int tasks = 1; tasks <= 32; tasks++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedTasks(tasks);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + tasks + '\t' + avgTime);
		}
	}
	
	private static void testOptimizedTasks(CombiGridAligned cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimizedTasks");
		System.out.println("Tasks" + '\t' + "Time in ms");
		
		//Run for 1-32 tasks
		for(int tasks = 1; tasks <= 32; tasks++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedTasks(tasks);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + tasks + '\t' + avgTime);
		}
	}
	
	private static void testRecursive(CombiGridAligned cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeRecursive");
		System.out.println("Run" + '\t' + "Time in ms");
		
		//Run for 1-32 threads
		for(int threads = 1; threads <= 1; threads++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeRecursive();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}
	
	private static void testRecursiveThread(CombiGridAligned cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeRecursiveThread");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-32 threads
		for(int threads = 1; threads <= 8; threads++) {
			totalTime = 0;
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeRecursiveThread(threads);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}
	
//	private static void testParallelStream(CombiGrid cg, int repititions) {
//		long start;
//		long end;
//		double time;
//		double totalTime;
//		double avgTime;
//
//		System.out.println("HierarchizeUnoptimizedParallelStream");
//		System.out.println("Run" + '\t' + "Time in ms");
//		
//		//Run 10 tests
//		for(int parallelStream = 1; parallelStream <= 16; parallelStream++) {
//			totalTime = 0;
//			//Run 10 times to average result
//			for(int i = 0; i < repititions; i++) {
//				cg.setValues(GridFunctions.ALLONE);
//				start = System.currentTimeMillis();
//				cg.hierarchizeUnoptimizedParallelStream();
//				end = System.currentTimeMillis();
//				time = end - start;
//				totalTime += time;
//			}
//			avgTime = totalTime / repititions;
//			System.out.println("" + parallelStream + '\t' + avgTime);
//		}
//	}
	
//	private static void testParallelStreamChunks(CombiGrid cg, int repititions) {
//		long start;
//		long end;
//		double time;
//		double totalTime;
//		double avgTime;
//
//		System.out.println("HierarchizeUnoptimizedParallelStreamChunks");
//		System.out.println("Chunks" + '\t' + "Time in ms");
//		
//		//Run for 1-100 chunks
//		for(int chunks = 1; chunks <= 16; chunks++) {
//			totalTime = 0;
//			//Run 10 times to average result
//			for(int i = 0; i < repititions; i++) {
//				cg.setValues(GridFunctions.ALLONES);
//				start = System.currentTimeMillis();
//				cg.hierarchizeUnoptimizedParallelStream(chunks);
//				end = System.currentTimeMillis();
//				time = end - start;
//				totalTime += time;
//			}
//			avgTime = totalTime / repititions;
//			System.out.println("" + chunks + '\t' + avgTime);
//		}
//	}
	
//	private static void testOptimizedParallelStreamChunks(CombiGridAligned cg, int repititions) {
//		long start;
//		long end;
//		double time;
//		double totalTime;
//		double avgTime;
//
//		System.out.println("HierarchizeOptimizedParallelStreamChunks");
//		System.out.println("Chunks" + '\t' + "Time in ms");
//		
//		//Run for 1-100 chunks
//		for(int chunks = 1; chunks <= 16; chunks++) {
//			totalTime = 0;
//			//Run 10 times to average result
//			for(int i = 0; i < repititions; i++) {
//				cg.setValues(GridFunctions.ALLONES);
//				start = System.currentTimeMillis();
//				cg.hierarchizeOptimizedParallelStream(16, chunks);
//				end = System.currentTimeMillis();
//				time = end - start;
//				totalTime += time;
//			}
//			avgTime = totalTime / repititions;
//			System.out.println("" + chunks + '\t' + avgTime);
//		}
//	}
	
	private static CombiGrid createGrid(int size) {		
		if(size < 4) {
			int[] levels = {size};
			return new CombiGrid(levels);
		}
		
		else if(size == 4) {
			int[] levels = {2, 2};
			return new CombiGrid(levels);
		}
		
		else if(size == 5) {
			int[] levels = {3, 2};
			return new CombiGrid(levels);
		}
		
		else if(size == 6) {
			int[] levels = {2, 2, 2};
			return new CombiGrid(levels);
		}
		
		else if(size == 7) {
			int[] levels = {3, 2, 2};
			return new CombiGrid(levels);
		}
		
		else if(size == 8) {
			int[] levels = {2, 2, 2, 2};
			return new CombiGrid(levels);
		}
		
		else if(size == 9) {
			int[] levels = {3, 2, 2, 2};
			return new CombiGrid(levels);
		}
		
		else {
			int dim = size / 5;
			int rem = size % 5;
			int[] levels = {dim + rem, dim, dim, dim, dim};
			return new CombiGrid(levels);
		}
	}
	
	private static CombiGridAligned createAlignedGrid(int size) {		
		if(size < 4) {
			int[] levels = {size};
			return new CombiGridAligned(levels, 32);
		}
		
		else if(size == 4) {
			int[] levels = {2, 2};
			return new CombiGridAligned(levels, 32);
		}
		
		else if(size == 5) {
			int[] levels = {3, 2};
			return new CombiGridAligned(levels, 32);
		}
		
		else if(size == 6) {
			int[] levels = {2, 2, 2};
			return new CombiGridAligned(levels, 32);
		}
		
		else if(size == 7) {
			int[] levels = {3, 2, 2};
			return new CombiGridAligned(levels, 32);
		}
		
		else if(size == 8) {
			int[] levels = {2, 2, 2, 2};
			return new CombiGridAligned(levels, 32);
		}
		
		else if(size == 9) {
			int[] levels = {3, 2, 2, 2};
			return new CombiGridAligned(levels, 32);
		}
		
		else {
			int dim = size / 5;
			int rem = size % 5;
			int[] levels = {dim + rem, dim, dim, dim, dim};
			for(int i : levels)
				System.out.print(i + " ");
			return new CombiGridAligned(levels, 32);
		}
	}
	
	private static void printInfo() {
		RuntimeMXBean runtimeMxBean = ManagementFactory.getRuntimeMXBean();
		List<String> arguments = runtimeMxBean.getInputArguments();
		System.out.printf("# OS: %s; %s; %s%n",
		System.getProperty("os.name"),
		System.getProperty("os.version"),
		System.getProperty("os.arch"));
		System.out.printf("# JVM: %s; %s%n",
		System.getProperty("java.vendor"),
		System.getProperty("java.version"));
		// This line works only on MS Windows:
		System.out.printf("# CPU: %s%n", System.getenv("PROCESSOR_IDENTIFIER"));
		System.out.println("# Available processors: " + Runtime.getRuntime().availableProcessors());
		System.out.println("# Maximum memory: " + Runtime.getRuntime().maxMemory() / (1024 * 1024) + " MB");
		java.util.Date now = new java.util.Date();
		System.out.printf("# Date: %s%n",
		new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ssZ").format(now));
		System.out.printf("# VM Arguments: ");
		for(String s : arguments)
			System.out.printf(s + " ");
		System.out.println();
	}
}
