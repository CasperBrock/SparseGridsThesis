package testing;
import grid.*;
import gridFunctions.*;


public class ScalingTest {

	public static void main(String[] args) {
		//int[] levels = {10, 6, 4, 3, 2};
		int[] levels = {5, 5, 5, 5, 5};
		CombiGrid cg = new CombiGrid(5, levels);
		testAll();
		//test(cg, 10);
		//testThreads(cg, 10);
		//testThreadsOnce(cg, 10);
		//testTasks(cg, 10);
		//testParallelStream(cg, 10);
		//testParallelStreamBlocks(cg, 10);
	}
	
	private static void testAll() {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;
		double avgTimeThreads;
		double avgTimeThreadsOnce;
		double avgTimeTasks;
		double avgTimeParallelStream;
		double avgTimeParallelStreamBlocks;
		
		System.out.println("Testing all");
		System.out.println("Gridsize" + '\t' + "Sequentiel" + '\t' + "Threads" + '\t' + "ThreadsOnce" + '\t' + "Tasks" + '\t' + "ParallelStream" + '\t' + "ParallelStreamBlock");
		
		for(int size = 19; size <= 26; size++) {
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
				cg.hierarchizeUnoptimizedThreads(3);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeThreads = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedThreadsOnce(4);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeThreadsOnce = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedTasks(22);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeTasks = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedParallelStream();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeParallelStream = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedParallelStream(80);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeParallelStreamBlocks = totalTime / 10;
			
			System.out.println("" + size + '\t' + avgTime + '\t' + avgTimeThreads + '\t' + avgTimeThreadsOnce + '\t' + avgTimeTasks + '\t' + avgTimeParallelStream + '\t' + avgTimeParallelStreamBlocks);
		}
	}

	private static void test(CombiGrid cg, int runs) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimized");
		System.out.println("Run" + '\t' + "Time in ms");
		
		//Run 10 times
		for(int run = 1; run <= 10; run++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < runs; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimized();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / runs;
			System.out.println("" + run + '\t' + avgTime);
		}
	}

	private static void testThreads(CombiGrid cg, int runs) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedThreads");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-10 threads
		for(int threads = 1; threads <= 10; threads++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < runs; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedThreads(threads);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / runs;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}

	private static void testThreadsOnce(CombiGrid cg, int runs) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedThreadsOnce");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-10 threads
		for(int threads = 1; threads <= 10; threads++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < runs; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedThreadsOnce(threads);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / runs;
			System.out.println("" + threads + '\t' + avgTime);
		}
	}
	
	private static void testTasks(CombiGrid cg, int runs) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedTasks");
		System.out.println("Tasks" + '\t' + "Time in ms");
		
		//Run for 1-100 tasks
		for(int tasks = 1; tasks <= 100; tasks++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < runs; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedTasks(tasks);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / runs;
			System.out.println("" + tasks + '\t' + avgTime);
		}
	}
	
	private static void testParallelStream(CombiGrid cg, int runs) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedParallelStream");
		System.out.println("Run" + '\t' + "Time in ms");
		
		//Run 10 tests
		for(int parallelStream = 1; parallelStream <= 10; parallelStream++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < runs; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedParallelStream();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / runs;
			System.out.println("" + parallelStream + '\t' + avgTime);
		}
	}
	
	private static void testParallelStreamBlocks(CombiGrid cg, int runs) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedParallelStreamBlocks");
		System.out.println("Blocks" + '\t' + "Time in ms");
		
		//Run for 1-100 blocks
		for(int blocks = 1; blocks <= 100; blocks++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < runs; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedParallelStream(blocks);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / runs;
			System.out.println("" + blocks + '\t' + avgTime);
		}
	}
	
	private static CombiGrid createGrid(int size) {		
		if(size < 4) {
			int[] levels = {size};
			return new CombiGrid(1, levels);
		}
		
		else if(size == 4) {
			int[] levels = {2, 2};
			return new CombiGrid(2, levels);
		}
		
		else if(size == 5) {
			int[] levels = {3, 2};
			return new CombiGrid(2, levels);
		}
		
		else if(size == 6) {
			int[] levels = {2, 2, 2};
			return new CombiGrid(3, levels);
		}
		
		else if(size == 7) {
			int[] levels = {3, 2, 2};
			return new CombiGrid(3, levels);
		}
		
		else if(size == 8) {
			int[] levels = {2, 2, 2, 2};
			return new CombiGrid(4, levels);
		}
		
		else if(size == 9) {
			int[] levels = {3, 2, 2, 2};
			return new CombiGrid(4, levels);
		}
		
		else if(size == 10) {
			int[] levels = {3, 3, 2, 2, 2};
			return new CombiGrid(5, levels);
		}
		
		/*else if(size >= 19) {
			int f = size - 5 - 4 - 3 - 2;
			int[] levels = {f, 6, 4, 3, 2};
			return new CombiGrid(5, levels);
		}*/
		
		else {
			int dim = size / 5;
			int rem = size % 5;
			int[] levels = {dim + rem, dim, dim, dim, dim};
			return new CombiGrid(5, levels);
		}
	}
}
