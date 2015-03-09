package testing;
import grid.*;
import gridFunctions.*;

import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.util.List;

public class ScalingTest {

	public static void main(String[] args) {
		//int[] levels = {10, 6, 4, 3, 2};
		int[] levels = {5, 5, 5, 5, 5};
		//CombiGrid cg = new CombiGrid(levels);
		//CombiGridAligned cga = new CombiGridAligned(levels, 32);
		//testAllUnoptimized();
		//testAllOptimized();
		//test(cg, 10);
		//testOptimized(cga, 10);
		//testThreads(cg, 10);
		//testOptimizedThreads(cga, 10);
		//testThreadsOnce(cg, 10);
		//testOptimizedThreadsOnce(cga, 10);
		//testTasks(cg, 10);
		//testOptimizedTasks(cga, 10);
		//testParallelStream(cg, 10);
		//testParallelStreamChunks(cg, 10);
		//testOptimizedParallelStreamChunks(cga, 10);
	}
	
	private static void testAllUnoptimized() {
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
		
		System.out.println("Testing all Unoptimized");
		System.out.println("Gridsize" + '\t' + "Sequentiel" + '\t' + "Threads" + '\t' + "ThreadsOnce" + '\t' + "Tasks" + '\t' + "ParallelStream" + '\t' + "ParallelStreamChunks");
		
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
				cg.hierarchizeUnoptimizedThreads(4);
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
				cg.hierarchizeUnoptimizedTasks(4);
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
				cg.hierarchizeUnoptimizedParallelStream(4);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeParallelStreamBlocks = totalTime / 10;
			System.out.println("" + size + '\t' + avgTime + '\t' + avgTimeThreads + '\t' + avgTimeThreadsOnce + '\t' + avgTimeTasks + '\t' + avgTimeParallelStream + '\t' + avgTimeParallelStreamBlocks);
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
		double avgTimeParallelStreamBlocks;
		
		System.out.println("Testing all Optimized");
		System.out.println("Gridsize" + '\t' + "Sequentiel" + '\t' + "Threads" + '\t' + "ThreadsOnce" + '\t' + "Tasks" + '\t' + "ParallelStreamChunks");
		
		for(int size = 19; size <= 26; size++) {
			CombiGridAligned cg = createAlignedGrid(size);
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimized(4);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreads(4, 4);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeThreads = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreadsOnce(4, 4);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeThreadsOnce = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedTasks(4, 4);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeTasks = totalTime / 10;
			
			totalTime = 0;
			for(int i = 0; i < 10; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedParallelStream(4, 4);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTimeParallelStreamBlocks = totalTime / 10;
			
			System.out.println("" + size + '\t' + avgTime + '\t' + avgTimeThreads + '\t' + avgTimeThreadsOnce + '\t' + avgTimeTasks + '\t' + avgTimeParallelStreamBlocks);
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
		
		//Run 10 times
		for(int run = 1; run <= 16; run++) {
			totalTime = 0;
			//Run 10 times to average result
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
	
	private static void testOptimized(CombiGridAligned cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimized");
		System.out.println("Run" + '\t' + "Time in ms");
		
		//Run 10 times
		for(int run = 1; run <= 16; run++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimized(16);
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
		
		//Run for 1-10 threads
		for(int threads = 1; threads <= 16; threads++) {
			totalTime = 0;
			//Run 10 times to average result
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
	
	private static void testOptimizedThreads(CombiGridAligned cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimizedThreads");
		System.out.println("Threads" + '\t' + "Time in ms");
		
		//Run for 1-10 threads
		for(int threads = 1; threads <= 16; threads++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreads(16, threads);
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
		
		//Run for 1-10 threads
		for(int threads = 1; threads <= 16; threads++) {
			totalTime = 0;
			//Run 10 times to average result
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
		
		//Run for 1-10 threads
		for(int threads = 1; threads <= 16; threads++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreadsOnce(16, threads);
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
		
		//Run for 1-100 tasks
		for(int tasks = 1; tasks <= 16; tasks++) {
			totalTime = 0;
			//Run 10 times to average result
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
		
		//Run for 1-100 tasks
		for(int tasks = 1; tasks <= 16; tasks++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedTasks(16, tasks);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + tasks + '\t' + avgTime);
		}
	}
	
	private static void testParallelStream(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedParallelStream");
		System.out.println("Run" + '\t' + "Time in ms");
		
		//Run 10 tests
		for(int parallelStream = 1; parallelStream <= 16; parallelStream++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedParallelStream();
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + parallelStream + '\t' + avgTime);
		}
	}
	
	private static void testParallelStreamChunks(CombiGrid cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeUnoptimizedParallelStreamChunks");
		System.out.println("Chunks" + '\t' + "Time in ms");
		
		//Run for 1-100 chunks
		for(int chunks = 1; chunks <= 16; chunks++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedParallelStream(chunks);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + chunks + '\t' + avgTime);
		}
	}
	
	private static void testOptimizedParallelStreamChunks(CombiGridAligned cg, int repititions) {
		long start;
		long end;
		double time;
		double totalTime;
		double avgTime;

		System.out.println("HierarchizeOptimizedParallelStreamChunks");
		System.out.println("Chunks" + '\t' + "Time in ms");
		
		//Run for 1-100 chunks
		for(int chunks = 1; chunks <= 16; chunks++) {
			totalTime = 0;
			//Run 10 times to average result
			for(int i = 0; i < repititions; i++) {
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedParallelStream(16, chunks);
				end = System.currentTimeMillis();
				time = end - start;
				totalTime += time;
			}
			avgTime = totalTime / repititions;
			System.out.println("" + chunks + '\t' + avgTime);
		}
	}
	
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
		
		else if(size == 10) {
			int[] levels = {3, 3, 2, 2, 2};
			return new CombiGrid(levels);
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
		
		else if(size == 10) {
			int[] levels = {3, 3, 2, 2, 2};
			return new CombiGridAligned(levels, 32);
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
