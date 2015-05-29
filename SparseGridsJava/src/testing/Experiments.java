package testing;

import grid.*;
import gridFunctions.GridFunctions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class Experiments {

	static String folderName = "";
	/**
	 * First parameter is the name of the folder to put the results in, this should be experiment specific
	 * Second parameter is max size. This must be at least 10, third parameter is max threads
	 * Fourth (optional) parameter is minSize, if undeclared will be set at maxSize - 1
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length < 3)
			throw new IllegalArgumentException("Declare experiment folder name, max size and max threads");
		folderName = args[0];
		int maxSize = Integer.parseInt(args[1]);
		int threads = Integer.parseInt(args[2]);
		int minSize;
		if(args.length > 3)
			minSize = Integer.parseInt(args[3]);
		else
			minSize = maxSize - 1;
		
		
		scalingExperiment(maxSize, 5, 10, 1, threads, true);
		scalingExperiment(maxSize, 5, 10, 1, threads, false);
		varyingSizeExperiment(minSize, maxSize, 10, threads, true);
		varyingSizeExperiment(minSize, maxSize, 10, threads, false);
	}

	public static void scalingExperiment(int size, int dim, int rep, int min, int max, boolean iso) {
		UnoptimizedThreadsScaling(size, dim, rep, min, max, iso);
		UnoptimizedThreadsOnceScaling(size, dim, rep, min, max, iso);
		UnoptimizedTasksScaling(size, dim, rep, min, max, iso);
		OptimizedThreadsScaling(size, dim, rep, min, max, iso);
		OptimizedThreadsNoUnrollScaling(size, dim, rep, min, max, iso);
		AlignedThreadsScaling(size, dim, rep, min, max, iso);
		AlignedThreadsOnceScaling(size, dim, rep, min, max, iso);
		AlignedTasksScaling(size, dim, rep, min, max, iso);
		RecursiveScaling(size, dim, rep, min, max, iso);
	}

	public static void varyingSizeExperiment(int minSize, int maxSize, int rep, int threads, boolean iso) {
		UnoptimizedThreadsVaryingSize(minSize, maxSize, rep, threads, iso);
		UnoptimizedThreadsOnceVaryingSize(minSize, maxSize, rep, threads, iso);
		UnoptimizedTasksVaryingSize(minSize, maxSize, rep, threads, iso);
		OptimizedThreadsVaryingSize(minSize, maxSize, rep, threads, iso);
		OptimizedThreadsNoUnrollVaryingSize(minSize, maxSize, rep, threads, iso);
		AlignedThreadsVaryingSize(minSize, maxSize, rep, threads, iso);
		AlignedThreadsOnceVaryingSize(minSize, maxSize, rep, threads, iso);
		AlignedTasksVaryingSize(minSize, maxSize, rep, threads, iso);
		RecursiveVaryingSize(minSize, maxSize, rep, threads, iso);
	}

	private static void UnoptimizedThreadsScaling(int size, int dimensions, int repititions, int minNumberOfThreads, int maxNumberOfThreads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeUnoptimizedThreads(16);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		if(isotropic)
			cg = CombiGridBuilder.isotropicGrid(size, dimensions);
		else
			cg = CombiGridBuilder.anisotropicGrid(size, dimensions);
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeUnoptimizedThreads");
		data.add("Grid: " + cg.getLevels());
		data.add("Threads" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		for(int thread = minNumberOfThreads; thread <= maxNumberOfThreads; thread++) {
			double[] times = new double[repititions];
			dev = 0;
			minTime = Double.MAX_VALUE;
			maxTime = Double.MIN_VALUE;
			for(int i = 0; i < repititions; i++) {
				System.gc();
				try {Thread.sleep(1000);} catch (InterruptedException e) {}
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedThreads(thread);
				end = System.currentTimeMillis();
				time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times[i] = time;
			}

			Arrays.sort(times);
			avgTime = (times[4] + times[5]) / 2;
			
			//Calculate standard deviation
			for(double d : times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));			
			
			data.add("" + thread + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}
		writeToFile("UnoptimizedThreads Scaling " + cg.getLevels(), data);
	}

	private static void UnoptimizedThreadsOnceScaling(int size, int dimensions, int repititions, int minNumberOfThreads, int maxNumberOfThreads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeUnoptimizedThreadsOnce(16);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		if(isotropic)
			cg = CombiGridBuilder.isotropicGrid(size, dimensions);
		else
			cg = CombiGridBuilder.anisotropicGrid(size, dimensions);
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeUnoptimizedThreadsOnce");
		data.add("Grid: " + cg.getLevels());
		data.add("Threads" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		for(int thread = minNumberOfThreads; thread <= maxNumberOfThreads; thread++) {
			double[] times = new double[repititions];
			dev = 0;
			minTime = Double.MAX_VALUE;
			maxTime = Double.MIN_VALUE;
			for(int i = 0; i < repititions; i++) {
				System.gc();
				try {Thread.sleep(1000);} catch (InterruptedException e) {}
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedThreadsOnce(thread);
				end = System.currentTimeMillis();
				time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times[i] = time;
			}

			Arrays.sort(times);
			avgTime = (times[4] + times[5]) / 2;

			//Calculate standard deviation
			for(double d : times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));

			data.add("" + thread + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}
		writeToFile("UnoptimizedThreadsOnce Scaling " + cg.getLevels(), data);
	}

	private static void UnoptimizedTasksScaling(int size, int dimensions, int repititions, int minNumberOfTasks, int maxNumberOfTasks, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//Reduced warmup for task based hierarchization, for some it consumes a lot of the machine's resources.
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeUnoptimizedTasks(1);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		if(isotropic)
			cg = CombiGridBuilder.isotropicGrid(size, dimensions);
		else
			cg = CombiGridBuilder.anisotropicGrid(size, dimensions);
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeUnoptimizedTasks");
		data.add("Grid: " + cg.getLevels());
		data.add("Tasks" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		for(int task = minNumberOfTasks; task <= maxNumberOfTasks; task++) {
			double[] times = new double[repititions];
			dev = 0;
			minTime = Double.MAX_VALUE;
			maxTime = Double.MIN_VALUE;
			for(int i = 0; i < repititions; i++) {
				System.gc();
				try {Thread.sleep(1000);} catch (InterruptedException e) {}
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeUnoptimizedTasks(task);
				end = System.currentTimeMillis();
				time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times[i] = time;
			}

			Arrays.sort(times);
			avgTime = (times[4] + times[5]) / 2;

			//Calculate standard deviation
			for(double d : times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));

			data.add("" + task + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}
		writeToFile("UnoptimizedTasks Scaling " + cg.getLevels(), data);
	}

	private static void OptimizedThreadsScaling(int size, int dimensions, int repititions, int minNumberOfThreads, int maxNumberOfThreads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedThreads(16);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		if(isotropic)
			cg = CombiGridBuilder.isotropicGrid(size, dimensions);
		else
			cg = CombiGridBuilder.anisotropicGrid(size, dimensions);
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeOptimizedThreads");
		data.add("Grid: " + cg.getLevels());
		data.add("Threads" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		for(int thread = minNumberOfThreads; thread <= maxNumberOfThreads; thread++) {
			double[] times = new double[repititions];
			dev = 0;
			minTime = Double.MAX_VALUE;
			maxTime = Double.MIN_VALUE;
			for(int i = 0; i < repititions; i++) {
				System.gc();
				try {Thread.sleep(1000);} catch (InterruptedException e) {}
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreads(thread);
				end = System.currentTimeMillis();
				time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times[i] = time;
			}

			Arrays.sort(times);
			avgTime = (times[4] + times[5]) / 2;

			//Calculate standard deviation
			for(double d : times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));

			data.add("" + thread + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}
		writeToFile("OptimizedThreads Scaling " + cg.getLevels(), data);
	}

	private static void OptimizedThreadsNoUnrollScaling(int size, int dimensions, int repititions, int minNumberOfThreads, int maxNumberOfThreads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedThreadsNoUnroll(16);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		if(isotropic)
			cg = CombiGridBuilder.isotropicGrid(size, dimensions);
		else
			cg = CombiGridBuilder.anisotropicGrid(size, dimensions);
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeOptimizedThreadsNoUnroll");
		data.add("Grid: " + cg.getLevels());
		data.add("Threads" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		for(int thread = minNumberOfThreads; thread <= maxNumberOfThreads; thread++) {
			double[] times = new double[repititions];
			dev = 0;
			minTime = Double.MAX_VALUE;
			maxTime = Double.MIN_VALUE;
			for(int i = 0; i < repititions; i++) {
				System.gc();
				try {Thread.sleep(1000);} catch (InterruptedException e) {}
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreadsNoUnroll(thread);
				end = System.currentTimeMillis();
				time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times[i] = time;
			}

			Arrays.sort(times);
			avgTime = (times[4] + times[5]) / 2;

			//Calculate standard deviation
			for(double d : times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));

			data.add("" + thread + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}
		writeToFile("OptimizedThreadsNoUnroll Scaling " + cg.getLevels(), data);
	}

	private static void AlignedThreadsScaling(int size, int dimensions, int repititions, int minNumberOfThreads, int maxNumberOfThreads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGridAligned grid = CombiGridBuilder.isotropicAlignedGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedThreads(16);
		}
		grid = null;
		System.gc();

		CombiGridAligned cg;
		if(isotropic)
			cg = CombiGridBuilder.isotropicAlignedGrid(size, dimensions);
		else
			cg = CombiGridBuilder.anisotropicAlignedGrid(size, dimensions);
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeAlignedOptimizedThreads");
		data.add("Grid: " + cg.getLevels());
		data.add("Threads" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		for(int thread = minNumberOfThreads; thread <= maxNumberOfThreads; thread++) {
			double[] times = new double[repititions];
			dev = 0;
			minTime = Double.MAX_VALUE;
			maxTime = Double.MIN_VALUE;
			for(int i = 0; i < repititions; i++) {
				System.gc();
				try {Thread.sleep(1000);} catch (InterruptedException e) {}
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreads(thread);
				end = System.currentTimeMillis();
				time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times[i] = time;
			}

			Arrays.sort(times);
			avgTime = (times[4] + times[5]) / 2;

			//Calculate standard deviation
			for(double d : times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));

			data.add("" + thread + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}
		writeToFile("Aligned optimizedThreads Scaling " + cg.getLevels(), data);
	}

	private static void AlignedThreadsOnceScaling(int size, int dimensions, int repititions, int minNumberOfThreads, int maxNumberOfThreads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGridAligned grid = CombiGridBuilder.isotropicAlignedGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedThreadsOnce(16);
		}
		grid = null;
		System.gc();

		CombiGridAligned cg;
		if(isotropic)
			cg = CombiGridBuilder.isotropicAlignedGrid(size, dimensions);
		else
			cg = CombiGridBuilder.anisotropicAlignedGrid(size, dimensions);
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeAlignedOptimizedThreadsOnce");
		data.add("Grid: " + cg.getLevels());
		data.add("Threads" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		for(int thread = minNumberOfThreads; thread <= maxNumberOfThreads; thread++) {
			double[] times = new double[repititions];
			dev = 0;
			minTime = Double.MAX_VALUE;
			maxTime = Double.MIN_VALUE;
			for(int i = 0; i < repititions; i++) {
				System.gc();
				try {Thread.sleep(1000);} catch (InterruptedException e) {}
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedThreadsOnce(thread);
				end = System.currentTimeMillis();
				time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times[i] = time;
			}
			
			Arrays.sort(times);
			avgTime = (times[4] + times[5]) / 2;

			//Calculate standard deviation
			for(double d : times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));

			data.add("" + thread + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}
		writeToFile("Aligned optimizedThreadsOnce Scaling " + cg.getLevels(), data);
	}

	private static void AlignedTasksScaling(int size, int dimensions, int repititions, int minNumberOfTasks, int maxNumberOfTasks, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//Reduced warmup as the task based methods consume a lot of memory for some reason
		CombiGridAligned grid = CombiGridBuilder.isotropicAlignedGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedTasks(1);
		}
		grid = null;
		System.gc();

		CombiGridAligned cg;
		if(isotropic)
			cg = CombiGridBuilder.isotropicAlignedGrid(size, dimensions);
		else
			cg = CombiGridBuilder.anisotropicAlignedGrid(size, dimensions);
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeAlignedOptimizedTasks");
		data.add("Grid: " + cg.getLevels());
		data.add("Threads" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		for(int tasks = minNumberOfTasks; tasks <= maxNumberOfTasks; tasks++) {
			double[] times = new double[repititions];
			dev = 0;
			minTime = Double.MAX_VALUE;
			maxTime = Double.MIN_VALUE;
			for(int i = 0; i < repititions; i++) {
				System.gc();
				try {Thread.sleep(1000);} catch (InterruptedException e) {}
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeOptimizedTasks(tasks);
				end = System.currentTimeMillis();
				time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times[i] = time;
			}

			Arrays.sort(times);
			avgTime = (times[4] + times[5]) / 2;

			//Calculate standard deviation
			for(double d : times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));

			data.add("" + tasks + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}
		writeToFile("Aligned optimizedTasks Scaling " + cg.getLevels(), data);
	}

	private static void RecursiveScaling(int size, int dimensions, int repititions, int minNumberOfThreads, int maxNumberOfThreads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGridAligned grid = CombiGridBuilder.isotropicAlignedGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeRecursiveThreads(3);
		}
		grid = null;
		System.gc();

		CombiGridAligned cg;
		if(isotropic)
			cg = CombiGridBuilder.isotropicAlignedGrid(size, dimensions);
		else
			cg = CombiGridBuilder.anisotropicAlignedGrid(size, dimensions);
		cg.recTile = cg.levels[0];
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeRecursiveThreads");
		data.add("Grid: " + cg.getLevels());
		data.add("RecLevel" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		for(int thread = minNumberOfThreads; thread <= maxNumberOfThreads; thread++) {
			double[] times = new double[repititions];
			dev = 0;
			minTime = Double.MAX_VALUE;
			maxTime = Double.MIN_VALUE;
			for(int i = 0; i < repititions; i++) {
				System.gc();
				try {Thread.sleep(1000);} catch (InterruptedException e) {}
				cg.setValues(GridFunctions.ALLONES);
				start = System.currentTimeMillis();
				cg.hierarchizeRecursiveThreads(thread);
				end = System.currentTimeMillis();
				time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times[i] = time;
			}

			Arrays.sort(times);
			avgTime = (times[4] + times[5]) / 2;

			//Calculate standard deviation
			for(double d : times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));

			data.add("" + thread + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}
		writeToFile("Aligned RecursiveThreads Scaling " + cg.getLevels(), data);
	}

	private static void UnoptimizedThreadsVaryingSize(int minSize, int maxSize, int repititions, int threads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeUnoptimizedThreads(16);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeUnoptimizedThreads");
		data.add("Grid" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		//Loop over dimensions
		for(int dim = 2; dim <= 5; dim++) {
			for(int size = minSize; size <= maxSize; size++) {
				if(isotropic)
					cg = CombiGridBuilder.isotropicGrid(size, dim);
				else
					cg = CombiGridBuilder.anisotropicGrid(size, dim);
				double[] times = new double[repititions];
				dev = 0;
				minTime = Double.MAX_VALUE;
				maxTime = Double.MIN_VALUE;
				for(int i = 0; i < repititions; i++) {
					System.gc();
					try {Thread.sleep(1000);} catch (InterruptedException e) {}
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeUnoptimizedThreads(threads);
					end = System.currentTimeMillis();
					time = end - start;
					if(time < minTime)
						minTime = time;
					if(time > maxTime)
						maxTime = time;
					times[i] = time;
				}

				Arrays.sort(times);
				avgTime = (times[4] + times[5]) / 2;

				//Calculate standard deviation
				for(double d : times) {
					dev += Math.pow(d - avgTime, 2);
				}
				dev = Math.sqrt(dev / (times.length - 1));

				data.add("" + cg.getLevels() + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
			}
		}
		String form;
		if(isotropic)
			form = "Isotropic";
		else
			form = "Anisotropic";
		writeToFile("UnoptimizedThreads VaryingSize " + minSize + " " + maxSize + " " + form, data);
	}
	
	private static void UnoptimizedThreadsOnceVaryingSize(int minSize, int maxSize, int repititions, int threads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeUnoptimizedThreadsOnce(16);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeUnoptimizedThreadsOnce");
		data.add("Grid" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		//Loop over dimensions
		for(int dim = 2; dim <= 5; dim++) {
			for(int size = minSize; size <= maxSize; size++) {
				if(isotropic)
					cg = CombiGridBuilder.isotropicGrid(size, dim);
				else
					cg = CombiGridBuilder.anisotropicGrid(size, dim);
				double[] times = new double[repititions];
				dev = 0;
				minTime = Double.MAX_VALUE;
				maxTime = Double.MIN_VALUE;
				for(int i = 0; i < repititions; i++) {
					System.gc();
					try {Thread.sleep(1000);} catch (InterruptedException e) {}
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeUnoptimizedThreadsOnce(threads);
					end = System.currentTimeMillis();
					time = end - start;
					if(time < minTime)
						minTime = time;
					if(time > maxTime)
						maxTime = time;
					times[i] = time;
				}

				Arrays.sort(times);
				avgTime = (times[4] + times[5]) / 2;

				//Calculate standard deviation
				for(double d : times) {
					dev += Math.pow(d - avgTime, 2);
				}
				dev = Math.sqrt(dev / (times.length - 1));

				data.add("" + cg.getLevels() + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
			}
		}
		String form;
		if(isotropic)
			form = "Isotropic";
		else
			form = "Anisotropic";
		writeToFile("UnoptimizedThreadsOnce VaryingSize " + minSize + " " + maxSize + " " + form, data);
	}
	
	private static void UnoptimizedTasksVaryingSize(int minSize, int maxSize, int repititions, int tasks, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//Reduced warmup as the task methods use a lot of memory if called many times
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 100; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeUnoptimizedTasks(1);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeUnoptimizedTasks");
		data.add("Grid" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		//Loop over dimensions
		for(int dim = 2; dim <= 5; dim++) {
			for(int size = minSize; size <= maxSize; size++) {
				if(isotropic)
					cg = CombiGridBuilder.isotropicGrid(size, dim);
				else
					cg = CombiGridBuilder.anisotropicGrid(size, dim);
				double[] times = new double[repititions];
				dev = 0;
				minTime = Double.MAX_VALUE;
				maxTime = Double.MIN_VALUE;
				for(int i = 0; i < repititions; i++) {
					System.gc();
					try {Thread.sleep(1000);} catch (InterruptedException e) {}
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeUnoptimizedTasks(tasks);
					end = System.currentTimeMillis();
					time = end - start;
					if(time < minTime)
						minTime = time;
					if(time > maxTime)
						maxTime = time;
					times[i] = time;
				}

				Arrays.sort(times);
				avgTime = (times[4] + times[5]) / 2;

				//Calculate standard deviation
				for(double d : times) {
					dev += Math.pow(d - avgTime, 2);
				}
				dev = Math.sqrt(dev / (times.length - 1));

				data.add("" + cg.getLevels() + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
			}
		}
		String form;
		if(isotropic)
			form = "Isotropic";
		else
			form = "Anisotropic";
		writeToFile("UnoptimizedTasks VaryingSize " + minSize + " " + maxSize + " " + form, data);
	}
	
	private static void OptimizedThreadsVaryingSize(int minSize, int maxSize, int repititions, int threads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedThreads(16);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeOptimizedThreads");
		data.add("Grid" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		//Loop over dimensions
		for(int dim = 2; dim <= 5; dim++) {
			for(int size = minSize; size <= maxSize; size++) {
				if(isotropic)
					cg = CombiGridBuilder.isotropicGrid(size, dim);
				else
					cg = CombiGridBuilder.anisotropicGrid(size, dim);
				double[] times = new double[repititions];
				dev = 0;
				minTime = Double.MAX_VALUE;
				maxTime = Double.MIN_VALUE;
				for(int i = 0; i < repititions; i++) {
					System.gc();
					try {Thread.sleep(1000);} catch (InterruptedException e) {}
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeOptimizedThreads(threads);
					end = System.currentTimeMillis();
					time = end - start;
					if(time < minTime)
						minTime = time;
					if(time > maxTime)
						maxTime = time;
					times[i] = time;
				}

				Arrays.sort(times);
				avgTime = (times[4] + times[5]) / 2;

				//Calculate standard deviation
				for(double d : times) {
					dev += Math.pow(d - avgTime, 2);
				}
				dev = Math.sqrt(dev / (times.length - 1));

				data.add("" + cg.getLevels() + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
			}
		}
		String form;
		if(isotropic)
			form = "Isotropic";
		else
			form = "Anisotropic";
		writeToFile("OptimizedThreads VaryingSize " + minSize + " " + maxSize + " " + form, data);
	}
	
	private static void OptimizedThreadsNoUnrollVaryingSize(int minSize, int maxSize, int repititions, int threads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGrid grid = CombiGridBuilder.isotropicGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedThreadsNoUnroll(16);
		}
		grid = null;
		System.gc();

		CombiGrid cg;
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeOptimizedThreadsNoUnroll");
		data.add("Grid" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		//Loop over dimensions
		for(int dim = 2; dim <= 5; dim++) {
			for(int size = minSize; size <= maxSize; size++) {
				if(isotropic)
					cg = CombiGridBuilder.isotropicGrid(size, dim);
				else
					cg = CombiGridBuilder.anisotropicGrid(size, dim);
				double[] times = new double[repititions];
				dev = 0;
				minTime = Double.MAX_VALUE;
				maxTime = Double.MIN_VALUE;
				for(int i = 0; i < repititions; i++) {
					System.gc();
					try {Thread.sleep(1000);} catch (InterruptedException e) {}
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeOptimizedThreadsNoUnroll(threads);
					end = System.currentTimeMillis();
					time = end - start;
					if(time < minTime)
						minTime = time;
					if(time > maxTime)
						maxTime = time;
					times[i] = time;
				}

				Arrays.sort(times);
				avgTime = (times[4] + times[5]) / 2;

				//Calculate standard deviation
				for(double d : times) {
					dev += Math.pow(d - avgTime, 2);
				}
				dev = Math.sqrt(dev / (times.length - 1));

				data.add("" + cg.getLevels() + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
			}
		}
		String form;
		if(isotropic)
			form = "Isotropic";
		else
			form = "Anisotropic";
		writeToFile("OptimizedThreadsNoUnroll VaryingSize " + minSize + " " + maxSize + " " + form, data);
	}
	
	private static void AlignedThreadsVaryingSize(int minSize, int maxSize, int repititions, int threads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGridAligned grid = CombiGridBuilder.isotropicAlignedGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedThreads(16);
		}
		grid = null;
		System.gc();

		CombiGridAligned cg;
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeOptimizedThreads");
		data.add("Grid" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		//Loop over dimensions
		for(int dim = 2; dim <= 5; dim++) {
			for(int size = minSize; size <= maxSize; size++) {
				if(isotropic)
					cg = CombiGridBuilder.isotropicAlignedGrid(size, dim);
				else
					cg = CombiGridBuilder.anisotropicAlignedGrid(size, dim);
				double[] times = new double[repititions];
				dev = 0;
				minTime = Double.MAX_VALUE;
				maxTime = Double.MIN_VALUE;
				for(int i = 0; i < repititions; i++) {
					System.gc();
					try {Thread.sleep(1000);} catch (InterruptedException e) {}
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeOptimizedThreads(threads);
					end = System.currentTimeMillis();
					time = end - start;
					if(time < minTime)
						minTime = time;
					if(time > maxTime)
						maxTime = time;
					times[i] = time;
				}

				Arrays.sort(times);
				avgTime = (times[4] + times[5]) / 2;

				//Calculate standard deviation
				for(double d : times) {
					dev += Math.pow(d - avgTime, 2);
				}
				dev = Math.sqrt(dev / (times.length - 1));

				data.add("" + cg.getLevels() + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
			}
		}
		String form;
		if(isotropic)
			form = "Isotropic";
		else
			form = "Anisotropic";
		writeToFile("Aligned OptimizedThreads VaryingSize " + minSize + " " + maxSize + " " + form, data);
	}

	private static void AlignedThreadsOnceVaryingSize(int minSize, int maxSize, int repititions, int threads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGridAligned grid = CombiGridBuilder.isotropicAlignedGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedThreadsOnce(16);
		}
		grid = null;
		System.gc();

		CombiGridAligned cg;
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeOptimizedThreadsOnce");
		data.add("Grid" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		//Loop over dimensions
		for(int dim = 2; dim <= 5; dim++) {
			for(int size = minSize; size <= maxSize; size++) {
				if(isotropic)
					cg = CombiGridBuilder.isotropicAlignedGrid(size, dim);
				else
					cg = CombiGridBuilder.anisotropicAlignedGrid(size, dim);
				double[] times = new double[repititions];
				dev = 0;
				minTime = Double.MAX_VALUE;
				maxTime = Double.MIN_VALUE;
				for(int i = 0; i < repititions; i++) {
					System.gc();
					try {Thread.sleep(1000);} catch (InterruptedException e) {}
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeOptimizedThreadsOnce(threads);
					end = System.currentTimeMillis();
					time = end - start;
					if(time < minTime)
						minTime = time;
					if(time > maxTime)
						maxTime = time;
					times[i] = time;
				}

				Arrays.sort(times);
				avgTime = (times[4] + times[5]) / 2;

				//Calculate standard deviation
				for(double d : times) {
					dev += Math.pow(d - avgTime, 2);
				}
				dev = Math.sqrt(dev / (times.length - 1));

				data.add("" + cg.getLevels() + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
			}
		}
		String form;
		if(isotropic)
			form = "Isotropic";
		else
			form = "Anisotropic";
		writeToFile("Aligned OptimizedThreadsOnce VaryingSize " + minSize + " " + maxSize + " " + form, data);
	}
	
	private static void AlignedTasksVaryingSize(int minSize, int maxSize, int repititions, int tasks, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//Reduced warmup as the task based methods use a lot of memory if called in succession
		CombiGridAligned grid = CombiGridBuilder.isotropicAlignedGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeOptimizedTasks(16);
		}
		grid = null;
		System.gc();

		CombiGridAligned cg;
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeOptimizedTasks");
		data.add("Grid" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		//Loop over dimensions
		for(int dim = 2; dim <= 5; dim++) {
			for(int size = minSize; size <= maxSize; size++) {
				if(isotropic)
					cg = CombiGridBuilder.isotropicAlignedGrid(size, dim);
				else
					cg = CombiGridBuilder.anisotropicAlignedGrid(size, dim);
				double[] times = new double[repititions];
				dev = 0;
				minTime = Double.MAX_VALUE;
				maxTime = Double.MIN_VALUE;
				for(int i = 0; i < repititions; i++) {
					System.gc();
					try {Thread.sleep(1000);} catch (InterruptedException e) {}
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeOptimizedTasks(tasks);
					end = System.currentTimeMillis();
					time = end - start;
					if(time < minTime)
						minTime = time;
					if(time > maxTime)
						maxTime = time;
					times[i] = time;
				}

				Arrays.sort(times);
				avgTime = (times[4] + times[5]) / 2;

				//Calculate standard deviation
				for(double d : times) {
					dev += Math.pow(d - avgTime, 2);
				}
				dev = Math.sqrt(dev / (times.length - 1));

				data.add("" + cg.getLevels() + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
			}
		}
		String form;
		if(isotropic)
			form = "Isotropic";
		else
			form = "Anisotropic";
		writeToFile("Aligned OptimizedTasks VaryingSize " + minSize + " " + maxSize + " " + form, data);
	}
	
	private static void RecursiveVaryingSize(int minSize, int maxSize, int repititions, int threads, boolean isotropic) {
		long start;
		long end;
		double time;
		double avgTime;
		double minTime;
		double maxTime;
		double dev = 0;

		//warmup
		CombiGridAligned grid = CombiGridBuilder.isotropicAlignedGrid(15, 3);
		for(int i = 0; i < 10000; i++) {
			grid.setValues(GridFunctions.ALLONES);
			grid.hierarchizeRecursiveThreads(3);
		}
		grid = null;
		System.gc();

		CombiGridAligned cg;
		List<String> data = new LinkedList<String>();
		data.add(PCInfo());
		data.add("HierarchizeRecursiveThread");
		data.add("Grid" + '\t' + "Mean" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev");
		//Loop over dimensions
		for(int dim = 2; dim <= 5; dim++) {
			for(int size = minSize; size <= maxSize; size++) {
				if(isotropic)
					cg = CombiGridBuilder.isotropicAlignedGrid(size, dim);
				else
					cg = CombiGridBuilder.anisotropicAlignedGrid(size, dim);
				cg.recTile = cg.levels[0];
				double[] times = new double[repititions];
				dev = 0;
				minTime = Double.MAX_VALUE;
				maxTime = Double.MIN_VALUE;
				for(int i = 0; i < repititions; i++) {
					System.gc();
					try {Thread.sleep(1000);} catch (InterruptedException e) {}
					cg.setValues(GridFunctions.ALLONES);
					start = System.currentTimeMillis();
					cg.hierarchizeRecursiveThreads(3);
					end = System.currentTimeMillis();
					time = end - start;
					if(time < minTime)
						minTime = time;
					if(time > maxTime)
						maxTime = time;
					times[i] = time;
				}

				Arrays.sort(times);
				avgTime = (times[4] + times[5]) / 2;

				//Calculate standard deviation
				for(double d : times) {
					dev += Math.pow(d - avgTime, 2);
				}
				dev = Math.sqrt(dev / (times.length - 1));

				data.add("" + cg.getLevels() + '\t' + avgTime + '\t' + minTime + '\t' + maxTime + '\t' + dev);
			}
		}
		String form;
		if(isotropic)
			form = "Isotropic";
		else
			form = "Anisotropic";
		writeToFile("Aligned Recursive VaryingSize " + minSize + " " + maxSize + " " + form, data);
	}

	private static String PCInfo() {
		RuntimeMXBean runtimeMxBean = ManagementFactory.getRuntimeMXBean();
		List<String> arguments = runtimeMxBean.getInputArguments();
		StringBuilder sb = new StringBuilder();
		sb.append("# OS: " +
				System.getProperty("os.name") + " " +
				System.getProperty("os.version") + " " +
				System.getProperty("os.arch"));
		sb.append('\n');
		sb.append("# JVM: " +
				System.getProperty("java.vendor") + " " +
				System.getProperty("java.version"));
		sb.append('\n');
		// This line works only on MS Windows:
		sb.append("# CPU: " + System.getenv("PROCESSOR_IDENTIFIER"));
		sb.append('\n');
		sb.append("# Available processors: " + Runtime.getRuntime().availableProcessors());
		sb.append('\n');
		sb.append("# Maximum memory: " + Runtime.getRuntime().maxMemory() / (1024 * 1024) + " MB");
		sb.append('\n');
		java.util.Date now = new java.util.Date();
		sb.append("# Date: " +
				new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ssZ").format(now));
		sb.append('\n');
		sb.append("# VM Arguments: ");
		for(String s : arguments)
			sb.append(s + " ");
		sb.append('\n');
		return sb.toString();
	}

	private static void writeToFile(String filename, List<String> data)
	{
		try {
			File dir = new File(folderName);
			if(!dir.exists())
				dir.mkdir();
			File file = new File(folderName + File.separator + filename + ".dat");

			// if file exists deletes it then creates it again
			if (file.exists()) {
				file.delete();
			}
			file.createNewFile();


			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			for(String s : data) {
				bw.write(s);
				bw.write("\n");
			}
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
