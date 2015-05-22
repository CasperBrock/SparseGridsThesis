package testingScala

import grid._
import java.lang.management.ManagementFactory
import java.lang.management.RuntimeMXBean
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter
import java.io.IOException
import scala.collection.immutable.List
import scala.collection.mutable.ListBuffer

object experiments {

  /**
   * First parameter is max size, second parameter is max threads
   * Third (optional) parameter is minSize, if undeclared will be set at maxSize - 1
   */
  var testName:String="UnsetTestName";
  

	def main(args: Array[String]) : Unit = {    
 
      //Input handling
			if(args.length < 2)
				throw new IllegalArgumentException("Declare test name and maxSize, and optionally minSize: name max min");


			var minSize:Int = 1;
			var maxSize: Int = 2;

      testName = args(0);
      maxSize = Integer.parseInt(args(1))
      if (args.length > 2)	{				//minSize = Integer.parseInt(args(0))
						minSize = Integer.parseInt(args(2))	
      } else {minSize=maxSize}
				
      
      if (minSize == maxSize) {
        println("Running test "+testName+" on single grid, of size " + minSize)
      } else 
			    println("Running test "+testName+" with minSize " + minSize + ", and maxSize " + maxSize);

			hierarchizeVarSizeExperiment(minSize, maxSize, 10, true); //run varying size on both iso and anisotropic grids.
			hierarchizeVarSizeExperiment(minSize, maxSize, 10, false);
			recursiveVarSizeExperiment(minSize, maxSize, 10, true);
			recursiveVarSizeExperiment(minSize, maxSize, 10, false);
      println("All tests finished.")
	}
  
  def recursiveVarSizeExperiment(minSize: Int, maxSize: Int, rep: Int, iso:Boolean) {
    ScalaRecursiveThreadedHierarchizeVaryingSize(minSize, maxSize, rep, iso);
    ScalaRecursiveNonThreadedVaryingSize(minSize, maxSize, rep, iso);
  }

	def hierarchizeVarSizeExperiment(minSize: Int, maxSize: Int, rep: Int, iso:Boolean) {
		ScalaNonThreadedVaryingSize(minSize, maxSize, rep, iso);
    ScalaParallelStreamVaryingSize(minSize, maxSize, rep, iso);
    ScalaHierarchizeSequentialVaryingSize(minSize, maxSize, rep, iso);
	}

	def buildLevelVector(size: Int, dimensionsInput: Int, isotropic: Boolean): Array[Int] =  {
			if(dimensionsInput < 2)
				throw new IllegalArgumentException("A grid must have atleast 2 dimensions");
			if(size / dimensionsInput < 2)
				throw new IllegalArgumentException("Size must be atleast twice the number of dimensions");
			var rem = size - (dimensionsInput * 2);
			var levels = new Array[Int](dimensionsInput);
			if (!isotropic) { //If anisotropic
				for(i <- 0 to dimensionsInput-1 by 1) {
					var used=0;
					if(rem == 1) {
						used = 1;
						rem = 0;
					}
					else if(i + 1 == dimensionsInput) {
						used = rem;
						rem = 0;
					}
					else {
						used = Math.ceil(rem / 2.0).toInt;
						rem = rem - used;
					}
					levels(i) = 2 + used;
				}
			} else { //if isotropic
				var level: Int = size / dimensionsInput;
				for(i <- 0 to dimensionsInput-1 by 1) {
					levels(i) = level;
				}
			}
			levels;
	}

	def PCInfo(): String = {
			val runtimeMxBean: RuntimeMXBean = ManagementFactory.getRuntimeMXBean();
	    var arguments = runtimeMxBean.getInputArguments().toArray();
	    val sb: StringBuilder = new StringBuilder();
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
	    val now: java.util.Date = new java.util.Date();
	    sb.append("# Date: " +
			new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ssZ").format(now));
    	sb.append('\n');
	    sb.append("# VM Arguments: ");
  	//for(s <- arguments)
	    for (i <- arguments) { //0 to arguments.size()-1 by 1) {
		    sb.append(i + " ");
		    sb.append('\n');
		    sb.toString();
	    }
	    sb.toString();
  }

	def ScalaNonThreadedVaryingSize(minSize: Int, maxSize: Int, repititions: Int, isotropic: Boolean) {
		var dev: Double = 0;
	  var data  = scala.collection.mutable.MutableList[String]();
	  data += PCInfo();
	  data += "ScalaNonThreadedVaryingSize";
	  data += "Grid" + '\t' + "Median" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev";

  	//Warmup
  	for(j <- 0 to 9999 by 1) {
	  	CombiGrid.CombiGrid(buildLevelVector(15, 2, true)); //Create
		  CombiGrid.FillCurrentGridWithOnes();                 //Fill
		  CombiGrid.hierarchizeOptimized();                  //Run
	  }
	  System.gc();

	  for(dim <- 2 to 5 by 1) {
		  for(size <- minSize to  maxSize by 1) {
			  var minTime=Long.MaxValue;
			  var maxTime=Long.MinValue;  
			  var times = new Array[Long](repititions);
			  var totalTime: Long = 0;

			  for(i <- 0 to repititions-1 by 1) {
				  System.gc();
  				try {Thread.sleep(1000);} catch {
	    			case e: InterruptedException => println("error: " + e)
			  	}
				  CombiGrid.CombiGrid(buildLevelVector(size, dim, isotropic)); //Create
				  CombiGrid.FillCurrentGridWithOnes();                 //Fill
				  var start = System.currentTimeMillis();
				  CombiGrid.hierarchizeOptimized();                  //Run
				  var end = System.currentTimeMillis();
				  var time = end - start;
				  if(time < minTime) {
					  minTime = time;
          }
				  if(time > maxTime){
  					maxTime = time;
          }
		  		times(i) = time;
			  	totalTime += time;
			  }
        val avgTime = totalTime / repititions;
        
        var median: Long=0;
        scala.util.Sorting.quickSort(times); //Sorts from smallest to largest
        if (repititions % 2 == 0) { //if even, take average over two middle values as median.
          median = (times(repititions/2) + times((repititions/2)+1))/2; //take average of the two median values, if even.
        } else { median=times(repititions/2)}
        
        
        var Time = totalTime / repititions;

			//Calculate standard deviation
        var dev: Double = 0;
		  	for(d <- times) {
			  	dev += Math.pow(d - avgTime, 2);
			  }
        
			  dev = Math.sqrt(dev / (times.length - 1));
			  data += ("" + CombiGrid.getLevelVector() + '\t' + median + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		  }    
	  }
        
	if(isotropic) {
		val form = "Isotropic";
		writeToFile("ScalaNonThreaded VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	} else {
		val form = "Anisotropic";
		writeToFile("ScalaNonThreaded VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	  }
	}

	def ScalaParallelStreamVaryingSize(minSize: Int, maxSize: Int, repititions: Int, isotropic: Boolean) {
		var dev: Double = 0;
	  var data = scala.collection.mutable.MutableList[String]();
	  data += PCInfo();
	  data += "ScalaParallelStreamVaryingSize";
	  data += "Grid" + '\t' + "Median" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev";

  	//Warmup
	  for(i <- 0 to 9999 by 1) {
		  CombiGrid.CombiGrid(buildLevelVector(15, 2, false)); //Create
		  CombiGrid.FillCurrentGridWithOnes();                 //Fill
		  CombiGrid.hierarchizeUnoptimizedParallelStream();                  //Run
	  }
	  System.gc();

	  for(dim <- 2 to 5 by 1) {
		  for(size <- minSize to  maxSize by 1) {
  			var minTime=Long.MaxValue;
	  		var maxTime=Long.MinValue;  
		  	var times: Array[Long] = new Array[Long](repititions);
			  var totalTime: Long = 0;

			  for(i <- 0 to repititions-1 by 1) {
  				System.gc();
	  			try {Thread.sleep(1000);} catch {
		    		case e: InterruptedException => println("error: " + e)
				  }
  				CombiGrid.CombiGrid(buildLevelVector(size, dim, isotropic)); //Create
	  			CombiGrid.FillCurrentGridWithOnes();                 //Fill
		  		var start = System.currentTimeMillis();
			  	CombiGrid.hierarchizeUnoptimizedParallelStream();                  //Run
				  val end = System.currentTimeMillis();
  				val time = end - start;
	  			if(time < minTime)
		  	 		minTime = time;
			  	if(time > maxTime)
		  			maxTime = time;
			  	times(i) = time;
				  totalTime += time;
			  }
  			val avgTime = totalTime / repititions;
        
        var median: Long=0;
        scala.util.Sorting.quickSort(times); //Sorts from smallest to largest
        if (repititions % 2 == 0) { //if even, take average over two middle values as median.
          median = (times(repititions/2) + times((repititions/2)+1))/2; //take average of the two median values, if even.
        } else { median=times(repititions/2)}
        

  			//Calculate standard deviation
        var dev: Double = 0;
	  		for(d <- times) {
		  		dev += Math.pow(d - avgTime, 2);
			  }
  			dev = Math.sqrt(dev / (times.length - 1));
	  		data += ("" + CombiGrid.getLevelVector() + '\t' + median + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		  }    
	  }
	  if(isotropic) {
		  val form = "Isotropic";
		  writeToFile("Scala ParallelStream VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	  } else {
		  val form = "Anisotropic";
		  writeToFile("Scala ParallelStream VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	  }
	}

	def ScalaHierarchizeSequentialVaryingSize(minSize: Int, maxSize: Int, repititions: Int, isotropic: Boolean) {
		var dev: Double = 0;
	var data = scala.collection.mutable.MutableList[String]();
	data += PCInfo();
	data += "ScalaHierarchizeUnoptimizedSequentialVaryingSize";
	data += "Grid" + '\t' + "Median" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev";

	//Warmup
	for(i <- 0 to 9999 by 1) {
		CombiGrid.CombiGrid(buildLevelVector(15, 2, false)); //Create
		CombiGrid.FillCurrentGridWithOnes();                 //Fill
		CombiGrid.hierarchizeOptimized();                  //Run
	}
	System.gc();
  
	for(dim <- 2 to 5 by 1) {
		for(size <- minSize to  maxSize by 1) {
			var minTime=Long.MaxValue;
			var maxTime=Long.MinValue;  
			var times: Array[Long] = new Array[Long](repititions);
			var totalTime: Long = 0;

			for(i <- 0 to repititions-1 by 1) {
				System.gc();
				try {Thread.sleep(1000);} catch {
				case e: InterruptedException => println("error: " + e)
				}
				CombiGrid.CombiGrid(buildLevelVector(size, dim, isotropic)); //Create
				CombiGrid.FillCurrentGridWithOnes();                 //Fill
				var start = System.currentTimeMillis();
				CombiGrid.hierarchizeOptimized();                  //Run
				val end = System.currentTimeMillis();
				val time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times(i) = time;
				totalTime += time;
			}
			val avgTime = totalTime / repititions;

      var median: Long=0;
        scala.util.Sorting.quickSort(times); //Sorts from smallest to largest
        if (repititions % 2 == 0) { //if even, take average over two middle values as median.
          median = (times(repititions/2) + times((repititions/2)+1))/2; //take average of the two median values, if even.
        } else { median=times(repititions/2)}
      
			//Calculate standard deviation
      var dev: Double = 0;
			for(d <- times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));
			data += ("" + CombiGrid.getLevelVector() + '\t' + median + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}    
	}
	if(isotropic) {
		val form = "Isotropic";
		writeToFile("Scala hierarchizeUnoptimized VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	} else {
		val form = "Anisotropic";
		writeToFile("Scala hierarchizeUnoptimized VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	}
	}

	def ScalaRecursiveThreadedHierarchizeVaryingSize(minSize: Int, maxSize: Int, repititions: Int, isotropic: Boolean) {
		
	var data = scala.collection.mutable.MutableList[String]();
	data += PCInfo();
	data += "ScalaRecursiveThreadedHierarchizeVaryingSize";
	data += "Grid" + '\t' + "Median" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev";

	//Warmup
	for(i <- 0 to 9999 by 1) {
		CombiGrid.CombiGrid(buildLevelVector(15, 2, false)); //Create
		CombiGrid.FillCurrentGridWithOnes();                 //Fill
		CombiGrid.hierarchizeRecursiveThreaded();                  //Run
	}
	System.gc();
  
	for(dim <- 2 to 5 by 1) {
		for(size <- minSize to  maxSize by 1) {
			var minTime=Long.MaxValue;
			var maxTime=Long.MinValue;  
			var times: Array[Long] = new Array[Long](repititions);
			var totalTime: Long = 0;

			for(i <- 0 to repititions-1 by 1) {
        
				System.gc();
				try {Thread.sleep(1000);} catch {
				case e: InterruptedException => println("error: " + e)
				}
				CombiGrid.CombiGrid(buildLevelVector(size, dim, isotropic)); //Create
				CombiGrid.FillCurrentGridWithOnes();                 //Fill
				var start = System.currentTimeMillis();
				CombiGrid.hierarchizeRecursiveThreaded();                  //Run
				val end = System.currentTimeMillis();
				val time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times(i) = time;
				totalTime += time;
       
			}
			val avgTime = totalTime / repititions;
      
      var median: Long=0;
      scala.util.Sorting.quickSort(times); //Sorts from smallest to largest
      if (repititions % 2 == 0) { //if even, take average over two middle values as median.
         median = (times(repititions/2) + times((repititions/2)+1))/2; //take average of the two median values, if even.
      } else { median=times(repititions/2)}
        
			//Calculate standard deviation
      var dev: Double = 0;
			for(d <- times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));
			data += ("" + CombiGrid.getLevelVector() + '\t' + median + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}    
	}
	if(isotropic) {
		val form = "Isotropic";
		writeToFile("Scala hierarchizeRecursiveThreaded VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	} else {
		val form = "Anisotropic";
		writeToFile("Scala hierarchizeRecursiveThreaded VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	}
	}

	def ScalaRecursiveNonThreadedVaryingSize(minSize: Int, maxSize: Int, repititions: Int, isotropic: Boolean) {
	
	var data = scala.collection.mutable.MutableList[String]();
	data += PCInfo();
	data += "ScalaRecursiveNonThreadedVaryingSize";
	data += "Grid" + '\t' + "Median" + '\t' + "Min" + '\t' + "Max" + '\t' + "StdDev";

	//Warmup
	for(i <- 0 to 9999 by 1) {
		CombiGrid.CombiGrid(buildLevelVector(15, 3, false)); //Create
		CombiGrid.FillCurrentGridWithOnes();                 //Fill
		CombiGrid.hierarchizeRecursive();                  //Run
	}
	System.gc();

	for(dim <- 2 to 5 by 1) {
		for(size <- minSize to  maxSize by 1) {
			var minTime=Long.MaxValue;
			var maxTime=Long.MinValue;  
			var times: Array[Long] = new Array[Long](repititions);
			var totalTime: Long = 0;

			for(i <- 0 to repititions-1 by 1) {
				System.gc();
				try {Thread.sleep(1000);} catch {
				case e: InterruptedException => println("error: " + e)
				}
				CombiGrid.CombiGrid(buildLevelVector(size, dim, isotropic)); //Create
				CombiGrid.FillCurrentGridWithOnes();                 //Fill
				var start = System.currentTimeMillis();
				CombiGrid.hierarchizeRecursive();                  //Run
				val end = System.currentTimeMillis();
				val time = end - start;
				if(time < minTime)
					minTime = time;
				if(time > maxTime)
					maxTime = time;
				times(i) = time;
				totalTime += time;
			}
			val avgTime: Long = totalTime / repititions;
      
      var median: Long=0;
        scala.util.Sorting.quickSort(times); //Sorts from smallest to largest
        if (repititions % 2 == 0) { //if even, take average over two middle values as median.
          median = (times(repititions/2) + times((repititions/2)+1))/2; //take average of the two median values, if even.
        } else { median=times(repititions/2)}
      
			//Calculate standard deviation
			var dev: Double = 0;
      for(d <- times) {
				dev += Math.pow(d - avgTime, 2);
			}
			dev = Math.sqrt(dev / (times.length - 1));
			data += ("" + CombiGrid.getLevelVector() + '\t' + median + '\t' + minTime + '\t' + maxTime + '\t' + dev);
		}    
	}
	if(isotropic) {
		val form = "Isotropic";
		writeToFile("Scala RecursiveNonThreaded VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	} else {
		val form = "Anisotropic";
		writeToFile("Scala RecursiveNonThreaded VaryingSize " + minSize + " " + maxSize + " " + form, data.toList);
	}
	}



	def writeToFile(filename: String, data: List[String])
	{
		try {
			val dir: File = new File(testName);
		  if(!dir.exists())
			dir.mkdir();
		  val file: File = new File(testName + File.separator + filename + ".dat");

	  	// if file exists deletes it then creates it again
		  if (file.exists()) {
			  file.delete();
		  }
		  file.createNewFile();
		  val fw: FileWriter = new FileWriter(file.getAbsoluteFile());
		  val bw: BufferedWriter = new BufferedWriter(fw);
		  for(s <- data) {
        println(s);
  			bw.write(s);
	  		bw.write('\n');
  		}
	  	bw.close();
      fw.close();
	    } catch {
		    case e: IOException => println("IOError: " + e.printStackTrace())
		  }
	}
}

