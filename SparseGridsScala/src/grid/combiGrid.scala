package grid



import scala.concurrent._
import scala.collection.Parallel
import scala.concurrent.ExecutionContext.Implicits.global
import scala.collection.parallel.mutable.ParSeq
import scala.collection.parallel.mutable.ParArray


object CombiGrid {
	//TODO thread the unidirectional as in parallelstream (divide each poleblock into a thread.)
	//global variables:
	var levels:Array[Int] = Array(0);
  var dimensions:Int = 0;
  var grid:ParArray[Double] = ParArray(0);
  var gridSize:Int = 0;
  var pointsPerDimension:Array[Int] = Array(0);
  var strides:Array[Int]=Array(0);


//Variables for the recursion:
var recTile=5; //Some values are hard coded in the original.
val recTallPar=0.3; 
var recMaxSpawn=9; //Public, for varying within the test.
var recMinSpawn=6; //Public, for varying within the test.
var rf:Array[Float]=Array(0);


//Main. Use to launch.
def main(args: Array[String]): Unit = {
    println("runs")
		levels = Array(8, 7)
    CombiGrid(levels);    
		FillArrayWithOnes(grid);
    
//    PrintArray(grid)
    
    hierarchizeRecursive();    
	   //hierarchizeRecursiveThreaded();
     //hierarchizeRecursiveThreadedWhileLoops()
 
//     PrintArray(grid)
     
	var grid2 = new ParArray[Double](grid.size)
	for(i <- 0 to grid.size - 1) {grid2(i) = grid(i)}
     
  CombiGrid(levels);    
  FillArrayWithOnes(grid);
  
  
    hierarchizePoleBlockParArray() //verified
//  hierarchizeUnoptimized() //Verified
//  hierarchizeUnoptimizedParArray() //Verified
//  hierarchizeUnoptimizedSequential() //Verified
  
  
   PrintArray(grid);

  if(compare(grid2))
    System.out.println("Grids are equal")
  else
    System.out.println("Grid are not equal!")
            
}
//The following aligns the grid. Only works with 
def transformAllOnesToAligned(alignment: Int) {
  val appendageLength: Int = (grid.length/alignment)
  var appendage = new ParArray[Double](appendageLength)
  
  FillArrayWithOnes(appendage)
  
  grid = grid ++ appendage
  
  for (i <- alignment-1 to grid.length by alignment){
    grid(i) = Integer.MIN_VALUE
  }  
}

def hierarchizeUnoptimizedTest() {

  var dimension:Int=0;
var start:Int=0;
var stride:Int = 1;
var pointsInDimension:Int=0;
var numberOfPoles:Int=0;
var jump:Int=0;


//dimension 1 separate as start of each pole is easier to calculate
pointsInDimension = pointsPerDimension(0);
numberOfPoles = gridSize / pointsInDimension;
for (i <- 0 to numberOfPoles-1){
  start = i * pointsInDimension;
  hierarchize1DUnoptimized(start, 1, pointsInDimension, 0);
}
// end dimension 1

for(dimension <- 1 to dimensions-1){ // hierarchize all dimensions
  stride *= pointsInDimension;
  pointsInDimension = pointsPerDimension(dimension);
  jump = stride * pointsInDimension;
  numberOfPoles = gridSize / pointsInDimension;
  val jumps = numberOfPoles / stride;
  //System.out.println("Dimension: " + dimension + '\t' + "Jumps: " + jumps);
  for (i <- 0 to numberOfPoles-1){ // integer operations form bottleneck here -- nested loops are twice as slow
    val div = i / stride;
    start = div * jump + (i % stride);
    hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
  }
} // end loop over dimension 2 to d
}

def getLevelVector(): String = {
  val s: StringBuilder  = new StringBuilder();
    var sum = 0;
    for(i <- levels) {
      sum += i;
      s.append("" +  i + " ");
    }
    s.append("[" + sum + "]");
    s.toString();
}

//Small helper functions
def myPow2(j: Int): Int = { //Doubles input
		return 1<<j;
}


def compare(grid2:ParArray[Double]): Boolean = { //Compares two grids
		return grid.equals(grid2)
}



//Create a Grid, based on levels, import scala.collection.Paralleland sets all global variables, based on the level-array. Will return with unset values.
def CombiGrid(levelsInput: Array[Int]){
  levels = levelsInput.clone();
  
  recTile = levels(0); //Important for the recursive algorithms speed.
  
  var maxlvl=0;
  var maxval=0;

  recMaxSpawn = levels(0)-2
  recMinSpawn = levels(0)-5

  
	dimensions = levelsInput.length;
	pointsPerDimension = new Array[Int](dimensions);
	pointsPerDimension(0) = myPow2(levels(0)) - 1;
	gridSize = pointsPerDimension(0);  
	strides = new Array[Int](dimensions+1);
	strides(0)=1;
	var size = 1;
	for (i <- 1 to dimensions-1) {
		pointsPerDimension(i) = myPow2(levels(i)) - 1;
		strides(i) = gridSize;
		gridSize *= pointsPerDimension(i);
	}
	strides(dimensions)=gridSize;
	grid = new ParArray[Double](gridSize);

	rf =  new Array[Float](dimensions);
	for(j <- 0 to dimensions-1 by 1){
		rf(j)=1.0f;
	}
  
}

def FillCurrentGridWithOnes(){
  FillArrayWithOnes(grid);
}


def FillArrayWithOnes(grid: ParArray[Double]){
	for (i <- 0 to grid.length-1){
		grid(i)=1;
	}
}


def printValues2DArr(size:Int, offset:Int, n0:Int){
	for (i <- 1 to size) {
		System.out.print(grid(offset + i - 1))
		System.out.print('\t')
		if (i % n0 == 0)
			System.out.print('\n')
	}
	System.out.print('\n')
}



def PrintArray(grid: ParArray[Double]){
	System.out.println()
	if (dimensions <= 2) printValues2DArr(gridSize, 0, pointsPerDimension(0))
	else {
		var currentLevels = new Array[Int](dimensions)
				val chunkSize = pointsPerDimension(0) * pointsPerDimension(1)
				val numberOfChunks = gridSize / chunkSize
				for (ctr <- 0 to numberOfChunks-1){
					printValues2DArr(chunkSize, ctr*chunkSize, pointsPerDimension(0))
					currentLevels(2)+=1
					for (dd <- 2 to dimensions-1){
						if (currentLevels(dd) == pointsPerDimension(dd)) {
							currentLevels(dd) = 0
									if(dd + 1 < dimensions)
										currentLevels(dd + 1)+=1
						}
					}
				}
	}
}  

def hierarchizeParArrayUnoptimized() {

	var dimension:Int=0;
  var start:Int=0;
  var stride:Int = 1;
  var pointsInDimension:Int=0;
  var numberOfPoles:Int=0;
  var jump:Int=0;


  //dimension 1 separate as start of each pole is easier to calculate
  pointsInDimension = pointsPerDimension(0);
  numberOfPoles = gridSize / pointsInDimension;
  for (i <- 0 to numberOfPoles-1){
	  start = i * pointsInDimension;
	  hierarchize1DUnoptimized(start, 1, pointsInDimension, 0);
  }
  // end dimension 1

  for(dimension <- 1 to dimensions-1){ // hierarchize all dimensions
	  stride *= pointsInDimension;
	  pointsInDimension = pointsPerDimension(dimension);
	  jump = stride * pointsInDimension;
	  numberOfPoles = gridSize / pointsInDimension;
	  val jumps = numberOfPoles / stride;
  	//System.out.println("Dimension: " + dimension + '\t' + "Jumps: " + jumps);
	  for (i <- 0 to numberOfPoles-1){ // integer operations form bottleneck here -- nested loops are twice as slow
  		val div = i / stride;
	  	start = div * jump + (i % stride);
		  hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
	  }
  } // end loop over dimension 2 to d
}

def hierarchize1DUnoptimized(start:Int, stride:Int, size:Int, dimension:Int) { //Check if Scala has recommended ways of parallelisation.

	var level:Int = levels(dimension);
  var steps = myPow2(level - 1);
  var offset = 0;
  var stepSize = 2;
  var parentOffset = 1;
  var left:Double = 0;
  var right:Double = 0;

  
  
  for(l <- level-1 to 2 by -1) {
	  var parOffsetStrided = parentOffset*stride;
	  grid(start + offset * stride) -= 0.5 * grid(start + offset * stride + parOffsetStrided);
	  offset += stepSize;
	  left = 0.5 * grid(start + offset * stride - parOffsetStrided);
	  for (ctr <- 1 to steps-2) {
		  var val1 = grid(start + offset * stride);
		  right = 0.5 * grid(start + offset * stride + parOffsetStrided);
		  var val2 = val1 - left;
		  var val3 = val2 - right;
		  grid(start+offset*stride) = val3;
		  left = right;
		  offset += stepSize;
	  } 

	  grid(start + offset * stride) -= right;
  	//steps = steps >> 1;
	  steps = steps / 2;
	  offset = myPow2(levels(dimension) - l) - 1;
	  parentOffset =  stepSize;
  	//stepSize = stepSize << 1;
	  stepSize = stepSize * 2;
  }

  right = 0.5 * grid(start + (offset + parentOffset) * stride);
  grid(start + offset * stride) -= right;
  offset += stepSize;
  grid(start + offset * stride) -= right;
}

def hierarchizeRecursive(){
	// this method starts the recursion, using the hierarchizeRec-call.

	var centerInd:Array[Int] = new Array[Int](dimensions);
  var fullInterval:Content = new Content();
  for(i <- 0 to dimensions-1 by 1) {
	  fullInterval.l(i) = levels(i) - 1; // boundary need not be split away
	  centerInd(i) = myPow2(levels(i) - 1) - 1;
  }

  fullInterval.l(6) = 0 ; // no predecessors to the left
  fullInterval.l(7) = 0 ; // no predecessors to the right
  var center = pos(centerInd);
  hierarchizeRec(0, dimensions, center, fullInterval, 1);
}

def hierarchizeRecursiveThreaded(){
	// this method starts the recursion, using the hierarchizeRec-call.

	var centerInd:Array[Int] = new Array[Int](dimensions);
  var fullInterval:Content = new Content();
  for(i <- 0 to dimensions-1 by 1) {
	  fullInterval.l(i) = levels(i) - 1; // boundary need not be split away
	  centerInd(i) = myPow2(levels(i) - 1) - 1;
  }

  fullInterval.l(6) = 0 ; // no predecessors to the left
  fullInterval.l(7) = 0 ; // no predecessors to the right
  var center = pos(centerInd);

  hierRecThreads(0, dimensions, center, fullInterval, 1);
}


def hierRecThreads(si:Int, ti:Int, centeri:Int, interval:Content, leveli:Int){
		val level = leveli;
		var s=si;
		val ic = new Content();
		ic.copy(interval.l)
		val t=ti;
		val center= centeri;
		var localSize = 0; //Set's sum of levels, to determine size.
		for (i <- 1 to dimensions-1) {
			if(ic.l(i) > 0) {
				localSize += ic.l(i);
			}
		}
    
		if(localSize == 0) { // singleton cache line
			if(ic.l(0) <= 0) { // real singletons, t, center, inputContent
				println("s = " + s + ", t = " + t)
				for(i <- s to t-1 by 1) { 
					var rmask = myPow2(i);
					var dist = myPow2(-ic.l(i));
					var lVal:Double=0; //don't define before using.
					var rVal:Double=0; //doubles - define at first use instead.
					if((ic.l(6) & rmask)!=0) { //Checks if the bitwise combination equals to 1.
						var posLeft = center - dist*strides(i);
						var lVal = grid(posLeft);
					}
					else {
						var posLeft = -1;
						var lVal = 0.0;
					}
					if((ic.l(7) & rmask) !=0) { //centerInd[i] + dist < n[i] ) { // it should be == n[i], but hey
						var posRight = center + dist*strides(i);
						var rVal = grid(posRight);
					}
					else {
						var posRight = -1;
						var rVal = 0.0;
					}
					grid(center) = stencil(grid(center), lVal, rVal, true, true);
				}
			} 
			else {
				if( s == 0 ) { // actually hierarchize in dir 0
					var rmask = (1 << 0); // replace by iterative?
					var dist = myPow2(ic.l(0));
					var leftBdVal:Double=0.0;
					var rightBdVal:Double=0.0;
					var leftBdPos:Int = 0;
					var rightBdPos:Int =0;
					// if we don't split in dim0, we know we are at the boundary...
					if((ic.l(6) & rmask) !=0) {
						leftBdPos = center - dist;
						leftBdVal = grid(leftBdPos);
					}
					else {
						leftBdPos = -1;
						leftBdVal = 0.0;
					}
					if((ic.l(7) & rmask) !=0) {
						rightBdPos = center + dist;
						rightBdVal = grid(rightBdPos);
					}
					else {
						rightBdPos = -1;
						rightBdVal = 0.0;
					}
					var step = 1;
					while(step < dist) {
						var start = center - dist + step;
						grid(start) = stencil(grid(start), leftBdVal, grid(start + step), true, true);
						start += 2*step;
						while( start < center + dist - step ) {
							grid(start) = stencil(grid(start),grid(start-step),grid(start+step), true, true);
							start += 2*step;
						}
						assert( start == center+dist-step );
						grid(start) = stencil(grid(start), grid(start-step), rightBdVal, true, true);
						step *= 2;
					}
					// while of levels
					grid(center) = stencil(grid(center), leftBdVal, rightBdVal, true, true);
					s = 1; // hierarchized in dim 0
				}
				var d0dist = myPow2(ic.l(0));
				var first = - d0dist +1;
				var last = + d0dist -1;
				for(dim <- s to t-1 by 1) {
					var rmask = (1 << dim); 
					var dist = myPow2(-ic.l(dim));
					if(((ic.l(6) & rmask)) !=0 && (ic.l(7) & rmask)!=0) {
						for(i <- first to last by 1) {
							var offset=dist*strides(dim);
							grid(center+i) = stencil(grid(center+i), grid(center-offset+i), grid(center+offset+i), true, true); 
						}
					}
					if( (ic.l(6) & rmask)!=0 && (ic.l(7) & rmask)==0 ) {
						for(i <- first to last by 1) {
							var offset=dist*strides(dim);
							grid(center+i) = stencil(grid(center+i), grid(center-offset+i), 0, true, false);
						}
					}
					if( (ic.l(6) & rmask)==0 && (ic.l(7) & rmask)!=0 ) {
						for(i <- first to last by 1) {
							var offset=dist*strides(dim);
							grid(center+i) = stencil(grid(center+i), 0, grid(center+offset+i), false, true);
						}
					} 
				}
			}
		}
		else {
			var r=0;

			// ic used for right
			var maxl = 0
					for(i <- 1 to dimensions-1) {
						if(ic.l(i) > maxl) {
							r = i
									maxl = ic.l(i)
						}
					}
      
			var midI:Content= new Content();
					midI.copy(ic.l);
					midI.l(r) = -midI.l(r);
					ic.l(r) = ic.l(r) - 1;
					var dist = myPow2(ic.l(r)); // already reduced!
					
					var leftI:Content=new Content();
					leftI.copy(ic.l);
					var rmask = myPow2(r);
					ic.l(6) |= rmask;
					leftI.l(7) |= rmask;
					dist *= strides(r); // already reduced!
					if(r < s) r = s; // avoid calls!
					if(r > t) r = t;

          val sOut = s;
          val tOut = t;
          val rOut = r;
          val cOut = center;
          
					if ((localSize >= recMinSpawn) && (localSize <= recMaxSpawn) && r!=0){
          
          if(r > s) {

                hierRecThreads(sOut, rOut, cOut, midI, level + 1);

						}
            
          val t2: Thread = new Thread { override def run(): Unit = {
            hierRecThreads(sOut, tOut, cOut - dist, leftI, level + 1);
          }};
          t2.start();
       
        
        val t3: Thread = new Thread { override def run(): Unit = {
          hierRecThreads(sOut, tOut, cOut + dist, ic, level + 1)
        }};
        t3.start();
              
        t2.join()
        t3.join()
        
        if(t > r) {
          hierRecThreads(rOut, tOut, cOut, midI, level + 1);
        }
        
			} else { //don't run in seperate threads.
				if(r > s) hierRecThreads(sOut, rOut, cOut, midI, 0);
        hierRecThreads(sOut, tOut, cOut - dist, leftI, 0);
				hierRecThreads(sOut, tOut, cOut + dist, ic, 0);
				if(t > r) hierRecThreads(rOut, tOut, cOut, midI, 0);
			}
   }
}
  

def hierarchizeRecursiveThreadedWhileLoops(){
  // this method starts the recursion, using the hierarchizeRec-call.

  var centerInd:Array[Int] = new Array[Int](dimensions);
  var fullInterval:Content = new Content();
  var i = 0;
  while (i < dimensions) {
    fullInterval.l(i) = levels(i) - 1; // boundary need not be split away
    centerInd(i) = myPow2(levels(i) - 1) - 1;
    i = i + 1
  }

  fullInterval.l(6) = 0 ; // no predecessors to the left
  fullInterval.l(7) = 0 ; // no predecessors to the right
  var center = pos(centerInd);
  
  hierRecThreadsWhileLoops(0, dimensions, center, fullInterval, 1);
}

   
def hierRecThreadsWhileLoops(si:Int, ti:Int, centeri:Int, interval:Content, leveli:Int){
  var level = leveli;
  var s=si;
  var ic = new Content();
  ic.copy(interval.l)
  var t=ti;
  var center= centeri;
  var localSize = 0; //Set's sum of levels, to determine size.
  var localsizefiller = 1;
  while (localsizefiller < dimensions) {
    if(ic.l(localsizefiller) > 0) {
      localSize += ic.l(localsizefiller);
    }
    localsizefiller = localsizefiller + 1
  }
  
  if(localSize == 0) { // singleton cache line
    if(ic.l(0) <= 0) { // real singletons, t, center, inputContent
      println("s = " + s + ", t = " + t)
      var i=s
      while(i < t) { // Counts from s to to 
        var rmask = myPow2(i)
        var dist = myPow2(-ic.l(i))
        var lVal:Double=0; 
        var rVal:Double=0; 
        if((ic.l(6) & rmask)!=0) { //Checks if the bitwise combination equals to 1.
          var posLeft = center - dist*strides(i);
          var lVal = grid(posLeft);
        }
        else {
          var posLeft = -1;
          var lVal = 0.0;
        }
        if((ic.l(7) & rmask) !=0) {
          var posRight = center + dist*strides(i);
          var rVal = grid(posRight);
        }
        else {
          var posRight = -1;
          var rVal = 0.0;
        }
        grid(center) = stencil(grid(center), lVal, rVal, true, true);
        i = i + 1
      }
    } 
    else {
      if( s == 0 ) { 
        var rmask = (1 << 0);
        var dist = myPow2(ic.l(0));
        var leftBdVal:Double=0.0;
        var rightBdVal:Double=0.0;
        var leftBdPos:Int = 0;
        var rightBdPos:Int =0;
        // if we don't split in dim0, we know we are at the boundary...
        if((ic.l(6) & rmask) !=0) {
          leftBdPos = center - dist;
          leftBdVal = grid(leftBdPos);
        }
        else {
          leftBdPos = -1;
          leftBdVal = 0.0;
        }
        if((ic.l(7) & rmask) !=0) {
          rightBdPos = center + dist;
          rightBdVal = grid(rightBdPos);
        }
        else {
          rightBdPos = -1;
          rightBdVal = 0.0;
        }
        var step = 1;
        while(step < dist) {
          var start = center - dist + step;
          grid(start) = stencil(grid(start), leftBdVal, grid(start + step), true, true);
          start += 2*step;
          while( start < center + dist - step ) {
            grid(start) = stencil(grid(start),grid(start-step),grid(start+step), true, true);
            start += 2*step;
          }
          grid(start) = stencil(grid(start), grid(start-step), rightBdVal, true, true);
          step *= 2;
        }
        // while of levels
        grid(center) = stencil(grid(center), leftBdVal, rightBdVal, true, true);
        s = 1; // hierarchized in dim 0
      }
      var d0dist = myPow2(ic.l(0));
      var first = - d0dist +1;
      var last = + d0dist -1;
      
      var dim  = s
      while(dim < t) {
        
        var rmask = (1 << dim); 
        var dist = myPow2(-ic.l(dim));
        
        if(((ic.l(6) & rmask)) !=0 && (ic.l(7) & rmask)!=0) {
          var ia = first
          while(ia <= last) {
            var offset=dist*strides(dim);
            grid(center+ia) = stencil(grid(center+ia), grid(center-offset+ia), grid(center+offset+ia), true, true);
            ia = ia + 1
          }
        }
        
        if( (ic.l(6) & rmask)!=0 && (ic.l(7) & rmask)==0 ) {
          var ib = first
          while(ib <= last) {
            var offset=dist*strides(dim);
            grid(center+ib) = stencil(grid(center+ib), grid(center-offset+ib), 0, true, false);
            ib = ib +1
          }
        }
        
        if( (ic.l(6) & rmask)==0 && (ic.l(7) & rmask)!=0 ) {
          var id = first
          while(id <= last) {
            var offset=dist*strides(dim);
            grid(center+id) = stencil(grid(center+id), 0, grid(center+offset+id), false, true);
            id = id + 1
          }
        } 
        
        dim = dim + 1
      }
    }
  }
  else {
    var r=0;
    var maxl = 0
    var dimcount = 1
    while (dimcount < dimensions) {
      if(ic.l(dimcount) > maxl) {
        r = dimcount
        maxl = ic.l(dimcount)
      }
      dimcount = dimcount + 1
    }
    
    var midI:Content= new Content();
    midI.copy(ic.l);
    midI.l(r) = -midI.l(r);
    ic.l(r) = ic.l(r) - 1;
    var dist = myPow2(ic.l(r)); // already reduced!
    var leftI:Content=new Content();
    leftI.copy(ic.l);
    var rmask = myPow2(r);
    ic.l(6) |= rmask;
    leftI.l(7) |= rmask;
    dist *= strides(r); // already reduced!
    if(r < s) r = s; // avoid calls!
    if(r > t) r = t;

    val sOut = s;
    val tOut = t;
    val rOut = r;
    val cOut = center;

    if ((localSize >= recMinSpawn) && (localSize <= recMaxSpawn) && r!=0){

      if(r > s) {

         hierRecThreadsWhileLoops(sOut, rOut, cOut, midI, level + 1);

      }

      val t2: Thread = new Thread { override def run(): Unit = {
          hierRecThreadsWhileLoops(sOut, tOut, cOut - dist, leftI, level + 1);
      }};
      t2.start();


      val t3: Thread = new Thread { override def run(): Unit = {
          hierRecThreadsWhileLoops(sOut, tOut, cOut + dist, ic, level + 1)
      }};
      t3.start();

      t2.join()
      t3.join()
          
      if(t > r) hierRecThreadsWhileLoops(rOut, tOut, cOut, midI, level+1);
     

    } else { //don't run in seperate threads.
      if(r > s) hierRecThreadsWhileLoops(sOut, rOut, cOut, midI, 0);
      hierRecThreadsWhileLoops(sOut, tOut, cOut - dist, leftI, 0);
      hierRecThreadsWhileLoops(sOut, tOut, cOut + dist, ic, 0);
      if(t > r) hierRecThreadsWhileLoops(rOut, tOut, cOut, midI, 0);
    }
  }
}



def hierarchizeRec(si:Int, t:Int, center:Int, interval:Content, level:Int) {//Actual recursion method. Not threaded.
	var s = si;
	var ic = new Content();
	ic.copy(interval.l)
	var localSize = 0; //Set's sum of levels, to determine size.
	for (i <- 1 to dimensions-1) {
		if(ic.l(i) > 0) {
			localSize += ic.l(i);
		}
	}

	if(localSize == 0) { // singleton cache line
		if(ic.l(0) <= 0) { // real singletons, t, center, inputContent
			for(i <- s to t-1 by 1) {
				var rmask = myPow2(i);
				var dist = myPow2(-ic.l(i));
				//System.out.println(dist);
				var lVal:Double=0; //don't define before using.
				var rVal:Double=0; //doubles - define at first use instead.
				if((ic.l(6) & rmask)!=0) { //Checks if the bitwise combination equals to 1.
					var posLeft = center - dist*strides(i);
					var lVal = grid(posLeft);
				}
				else {
					var posLeft = -1;
					var lVal = 0.0;
				}
				if((ic.l(7) & rmask) !=0) { //centerInd[i] + dist < n[i] ) { // it should be == n[i], but hey
					var posRight = center + dist*strides(i);
					var rVal = grid(posRight);
				}
				else {
					var posRight = -1;
					var rVal = 0.0;
				}
				grid(center) = stencil(grid(center), lVal, rVal, true, true);
			}
		} 
		else {
			if( s == 0 ) { // actually hierarchize in dir 0
				var rmask = (1 << 0); // replace by iterative?
				var dist = myPow2(ic.l(0));
				var leftBdVal:Double=0.0;
				var rightBdVal:Double=0.0;
				var leftBdPos:Int = 0;
				var rightBdPos:Int =0;
				// if we don't split in dim0, we know we are at the boundary...
				if((ic.l(6) & rmask) !=0) { //centerInd[i] - dist >= 0 ) { // it should be == -1, but hey
					leftBdPos = center - dist;
					leftBdVal = grid(leftBdPos);
				}
				else {
					leftBdPos = -1;
					leftBdVal = 0.0;
				}
				if((ic.l(7) & rmask) !=0) { //centerInd[i] + dist < n[i] ) { // it should be == n[i], but hey
					rightBdPos = center + dist;
					rightBdVal = grid(rightBdPos);
				}
				else {
					rightBdPos = -1;
					rightBdVal = 0.0;
				}
				var step = 1;
				while(step < dist) {
					var start = center - dist + step;
					//System.out.println("Center: " + center + '\t' + "Dist: " + dist + '\t' + "Step: " + step + '\t' + "Start: " + start);
					grid(start) = stencil(grid(start), leftBdVal, grid(start + step), true, true);
					start += 2*step;
					while( start < center + dist - step ) {
						grid(start) = stencil(grid(start),grid(start-step),grid(start+step), true, true);
						start += 2*step;
					}
					assert( start == center+dist-step );
					grid(start) = stencil(grid(start), grid(start-step), rightBdVal, true, true);
					step *= 2;
				}
				// while of levels
				grid(center) = stencil(grid(center), leftBdVal, rightBdVal, true, true);
				s = 1; // hierarchized in dim 0
			}
			var d0dist = myPow2(ic.l(0));
			var first = - d0dist +1;
			var last = + d0dist -1;
			for(dim <- s to t-1 by 1) {
				var rmask = (1 << dim); // replace by iterative?
				var dist = myPow2(-ic.l(dim));
				//assert(0== (center+first) %4 );
				//assert(2== (center+last) %4 );
				if(((ic.l(6) & rmask)) !=0 && (ic.l(7) & rmask)!=0) {
					for(i <- first to last by 1) {
						var offset=dist*strides(dim);
						//System.out.println("Center: " + center + '\t' + "Offset: " + offset + '\t' + "i: " + i + '\t' + "Left " + '\t' + "Right");
						//System.out.println("Grid before: " + grid(center+i));
						grid(center+i) = stencil(grid(center+i), grid(center-offset+i), grid(center+offset+i), true, true); 
						//System.out.println("Grid after: " + grid(center+i));
					}
				}
				if( (ic.l(6) & rmask)!=0 && (ic.l(7) & rmask)==0 ) {
					for(i <- first to last by 1) {
						var offset=dist*strides(dim);
						//System.out.println("Center: " + center + '\t' + "Offset: " + offset + '\t' + "i: " + i + '\t' + "Left");
						//System.out.println("Grid before: " + grid(center+i));
						grid(center+i) = stencil(grid(center+i), grid(center-offset+i), 0, true, false);
						//System.out.println("Grid after: " + grid(center+i));
					}
				}
				if( (ic.l(6) & rmask)==0 && (ic.l(7) & rmask)!=0 ) {
					for(i <- first to last by 1) {
						var offset=dist*strides(dim);
						//System.out.println("Center: " + center + '\t' + "Offset: " + offset + '\t' + "i: " + i + '\t' + "Right");
						//System.out.println("Grid before: " + grid(center+i));
						//System.out.println("Grid pos " + (center+offset+i) + ": " + grid(center+offset+i));
						grid(center+i) = stencil(grid(center+i), 0, grid(center+offset+i), false, true);
						//System.out.println("Grid after: " + grid(center+i));
					}
				} 
			}
		}
	}
	else {
		var r=0;

		var maxl = 0
				for(i <- 1 to dimensions-1) {
					if(ic.l(i) > maxl) {
						r = i
								maxl = ic.l(i)
					}
				}

		var midI:Content= new Content();
		midI.copy(ic.l);
		midI.l(r) = -midI.l(r);
		ic.l(r) = ic.l(r) - 1;
		var dist = myPow2(ic.l(r)); // already reduced!

		var leftI:Content=new Content();
		leftI.copy(ic.l);
		var rmask = myPow2(r);
		ic.l(6) |= rmask;
		leftI.l(7) |= rmask;
		dist *= strides(r); // already reduced!
		if(r < s) r = s; // avoid calls!
		if(r > t) r = t;

		if (r!=0){
   		if(r > s) {
	    		hierarchizeRec(s, r, center, midI, level + 1);
			}

		hierarchizeRec(s, t, center - dist, leftI, level + 1);
		hierarchizeRec(s, t, center + dist, ic, level + 1);

  	if (t > r) {
			hierarchizeRec(r, t, center, midI, level + 1);
		}
  }
 }
}

//HELPER METHODS FOR RECURSION
//this subtracts half of right and left from center.
def stencil( center:Double, left:Double, right:Double, goLeft:Boolean, goRight:Boolean): Double = { 
		if (!goLeft) {
			return (center -.5* right);
		} else if (!goRight) {
			return (center -.5* left);
		}
		return (center -.5*left -.5* right);
}

def pos( index:Array[Int]): Int = {
		var retPos = 0;
		for(i <- 0 to dimensions-1 by 1) {
			//assert( index[i] < pointsPerDimension[i] );
			//assert( 0 <= index[i] );
			retPos += index(i)*strides(i);
		}
		return retPos;
};

class Content { //object for holding both integer and byte-array.
	var l:Array[Int]=Array[Int](0,0,0,0,0,0,0,0);

  def copy(lInput:Array[Int]){
	  System.arraycopy(lInput, 0, l, 0, lInput.length);
  }
  
  def printInterval() {
	  System.out.print("" + l(0));
	  for(i <- 1 to 7)
		  System.out.print(", " + l(i))
		  System.out.println()
  }
}

/***
   * Hierarchizes the grid using a parallel unoptimized algorithm using the Java Parallel Stream framework.
   * This method creates a pole object for each pole, and then uses the parallelStream() method to run all the poles
   * in parallel. All the parallelism is done by the Java framework.
   */
    def hierarchizeUnoptimizedParArray() {
      var dimension:Int=0;
      var start:Int=1;
      var stride = 1;
      var polesBuilder: scala.collection.mutable.ListBuffer[Pole]=scala.collection.mutable.ListBuffer[Pole]();// Must be parArray, as we expect this will start the operations in parallel.
 
  
      //dimension 1 separate as start of each pole is easier to calculate
      var pointsInDimension = pointsPerDimension(0);
      var numberOfPoles = gridSize / pointsInDimension;
      for (i <- 0 to numberOfPoles-1 by 1){
        start = i * pointsInDimension;
        var p: Pole= new Pole; 
        p.Pole(start, 1, pointsInDimension, 0);
        polesBuilder += p;
      }
      // end dimension 1
      
      var poles = polesBuilder.toParArray
      
      poles //Poles is parArray, so this should happen in parallel.
        .foreach { pole => hierarchize1DUnoptimized(pole.start, pole.stride, pole.pointsInDimension, pole.dimension) };
  
      for(dimension <- 1 to dimensions-1 by 1){ // hierarchize all dimensions
        polesBuilder = scala.collection.mutable.ListBuffer[Pole]();
        
        stride *= pointsInDimension;
        pointsInDimension = pointsPerDimension(dimension);
        var jump = stride * pointsInDimension;
        numberOfPoles = gridSize / pointsInDimension;
        for (i <- 0 to numberOfPoles-1 by 1){ // integer operations form bottleneck here -- nested loops are twice as slow
          val div = i / stride;
          start = div * jump + (i % stride);
          var p: Pole = new Pole();
          p.Pole(start, stride, pointsInDimension, dimension);
          polesBuilder += p;
        }
        
        poles = polesBuilder.toParArray
        poles
        .foreach{pole => hierarchize1DUnoptimized(pole.start, pole.stride, pole.pointsInDimension, pole.dimension)};
   } // end loop over dimension 2 to d
}
    
def hierarchizeUnoptimized() {

  var dimension:Int=0;
  var start:Int=0;
  var stride:Int = 1;
  var pointsInDimension:Int=0;
  var numberOfPoles:Int=0;
  var jump:Int=0;


//dimension 1 separate as start of each pole is easier to calculate
  pointsInDimension = pointsPerDimension(0);
  numberOfPoles = gridSize / pointsInDimension;
  for (i <- 0 to numberOfPoles-1){
    start = i * pointsInDimension;
    hierarchize1DUnoptimized(start, 1, pointsInDimension, 0);
  }
// end dimension 1

  for(dimension <- 1 to dimensions-1){ // hierarchize all dimensions
    stride *= pointsInDimension;
    pointsInDimension = pointsPerDimension(dimension);
    jump = stride * pointsInDimension;
    numberOfPoles = gridSize / pointsInDimension;
    val jumps = numberOfPoles / stride;
    //System.out.println("Dimension: " + dimension + '\t' + "Jumps: " + jumps);
    for (i <- 0 to numberOfPoles-1){ // integer operations form bottleneck here -- nested loops are twice as slow
      val div = i / stride;
      start = div * jump + (i % stride);
      hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
    }
  } // end loop over dimension 2 to d
}
    
def hierarchizeUnoptimizedSequential() {
    var start:Int=1;
    var stride = 1;
    var poles=scala.collection.mutable.ArrayBuffer[Pole]();// Must be parArray, as we expect this will start the operations in parallel.
  
  
    //dimension 1 separate as start of each pole is easier to calculate
    var pointsInDimension = pointsPerDimension(0);
    var numberOfPoles = gridSize / pointsInDimension;
    
    for (i <- 0 to numberOfPoles-1){
      start = i * pointsInDimension;
      var p: Pole= new Pole; 
      p.Pole(start, 1, pointsInDimension, 0);
      poles += p;      
    }
    // end dimension 1
    var index = 0
    while (index < poles.length) {
       hierarchize1DUnoptimized(poles(index).start, poles(index).stride, poles(index).pointsInDimension, poles(index).dimension) 
       index = index + 1
    }
    var dimension = 1;
//    while (dimension < dimensions) { // hierarchize all dimensions
    for (dimension <- 1 to dimensions-1){
      poles=scala.collection.mutable.ArrayBuffer[Pole](); //this makes a new, empty array.
      stride *= pointsInDimension;
      pointsInDimension = pointsPerDimension(dimension);
      var jump = stride * pointsInDimension;
      numberOfPoles = gridSize / pointsInDimension;
      var i = 0
      while (i < numberOfPoles){ // integer operations form bottleneck here -- nested loops are twice as slow
        val div = i / stride;
        start = div * jump + (i % stride);
        var p: Pole = new Pole();
        p.Pole(start, stride, pointsInDimension, dimension);
        poles += p;
        i = i + 1
      }
      
      var index_2 = 0
      while (index_2 < poles.length) {
       hierarchize1DUnoptimized(poles(index_2).start, poles(index_2).stride, poles(index_2).pointsInDimension, poles(index_2).dimension) 
       index_2 = index_2 + 1
    }
//    dimension = dimension + 1
  } // end loop over dimension 2 to d
}

def hierarchizePoleBlockParArray() {
  var stride = 1;

  //dimension 1 separate as start of each pole is easier to calculate
  var pointsInDimension = pointsPerDimension(0);
  var numberOfPoles = gridSize / pointsInDimension;
  
  for ( i <- 0 to numberOfPoles-1 by 1){
      val start = i * pointsInDimension;
      hierarchize1DUnoptimized(start, 1, pointsInDimension, 0);
  }
  
  for (dimension <-1 to dimensions-1 by 1) { // hierarchize d >=1
      stride *= pointsInDimension;
      pointsInDimension = pointsPerDimension(dimension);
      val jump = stride * pointsInDimension;
      numberOfPoles = gridSize / pointsInDimension;
      var blockSize = 0
      
      if (stride / 16 < 4) {
        blockSize = 4
      } else { blockSize = stride / 16}

      var i = 0;

      while(i < numberOfPoles) {
        val div = i / stride;
        var start = div * jump + (i % stride);
        var j=0;
        while (j < stride) {
          hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
          start = start + 1;
          j = j + 1
        }
        i = i + stride;
      }
    }   // end loop over dimension 2 to d
}

/***
 * Object to contain the variables needed for a pole.
 * Used in the hierarchizeUnoptimizedParallelStream method.
 */

    
class Pole {
  var start: Int=0; //public int start;
  var stride: Int=0;  //public int stride;
  var pointsInDimension: Int=0; //public int pointsInDimension;
  var dimension: Int=0;  //public int dimension;

  def Pole(start: Int, stride:Int, pointsInDimension:Int, dimension:Int) {
    this.start = start;
    this.stride = stride;
    this.pointsInDimension = pointsInDimension;
    this.dimension = dimension;
  }
}

class PoleBlock {
  
  var grid:ParArray[Double]= ParArray[Double](0);
  var stride: Int=0;
  var pointsInDimension: Int=0;
  var dimension: Int=0;
  var jump: Int=0;
  var from: Int=0;
  var to: Int=0;

  def PoleBlock(stride: Int, pointsInDimension: Int, dimension: Int, jump: Int, from: Int, to: Int) {
    this.stride = stride;
    this.pointsInDimension = pointsInDimension;
    this.dimension = dimension;
    this.jump = jump;
    this.from = from;
    this.to = to;
  }

  def hierarchize() {
    for ( i <- from to to-1 by 1) { 
      var div = i / stride;
      var start = div * jump + (i % stride);
      hierarchize1DUnoptimized(start, stride, pointsInDimension, dimension);
    }
  }
}

}




