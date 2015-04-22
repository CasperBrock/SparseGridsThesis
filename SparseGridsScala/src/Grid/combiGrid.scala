package Grid



import scala.concurrent._
import scala.collection.Parallel
import scala.concurrent.ExecutionContext.Implicits.global
import scala.collection.parallel.mutable.ParSeq
import scala.collection.parallel.mutable.ParArray


object CombiGridScala {
	//TODO thread the unidirectional as in parallelstream (divide each poleblock into a thread.)

	//global variables:
	var levels:Array[Int] = Array(0);
var dimensions:Int = 0;
var grid:ParArray[Double] = ParArray(0);
var gridSize:Int = 0;
var pointsPerDimension:Array[Int] = Array(0);
var strides:Array[Int]=Array(0);


//Variables for the recursion:
val recTile=5; //Some values are hard coded in the original.
val recTallPar=0.3; 
val recMaxSpawn=3; //Public, for varying within the test.
val recMinSpawn=2; //Public, for varying within the test.
var rf:Array[Float]=Array(0);


//Main. Use to launch.
def main(args: Array[String]): Unit = {

		levels = Array(2,2,2);
		CombiGrid(levels);
		FillArrayWithOnes(grid);
		//hierarchizeRecursive();
		hierarchizeRecursiveThreaded();
		var grid2 = new ParArray[Double](grid.size)
				for(i <- 0 to grid.size - 1)
					grid2(i) = grid(i)
					//hierarchizeUnoptimized();  
					PrintArray(grid);
						//println(gridSize)
						//pointsPerDimension.foreach(l => print(""+l+'\t'))


						//levels= Array(5, 5, 5, 5, 5);
						CombiGrid(levels);
						FillArrayWithOnes(grid);
						hierarchizeUnoptimized();  
						PrintArray(grid);
						//println(gridSize)
						//pointsPerDimension.foreach(l => print(""+l+'\t'))

						if(compare(grid2))
							System.out.println("Grids are equal")
							else
								System.out.println("Grid are not equal!")

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

def hierarchize1DUnoptimized(start:Int, stride:Int, size:Int, dimension:Int) { //Check if Scala has recommended ways of parallelisation.
	//var level, steps, ctr, offset, parentOffset, stepSize, parOffsetStrided
	//var val1, val2, val3, left, right

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
val H = new hierRecThreads(0, dimensions, center, fullInterval, 1);
H.run()
}

class hierRecThreads(si:Int, ti:Int, centeri:Int, interval:Content, leveli:Int) extends Runnable {

   
	//For running in threads
	 override def run(){
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
					System.out.println(dist);
					var lVal:Double=0; //don't define before using.
					var rVal:Double=0; //doubles - define at first use instead.
					//var posLeft, posRight; //integers
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
					// hierarchize1DUnoptimized(CGIndex start, CGIndex stride, CGIndex size, int dim) does not fit because it never uses boundary
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
				println(s + " " + t)
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
							//            System.out.println("Center: " + center + '\t' + "Offset: " + offset + '\t' + "i: " + i + '\t' + "Right");
							//            System.out.println("Grid before: " + grid(center+i));
							//            System.out.println("Grid pos " + (center+offset+i) + ": " + grid(center+offset+i));
							grid(center+i) = stencil(grid(center+i), 0, grid(center+offset+i), false, true);
							//            System.out.println("Grid after: " + grid(center+i));
						}
					} 
				}
			}
		}
		else {
			var r=0;
			//We added the cast to int in the following line.
			//var maxl=ic.l(0) - recTile - ((recTallPar * localSize).toInt); // block size of pseudo singletons, tall cache assumption. 
			//for(i <- 1 to dimensions-1 by 1){
			//  if( rf(i)*ic.l(i) > maxl ) {
			//    r = i;
			//    maxl = (rf(i)*ic.l(i)).toInt;
			//  }
			//}
			// ic used for right
			var maxl = 0
					for(i <- 1 to dimensions-1) {
						if(ic.l(i) > maxl) {
							r = i
									maxl = ic.l(i)
						}
					}
			//System.out.println();
			//System.out.println("r: " + r);
			//ic.printInterval()
			//System.out.println();
			var midI:Content= new Content();
					midI.copy(ic.l);
					midI.l(r) = -midI.l(r);
					ic.l(r) = ic.l(r) - 1;
					var dist = myPow2(ic.l(r)); // already reduced!
					//leftI.asInt = ic.asInt;
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
							val t1: hierRecThreads = new hierRecThreads(s, r, center, midI, level + 1);
						t1.run();
						}

						val t2: hierRecThreads = new hierRecThreads(s, t, center - dist, leftI, level + 1);
						t2.run();

						val t3: hierRecThreads = new hierRecThreads(s, t, center + dist, ic, level + 1);
						t3.run(); 

						if (t > r) {
							val t4: hierRecThreads = new hierRecThreads(r, t, center, midI, level + 1);
						t4.run();
						} 
					} else { //don't run in seperate threads, on the last level.
						if(r > s) hierarchizeRec(s, r, center, midI, 0);
						hierarchizeRec(s, t, center - dist, leftI, 0);
						hierarchizeRec(s, t, center + dist, ic, 0);
						if(t > r) hierarchizeRec(r, t, center, midI, 0);
					}




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
				System.out.println(dist);
				var lVal:Double=0; //don't define before using.
				var rVal:Double=0; //doubles - define at first use instead.
				//var posLeft, posRight; //integers
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
				// hierarchize1DUnoptimized(CGIndex start, CGIndex stride, CGIndex size, int dim) does not fit because it never uses boundary
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
		//We added the cast to int in the following line.
		//var maxl=ic.l(0) - recTile - ((recTallPar * localSize).toInt); // block size of pseudo singletons, tall cache assumption. 
		//for(i <- 1 to dimensions-1 by 1){
		//	if( rf(i)*ic.l(i) > maxl ) {
		//		r = i;
		//		maxl = (rf(i)*ic.l(i)).toInt;
		//	}
		//}
		// ic used for right
		var maxl = 0
				for(i <- 1 to dimensions-1) {
					if(ic.l(i) > maxl) {
						r = i
								maxl = ic.l(i)
					}
				}
		//System.out.println();
		//System.out.println("r: " + r);
		//ic.printInterval()
		//System.out.println();
		var midI:Content= new Content();
				midI.copy(ic.l);
				midI.l(r) = -midI.l(r);
				ic.l(r) = ic.l(r) - 1;
				var dist = myPow2(ic.l(r)); // already reduced!
				//leftI.asInt = ic.asInt;
				var leftI:Content=new Content();
				leftI.copy(ic.l);
				var rmask = myPow2(r);
				ic.l(6) |= rmask;
				leftI.l(7) |= rmask;
				dist *= strides(r); // already reduced!
				if(r < s) r = s; // avoid calls!
				if(r > t) r = t;

				//if (r!=0){
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
	//variable asInt:Int=0;
	var l:Array[Int]=Array[Int](0,0,0,0,0,0,0,0);

  def copy(lInput:Array[Int]){
	  //asInt = asIntInput;
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
    def hierarchizeUnoptimizedParallelStream() {
      var dimension:Int=0;
      var start:Int=1;
      var stride = 1;
      var poles: ParArray[Pole]=ParArray[Pole]();// Must be parArray, as we expect this will start the operations in parallel.
  
  
      //dimension 1 separate as start of each pole is easier to calculate
      var pointsInDimension = pointsPerDimension(0);
      var numberOfPoles = gridSize / pointsInDimension;
      for (i <- 0 to numberOfPoles-1 by 1){
        start = i * pointsInDimension;
        var p: Pole= new Pole; 
        p.Pole(start, 1, pointsInDimension, 0);
        poles:+ p;
      }
      // end dimension 1
      poles //Poles is parArray, so this should happen in parallel.
      .foreach { pole => hierarchize1DUnoptimized(pole.start, pole.stride, pole.pointsInDimension, pole.dimension) };
  
      for(dimension <- 1 to dimensions-1 by 1){ // hierarchize all dimensions
        poles=ParArray[Pole](); //this makes a new, empty array.
        stride *= pointsInDimension;
        pointsInDimension = pointsPerDimension(dimension);
        var jump = stride * pointsInDimension;
        numberOfPoles = gridSize / pointsInDimension;
        for (i <- 0 to numberOfPoles-1 by 1){ // integer operations form bottleneck here -- nested loops are twice as slow
          val div = i / stride;
          start = div * jump + (i % stride);
          var p: Pole = new Pole();
          p.Pole(start, stride, pointsInDimension, dimension);
          poles:+p;
        }
        poles
        .foreach{pole => hierarchize1DUnoptimized(pole.start, pole.stride, pole.pointsInDimension, pole.dimension)};
      } // end loop over dimension 2 to d
    }

  /***
   * Hierarchizes the grid using a parallel unoptimized algorithm using the Java Parallel Stream framework.
   * This method creates numberOfBlocks poleBlock objects per dimension to hierarchize a subset of the poles.
   * Then is uses the parallelStream() method to run the hierarchization in parallel.
   * All parallism is handled by the Java framework.
   *
   * @param numberOfChunks The number of block objects to use
   */
    def hierarchizeUnoptimizedParallelStream(numberOfChunks:Int) {
        var stride = 1;
        var blocks: ParArray[PoleBlock]=ParArray[PoleBlock]();
  
  
      //dimension 1 separate as start of each pole is easier to calculate
      var pointsInDimension = pointsPerDimension(0);
      var numberOfPoles = gridSize / pointsInDimension;
      var polesPerBlock = numberOfPoles / numberOfChunks;
  
      for(i <- 0 to numberOfChunks-1 by 1) {
        var from = i * polesPerBlock;
        var to:Int=1; //Must be set.
        if (i+1 == numberOfChunks){
          to = numberOfPoles;
        } else {
          to = polesPerBlock * (i-1);
        }
        var pb:PoleBlock = new PoleBlock;
        pb.PoleBlock(1, pointsInDimension, 0, pointsInDimension, from, to);
        blocks:+ pb;
      }
  
      blocks
        .foreach(block => block.hierarchize());
  
      for(dimension <- 1 to dimensions-1 by 1) { // hierarchize all dimensions
        blocks=ParArray[PoleBlock](); //this makes a new, empty array.
        stride *= pointsInDimension;
        pointsInDimension = pointsPerDimension(dimension);
        val jump = stride * pointsInDimension;
        numberOfPoles = gridSize / pointsInDimension;
        polesPerBlock = numberOfPoles / numberOfChunks;
        for(i <- 0 to numberOfChunks-1 by 1) {
          var from = i * polesPerBlock;
          var to: Int = 0;
          if (i+1==numberOfChunks) {
            to = numberOfPoles;
          } else {
            to = polesPerBlock * (i+1);
          }
          var pb:PoleBlock = new PoleBlock;
          pb.PoleBlock(stride, pointsInDimension, dimension, jump, from, to);
          blocks:+pb;
        }
  
        blocks
          .foreach(block => block.hierarchize());
      }
      // end loop over dimension 2 to d
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



