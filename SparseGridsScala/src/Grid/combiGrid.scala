package Grid

import scala.collection.Parallel
import scala.collection.parallel.mutable.ParSeq
import scala.collection.parallel.mutable.ParArray

object CombiGrid {

	//global variables
	var levels:Array[Int] = Array(0);
  var dimensions:Int = 0;
  var grid:ParArray[Double] = ParArray(0);
  var gridSize:Int = 0;
  var pointsPerDimension:Array[Int] = Array(0);


//Main. Use to launch.
def main(args: Array[String]): Unit = {
		levels= Array(5,5,5,5);
		CombiGrid(levels);
		FillArrayWithOnes(grid);
    hierarchizeUnoptimized();  
		PrintArray(grid);
		println(gridSize)
		//pointsPerDimension.foreach(l => print(""+l+'\t'))
}



//Small helper functions
def myPow2(j: Int): Int = { //Doubles input
		return 1<<j;
}

  def compare(grid2:ParArray[Double]): Boolean = { //Compares two grids
    return grid.equals(grid2)
  }

//Create a Grid, based on levels. Will return with unset values.
def CombiGrid(levels: Array[Int]){
	dimensions = levels.length
			pointsPerDimension = new Array[Int](dimensions)
			this.levels = new Array[Int](dimensions)
			var size = 1
			for(i <- 0 to dimensions-1) {//a <- 1 to 3;
				pointsPerDimension(i) = myPow2(levels(i)) - 1
						size *= pointsPerDimension(i)
			}
	println(dimensions)
	gridSize = size
	grid = new ParArray[Double](gridSize)
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
        //System.out.println("Start: " + start + '\t' + "Stride: " + stride);
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
      offset = myPow2(levels(dimension) - level) - 1;
      parentOffset =  stepSize;
      //stepSize = stepSize << 1;
      stepSize = stepSize * 2;
    }

    right = 0.5 * grid(start + (offset + parentOffset) * stride);
    grid(start + offset * stride) -= right;
    offset += stepSize;
    grid(start + offset * stride) -= right;
  }



}