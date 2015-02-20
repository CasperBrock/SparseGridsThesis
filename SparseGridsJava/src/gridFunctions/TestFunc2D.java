package gridFunctions;

public class TestFunc2D implements GridFunction {
	public double call(double[] x) {
		return 1 + Math.sin(x[0] * x[1]) + Math.cos(x[0] * x[0]) + x[1] * x[1] * x[1];
	}
}