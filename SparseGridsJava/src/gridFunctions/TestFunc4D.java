package gridFunctions;

public class TestFunc4D implements GridFunction {
	public double call(double[] x) {
		return 1 + 8 * (x[0] * x[1]) * Math.cos(x[2] + x[3]) * (4 * x[3]) + Math.cos(x[0] * x[0] + x[2]) * Math.sin(x[3]) * Math.sin(x[3]) * Math.cos(2 * x[2] + x[0]) * Math.sin(x[3]);
	}
}
