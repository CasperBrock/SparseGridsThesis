package gridFunctions;

public class TestFunc5D implements GridFunction {
	public double call(double[] x) {
		return 1 + Math.sin(x[0] * x[1]) * Math.cos(x[2] + x[3] * x[4]) * Math.sin(4 * x[3]) + Math.cos(x[0] * x[0] + x[2]) * Math.sin(x[4]) * Math.sin(x[4]) * Math.cos(2 *x[2] + x[0]) * Math.sin(x[3]) * Math.cos(x[4]);
	}
}
