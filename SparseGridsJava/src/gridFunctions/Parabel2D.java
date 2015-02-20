package gridFunctions;

public class Parabel2D implements GridFunction {
	public double call(double[] x) {
		return 4 * 4 * x[0] * (1 - x[0]) * x[1] * (1 - x[1]);
	}
}
