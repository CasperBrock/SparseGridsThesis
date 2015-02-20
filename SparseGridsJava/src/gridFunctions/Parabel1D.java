package gridFunctions;

public class Parabel1D implements GridFunction {
	public double call(double[] x) {
		return 4 * x[0] * (1 - x[0]);
	}
}
