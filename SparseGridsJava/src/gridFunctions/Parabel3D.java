package gridFunctions;

public class Parabel3D implements GridFunction {
	public double call(double[] x) {
		return 64 * x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]);
	}
}
