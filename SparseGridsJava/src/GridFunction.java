public interface GridFunction {
	abstract double call(double[] x);
}

class allOnes implements GridFunction {
	public double call(double[] x) {
		return 1;
	}
}

class testFunc1D implements GridFunction {
	public double call(double[] x) {
		return 1 + x[0] - Math.abs(x[0] - 0.5) + 3 * x[0] * (x[0] - 0.3);
	}
}

class testFunc2D implements GridFunction {
	public double call(double[] x) {
		return 1 + Math.sin(x[0] * x[1]) + Math.cos(x[0] * x[0]) + x[1] * x[1] * x[1];
	}
}

class testFunc4D implements GridFunction {
	public double call(double[] x) {
		return 1 + 8 * (x[0] * x[1]) * Math.cos(x[2] + x[3]) * (4 * x[3]) + Math.cos(x[0] * x[0] + x[2]) * Math.sin(x[3]) * Math.sin(x[3]) * Math.cos(2 * x[2] + x[0]) * Math.sin(x[3]);
	}
}

class testFunc5D implements GridFunction {
	public double call(double[] x) {
		return 1 + Math.sin(x[0] * x[1]) * Math.cos(x[2] + x[3] * x[4]) * Math.sin(4 * x[3]) + Math.cos(x[0] * x[0] + x[2]) * Math.sin(x[4]) * Math.sin(x[4]) * Math.cos(2 *x[2] + x[0]) * Math.sin(x[3]) * Math.cos(x[4]);
	}
}

class parabel1D implements GridFunction {
	public double call(double[] x) {
		return 4 * x[0] * (1 - x[0]);
	}
}

class parabel2D implements GridFunction {
	public double call(double[] x) {
		return 4 * 4 * x[0] * (1 - x[0]) * x[1] * (1 - x[1]);
	}
}

class parabel3D implements GridFunction {
	public double call(double[] x) {
		return 64 * x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]);
	}
}

class parabel4D implements GridFunction {
	public double call(double[] x) {
		return 64 * 4 * x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * x[3] * (1 - x[3]);
	}
}
