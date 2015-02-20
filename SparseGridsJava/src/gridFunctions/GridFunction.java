package gridFunctions;
/*public interface GridFunction {
	abstract double call(double[] x);
}*/

public class GridFunction {
	public static double call(double[] x, GridFunctions gf) {
		switch (gf) {
		case ALLONES:
			return 1;

		case TESTFUNC1D:
			return 1 + x[0] - Math.abs(x[0] - 0.5) + 3 * x[0] * (x[0] - 0.3);

		case TESTFUNC2D:
			return 1 + Math.sin(x[0] * x[1]) + Math.cos(x[0] * x[0]) + x[1] * x[1] * x[1];
		
		case TESTFUNC4D:
			return 1 + 8 * (x[0] * x[1]) * Math.cos(x[2] + x[3]) * (4 * x[3]) + Math.cos(x[0] * x[0] + x[2]) * Math.sin(x[3]) * Math.sin(x[3]) * Math.cos(2 * x[2] + x[0]) * Math.sin(x[3]);

		case TESTFUNC5D:
			return 1 + Math.sin(x[0] * x[1]) * Math.cos(x[2] + x[3] * x[4]) * Math.sin(4 * x[3]) + Math.cos(x[0] * x[0] + x[2]) * Math.sin(x[4]) * Math.sin(x[4]) * Math.cos(2 *x[2] + x[0]) * Math.sin(x[3]) * Math.cos(x[4]);
		
		case PARABEL1D:
			return 4 * x[0] * (1 - x[0]);
			
		case PARABEL2D:
			return 4 * 4 * x[0] * (1 - x[0]) * x[1] * (1 - x[1]);
			
		case PARABEL3D:
			return 64 * x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]);
			
		case PARABEL4D:
			return 64 * 4 * x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * x[3] * (1 - x[3]);
			
		case PARABELANYD: 
					int dim = x.length; //TODO no way to set this in current architecture.
					double result=Math.pow(2, dim);
					for (int i=0; i<dim;i++){
						result=result*x[i]*(1-x[i]);
					}
					return result;
		
		default:
			return 1;
		}
	}
}
