package grid;

public class CombiGridBuilder {

	/***
	 * Returns an anisotropic CombiGrid with the given number of dimensions and
	 * the given size. The grid will be anisotropic with higher refinement in the lower dimensions.
	 * 
	 * @param size The size of the grid as the sum of the level vector. The size must be atleast 2 * dimensions.
	 * @param dimensions The number of dimensions for the grid. Must be atleast 2.
	 * @return A CombiGrid with the given number of dimensions and size.
	 */
	public static CombiGrid anisotropicGrid(int size, int dimensions) {
		if(dimensions < 2)
			throw new IllegalArgumentException("A grid must have atleast 2 dimensions");
		if(size / dimensions < 2)
			throw new IllegalArgumentException("Size must be atleast twice the number of dimensions");
		int rem = size - (dimensions * 2);
		int[] levels = new int[dimensions];
		for(int i = 0; i < dimensions; i++) {
			int used;
			if(rem == 1) {
				used = 1;
				rem = 0;
			}
			else if(i + 1 == dimensions) {
				used = rem;
				rem = 0;
			}
			else {
				used = (int)Math.ceil(rem / 2.0);
				rem = rem - used;
			}
			levels[i] = 2 + used;
		}

		return new CombiGrid(levels);
	}
	
	/***
	 * Returns an anisotropic CombiGrid with the given number of dimensions and
	 * the given size. The grid will be anisotropic with higher refinement in the lower dimensions.
	 * 
	 * @param size The size of the grid as the sum of the level vector. The size must be atleast 2 * dimensions.
	 * @param dimensions The number of dimensions for the grid. Must be atleast 2.
	 * @return A CombiGrid with the given number of dimensions and size.
	 */
	public static CombiGridAligned anisotropicAlignedGrid(int size, int dimensions) {
		if(dimensions < 2)
			throw new IllegalArgumentException("A grid must have atleast 2 dimensions");
		if(size / dimensions < 2)
			throw new IllegalArgumentException("Size must be atleast twice the number of dimensions");
		int rem = size - (dimensions * 2);
		int[] levels = new int[dimensions];
		for(int i = 0; i < dimensions; i++) {
			int used;
			if(rem == 1) {
				used = 1;
				rem = 0;
			}
			else if(i + 1 == dimensions) {
				used = rem;
				rem = 0;
			}
			else {
				used = (int)Math.ceil(rem / 2.0);
				rem = rem - used;
			}
			levels[i] = 2 + used;
		}

		return new CombiGridAligned(levels, 32);
	}

	/***
	 * Returns an isotropic CombiGrid with the given number of dimensions
	 * close to the given size. The grid will be fully isotropic so a 
	 * mismatch in size and dimensions will result in a smaller grid than requested.
	 * 
	 * @param size The size of the grid as the sum of the level vector. The size must be atleast 2 * dimensions.
	 * @param dimensions The number of dimensions for the grid. Must be 2 or higher.
	 * @return A CombiGrid with the given number of dimensions and a size of
	 * (size / dimensions) * dimensions
	 */
	public static CombiGrid isotropicGrid(int size, int dimensions) {
		if(dimensions < 2)
			throw new IllegalArgumentException("A grid must have atleast 2 dimensions");
		if(size / dimensions < 2)
			throw new IllegalArgumentException("Size must be atleast twice the number of dimensions");
		int level = size / dimensions;
		int[] levels = new int[dimensions];
		for(int i = 0; i < dimensions; i++) {
			levels[i] = level;
		}

		return new CombiGrid(levels);
	}

	/***
	 * Returns a isotropic CombiGridAligned with the given number of dimensions
	 * close to the given size. The grid will be fully isotropic so a 
	 * mismatch in size and dimensions will result in a smaller grid than requested.
	 * 
	 * @param size The size of the grid as the sum of the level vector. The size must be atleast 2 * dimensions
	 * @param dimensions The number of dimensions for the grid. Must be atleast 2
	 * @return A CombiGridAligned with the given number of dimensions and a size of
	 * (dimensions / size) * dimensions
	 */
	public static CombiGridAligned isotropicAlignedGrid(int size, int dimensions) {
		if(dimensions < 2)
			throw new IllegalArgumentException("A grid must have atleast 2 dimensions");
		if(size / dimensions < 2)
			throw new IllegalArgumentException("Size must be atleast twice the number of dimensions");
		int level = size / dimensions;
		int[] levels = new int[dimensions];
		for(int i = 0; i < dimensions; i++) {
			levels[i] = level;
		}

		return new CombiGridAligned(levels, 32);
	}
}
