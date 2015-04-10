package grid;

public class Content { //object for holding both int and array.

	public int asInt;
	public int[] l = new int[8]; //This will keep 0 or 1 in all positions, for comparing.
	
	public Content(){
		
	}
	
	public Content(int asIntInput, int[] lInput){
		asInt = asIntInput;
		System.arraycopy(lInput, 0, l, 0, lInput.length);
	}
}
