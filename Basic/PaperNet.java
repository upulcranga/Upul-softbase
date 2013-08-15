package Basic;

import java.math.*;

public class PaperNet extends Network {

	public PaperNet() {
		super();
		double impactFactor = 0;
		int numberOfCitations = generateRandom(10,50);
		int [] authorArray = new int[generateRandom(1,5)];
		
	}
	
	public int generateRandom(int min, int max){
		return min + (int)(Math.random() * ((max - min) + 1)); //return a random number of citations to be initiated.
	}
	
	

}
