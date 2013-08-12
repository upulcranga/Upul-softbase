package Basic;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

public class BatchNetGenenator {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
	 	Random random = new Random(System.currentTimeMillis());
	    	BufferedWriter out = null;
	    	try{
	    	
	    	out = new BufferedWriter(new FileWriter(".\\LA-Extended\\NET.xls"));
	    	}
	    	
	    	 catch (IOException e) {
	 			System.out.println("cant write to file");
	 			e.printStackTrace();
	 		}
	    	 NetworkGenerator gen =  new NetworkGenerator();
	    	 
	    	 for(int i = 0; i< 100; i++)
	    	 {
	    	 
	    	
		    	 Network net =  gen.generateERRandom(500, 1500);
		    	 net =  gen.growPARGNet(500);
		    	 net =  gen.growNetworkPA(500, 0.5 + (0.01*i));
		    	 net.createTextFile("net" + i);
	    	 }

	}

}
