package Basic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Random;

public class NetFileParser {

	/**
	 * @param args
	 * This class exists to do some parsing operations to modify network files
	 */
	public static void main(String[] args) {
		
		System.out.println("Executed class is Net File Parser");
		
		Random random = new Random(System.currentTimeMillis());
		NetworkGenerator gen = new NetworkGenerator();
		


		
		BufferedWriter out2 = null;
		BufferedReader in2 = null;
	
	
	
		String filename = ".\\..\\Data\\schoolnet\\sd02.txt";
		
		Network net = new Network();
		FileProcessor  x= new FileProcessor();
	    x.SetFile(filename);
	    x.Process(); 
	    x.showall();
	    x.buildlinks_conditional(30);
	    net = x.getNetwork();
	    
		net.createTextFile(".\\..\\Data\\schoolnet\\Schoolnet I above 30");
		
		System.out.println("File parsed and written");
		
	
	}

}
