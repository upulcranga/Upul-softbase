package Basic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;


import oldClassesArchive.AssortativenessCalculator;
import Basic.Network;
import Basic.NetworkGenerator;
import Basic.Node;

public class NetGenCPA {

	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		NetworkGeneratorRefined gen = new NetworkGeneratorRefined();

		Network net = null;

		long startTime = System.currentTimeMillis();
		
			try {

				net = gen.growNetwork(10, 0.5);
				

				long endTime = System.currentTimeMillis();
				System.out.println("Took " + (endTime - startTime) + " ms");

				net.createTextFile("f://PA_100.txt");

				//System.out.println("Net size nodes " + net.getSize() + " Links " + net.getNoOfLinks());

			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		
		System.out.println("DONE!");
	}

}



