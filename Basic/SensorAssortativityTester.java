package Basic;


import java.io.BufferedWriter;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.Vector;

import oldClassesArchive.AssortativenessCalculator;

import Distributions.EDistribution;
import Distributions.QDistribution;
import Distributions.QDistribution_General;
import Distributions.QDistribution_PowerLaw;
public class SensorAssortativityTester {

	/**
	 * @param args
	 * 
	 * This class experiments with some basic concepts regarding 
	 * assortativity of non-degree quantities
	 * For example, sensor values of nodes in sensor networks
	 * This class will use the existing classes for q(k) and e(j,k),
	 * but with the different interpretation, where the values j,k etc
	 * represent sensor values, and NOT degrees or remaining degrees
	 * Created first: 18th September 2009
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		SensorAssortativityTester sim = new SensorAssortativityTester();
		sim.simulate();

	}
	
	
	public void simulate()
	{
		
		System.out.println("start simulation");
		Random random = new Random();
		NetworkGenerator gen = new NetworkGenerator();
		Network network;

		BufferedWriter out1 = null;
		BufferedWriter out2 = null;
		BufferedWriter out3 = null;
	
		
		try {
		
			out1 = new BufferedWriter(new FileWriter(".\\..\\Output\\CONG\\PA 500 cong e beginning.xls"));
			out2 = new BufferedWriter(new FileWriter(".\\..\\Output\\CONG\\PA 500 cong e  MI.xls"));
			out3 = new BufferedWriter(new FileWriter(".\\..\\Output\\CONG\\PA 500 cong e end .xls"));
			
		}
			
		 catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}
		
	
		BufferedReader in = null;
		
		//Grow a standard preferential attachment network
		network = gen.growNetworkPA(500, 1.5);
		//network = gen.generateFromFile();
		network.assignRandomBinaries();
		network.initialize_pagerank_random(); //needed only if we are going to use pageRank
		//network = gen.growNetworkAPA(3000, 200, 2.1, 0.0);
		//network.createTextFile("APAnet r=0.0 real  ");
		//network = gen.generateFullyConnected(100);
		
		//assign sensor readings to nodes
		//Lets assign randomly
//		int minReading = 500;
//		
//		int maxReading = 550;
//		
//		for(int i=0; i< network.getSize(); i++)
//		{
//			Node thisNode = network.getAllNodes().elementAt(i);
//			thisNode.setReading(minReading+ (maxReading - minReading)*Math.abs(random.nextDouble()));
//			//thisNode.setReading(minReading+ (maxReading - minReading)*0.5);
//			
//		}
//		
//		/* Measure r for the grown network */
		AssortativenessCalculator ac = new AssortativenessCalculator(network);
		AssortativenessCalculator ac2 = null;
		double calculated_r = Math.round(1000 * ac.calculate_nondegree_r()) / 1000.0;
		calculated_r = Math.round(1000 * ac.calculate_r_based_on_ejk()) / 1000.0;
		//double calculated_rOUT = Math.round(1000 * ac.calculate_theoritical_directed_r_OO()) / 1000.0;
		//double calculated_rIN = Math.round(1000 * ac.calculate_theoritical_directed_r_InIn()) / 1000.0;
//	
//
		System.out.println("CALCULATED degree based r: " + calculated_r);
		//System.out.println("CALCULATED degree based r out: " + calculated_rOUT);
		//System.out.println("CALCULATED degree based r in: " + calculated_rIN);
		System.out.println("CALCULATED Np: " + network.getCutoff());
		System.out.println("CALCULATED gamma: " + network.getGamma());
		System.out.println("CALCULATED N: " +network.getSize());
		System.out.println("CALCULATED links: " +network.getNoOfLinks());
		//network.createTextFile(".\\SensorAssort\\APA -0.8");
		

		
	
	calculated_r = Math.round(1000 * ac.calculate_theoritical_r_sensornet()) / 1000.0;
//		
//
	System.out.println("CALCULATED NON-DEGREE r (theoritical): " + calculated_r);
//		
		double[] rl_dist = ac.get_rl_e_dist_sensor();
		
//		
	try {
//			
//			
//			
//			
		for(int i=0; i< rl_dist.length;i++)
		{
			if( true)
			{
//				
				out1.write(i+"\t"+rl_dist[i]+"\t");
					out1.flush();
					out1.write( System.getProperty(("line.separator")));
				out1.flush();
			}
			}
//			
		out1.flush();
			System.out.println("FINISHED WRITING");
//			
	} 
	catch (IOException e)
		{
			System.out.println("cant write");
			e.printStackTrace();
		}
		

		
	
		for(int i = 0; i < 20; i++)
		{
		
		double oneProb = 0.2 + 0.8*((double)(i%80)/80.0) ;
		oneProb = 0.5;
		
		

		network.simNet(1);

		
	

		double readingSD =  network.getpMean_sensor(); //TODO: Printing MEAN: CHANGE BACK TO STANDARD DEVIATION
		
		//double readingSD =  network.countOnes() / (double)network.getSize();
		
		/* Measure r for the grown network */
		ac = new AssortativenessCalculator(network);
		calculated_r = Math.round(1000 * ac.calculate_nondegree_r()) / 1000.0;
	

		//System.out.println("CALCULATED NON-DEGREE r: " + calculated_r);
		
		calculated_r = ac.calculate_theoritical_r_sensornet(); //multiplication for Excel
		
		double mi_sensor = network.calculate_MI_straight_sensor();
		
		//double h1 = network.calculate_H1_straight_sensor();
		//double h2 = network.calculate_H2_straight_sensor();
		
		int pos =  ac.getPositiveLinks();
		int neg =  network.getNoOfLinks() - pos;
		
		System.out.println("CALCULATED NON-DEGREE r (theoritical): " + readingSD+ " "+calculated_r+ " "+mi_sensor+" "+pos+" "+neg);
		


		try{
			out2.write(i+"\t"+calculated_r+"\t"+readingSD + "\t"+mi_sensor+"\t"+pos+"\t"+neg);
			out2.flush();
			out2.write( System.getProperty(("line.separator")));
			out2.flush();
		}
		catch (IOException e)
		{
			System.out.println("cant write to file 2");
			e.printStackTrace();
		}
		
		}
		
		System.out.println("FIUNISJED PRINTING");
		
		int positiveLinks = 0;
		int negativeLinks = 0;
		
		
		 rl_dist = ac.get_rl_e_dist_sensor();
		
		 double[][] rl =  ac.get_rl_e_dist_sensor_MULTI();
		try {
			
			
			
			
			for(int i=0; i< rl.length;i++)
			{
				if( true)
				{
				
					out3.write(i+"\t"+rl[i][0]+"\t"+rl[i][1]+"\t"+rl[i][2]+"\t");
					out3.flush();
					out3.write( System.getProperty(("line.separator")));
					out3.flush();
					if(rl[i][0] > 0.0)
					{
						positiveLinks++;
					}
					else
					{
						negativeLinks++;
					}
				}
			}
			
			out3.flush();
			System.out.println("FINISHED WRITING");
			
			System.out.println(" pos neg links "+ positiveLinks +" "+ negativeLinks);
			
		} 
		catch (IOException e)
		{
			System.out.println("cant write");
			e.printStackTrace();
		}
		
		
	}
	

	
	
	

}
