


package Basic;

import java.awt.Color;
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
import oldClassesArchive.AssortativenessCalculatorLight;

import Distributions.EDistribution;
import Distributions.QDistribution;
import Distributions.QDistribution_General;
import Distributions.QDistribution_PowerLaw;

public class NodeAssortativityTester {
	
	public static void main(String[] args) {

	
		NodeAssortativityTester sim = new NodeAssortativityTester();
		sim.simulate();
	}
	
	
	public void simulate()
	{
		System.out.println("NodeAssort sim started");
		Random random = new Random();
		NetworkGenerator gen = new NetworkGenerator();
		Network network;

		
		BufferedWriter out2 = null;
		BufferedReader in2 = null;
	
	
	
		network =  gen.generateFromFile();
		//network = gen.generateFromMultipleFiles();
		//network = gen.generateERRandom(2000, 4000);
	
	
	
		double gamma;
		int Np; 
		
		BufferedWriter out = null;
	
		
		try {
		
			out = new BufferedWriter(new FileWriter(".\\..\\Output\\NODERHO\\IEEERESULTS\\uuu5.xls"));
			
		}
			
		 catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}

		
			
				try
				{
					System.out.println("number of nodes: " + network.getSize());
					System.out.println("number of links: " + network.getNoOfLinks());
				
					System.out.println("Maximum degree " + network.realMaxDegree());
					
					/* Measure r for the grown network */
					AssortativenessCalculatorLight ac = new AssortativenessCalculatorLight(network);
					double calculated_r = Math.round(1000 * ac.calculate_r_based_on_links()) / 1000.0;
					
					
					System.out.println("CALCULATED UNDIRECTED r: " + calculated_r);
					
					AssortativenessCalculator ac1 = new AssortativenessCalculator(network);
					calculated_r = Math.round(1000 * ac1.calculate_r_based_on_ejk()) / 1000.0;
					
					
					System.out.println("CALCULATED UNDIRECTED r: " + calculated_r);
				
					double GAMMA =  network.getGamma();
					System.out.println("GAMMA=" + GAMMA);
					
				
				
			
				
			

				//double[] rl_dist = ac.get_rdelta_dist_undirected();
				double[] rl_dist = ac.get_rnode_dist_undirected();
				
				boolean ZERO_OFF = false;
				

		

					try {
					
						
						for(int i=0; i< rl_dist.length;i++)
						{
						
							if( (!ZERO_OFF) || (ZERO_OFF && (rl_dist[i] !=0) ))
							{
							if(rl_dist[i]!=0)
							 {
								out.write(i+"\t"+rl_dist[i]+"\t");
								out.flush();
								out.write( System.getProperty(("line.separator")));
								out.flush();
							 }
							}
						}
						
				
						   
						out.flush();
						System.out.println("FINISHED WRITING");
						
					} 
					catch (IOException e)
					{
						System.out.println("cant write");
						e.printStackTrace();
					}
					
					//System.exit(0);
	
					//double[][] rl_distX = ac.get_rdelta_dist_undirected_xnode_EXTENDED();
					double[][] rl_distX = ac.get_rnode_dist_undirected_xnode_EXTENDED();
					
		

						try {
						
							
							
							for(int i=0; i< rl_distX.length;i++)
							{
								
								
									out.write(i+"\t"+rl_distX[i][0]+"\t"+rl_distX[i][1]+"\t"+rl_distX[i][2]+"\t"+rl_distX[i][3]+"\t");
									out.flush();
									out.write( System.getProperty(("line.separator")));
									out.flush();
			
							}
							
							out.flush();
							System.out.println("FINISHED WRITING");
							
						
							
						} 
						
						catch (IOException e) {
							System.out.println("cant write");
							e.printStackTrace();
						}
				

					
				}
				
				
					
					
				//Catch Exception in the Loop!!!
				catch(Exception e)
				{
					System.out.println("Exception Occured");
					e.printStackTrace();
					System.exit(0);
					
				}
				
				
			
				System.out.println("Writing ALL COMPLETE");
		
	
	}
	
	
	
	

	
	
	/** utility functions *
	 * @param x
	 * @return
	 */
	public static double weightedmean(double[] x) {
		double m = 0;
		for (int i = 0; i < x.length; i++) {
			m = m + ((i+1) * x[i]);
		}

		return m-1;
	}

	public static double variance(double[] x) {
		double m1 = 0;
		double m2 = 0;
		for (int i = 0; i < x.length; i++) {
			m1 = m1 + (i * i * x[i]);
			m2 = m2 + (i * x[i]);
		}

		return m1 - (m2 * m2);
	}
	
	public static double calculate_r(EDistribution eDist, QDistribution_PowerLaw qDist) {
		
		double[][] ejk = eDist.getDistribution();
		double[] qk = qDist.getDist();
		double r = 0;
		
		for (int j = 0; j < qk.length; j++)
		{
			for (int k = 0; k < qk.length; k++)
			{
				r = r + (j*k*( ejk[j][k] -(qk[j]*qk[k]) ));
			}

		}
		r = r /  variance(qk);
		
		return r;
	}
	
	public static double calculate_I(EDistribution eDist, QDistribution_PowerLaw qDist) {
		
		double[][] ejk = eDist.getDistribution();
		double[] qk = qDist.getDist();
		double i = 0;
		
		for (int j = 0; j < qk.length; j++)
		{
			for (int k = 0; k < qk.length; k++)
			{  double addition =  ejk[j][k] / (qk[j]*qk[k]);
				if(addition != 0)
				{
					i = i + (ejk[j][k])*Math.log(  ejk[j][k] / (qk[j]*qk[k])  );
				}
			}

		}
		 i  =  i /Math.log(2);
		
		return i;
	}
	
	
	public double calculateMIDifference(Network network1, Network network2)
	{
		double m = 0;
		double[][] ejk1 =  network1.getLinkDistribution();
		double[][] ejk2 =  network2.getLinkDistribution();
		double[] q1 = network1.getRemainingDegreeDistribution();
		double[] q2 = network2.getRemainingDegreeDistribution();
		double i1 = calculate_I(ejk1, q1);
		double i2 = calculate_I(ejk2, q2);
		
		return Math.abs(i1-i2);
	}
	
public static double calculate_I(double[][] ejk, double[] qk) {
		
		
		double i = 0;
		
		for (int j = 0; j < qk.length; j++)
		{
			for (int k = 0; k < qk.length; k++)
			{  double addition =  ejk[j][k] / (qk[j]*qk[k]);
				if(addition != 0)
				{
					i = i + (ejk[j][k])*Math.log(  ejk[j][k] / (qk[j]*qk[k])  );
				}
			}

		}
		 i  =  i /Math.log(2);
		
		return i;
	}

}
