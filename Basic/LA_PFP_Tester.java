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

import Distributions.EDistribution;
import Distributions.QDistribution;
import Distributions.QDistribution_General;
import Distributions.QDistribution_PowerLaw;

public class LA_PFP_Tester {
	
	public static void main(String[] args) {

	
		LA_PFP_Tester sim = new LA_PFP_Tester();
		sim.simulate();
	}
	
	
	public void simulate()
	{
		
		System.out.println("start simulation");
		Random random = new Random();
		NetworkGenerator gen = new NetworkGenerator();
		Network network;

		
		
		BufferedWriter out2 = null;
		BufferedReader in2 = null;

	network = gen.generateFromFile();
	//network = gen.growPFPNet(500, 0.042);
	//network = gen.growPARGNet(1000);
	//network = gen.growNetworkPA(3000, 1.1);
	//network = gen.growNetworkAPA(3000, 200, 2, 0.5);
	//network = gen.generateERRandom(700, 1500);
	//network.createTextFile("PARG  Net3000 0.4p");
	//network = gen.growIGNet(3000);
	
	//network.createTextFile(".\\PARGTEST\\APA ex.txt");
	System.out.println("net created /read");
				
		BufferedWriter out = null;
	
		
		try {
		
			out = new BufferedWriter(new FileWriter(".\\PARGTEST\\hpylori OLD LA DIST.xls"));
			
		}
			
		 catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}

		
			
			
			/*
			 * Build up a network of N nodes the APA way, for the given
			 * EDistribution and Q distribution
			 */
			
			//System.out.println("BUILDING NETWORK");
			
				
				try
				{
					
				
				
		
					/* Measure r for the grown network */
					SimpleAC ac = new SimpleAC(network);
					//AssortativenessCalculator acc =  new AssortativenessCalculator(network);
			
				
					System.out.println("Links " + ((network.getAllLinks().size())));
					System.out.println("Nodes " + ((network.getAllNodes().size())));
					System.out.println("Rich Club Coeffcient = " + network.richClubCoefficient(1.0));
					System.out.println("maximum degree = " + network.getCutoff());
					System.out.println("gamma = " + network.getGamma());
					System.out.println("r= " + ac.calculate_r() );
					//System.out.println("th r= " + acc.calculate_r_based_on_ejk() );
			
				
					
					
					
				double[] rl_dist = ac.get_rl_dist();
				//double[][] rlx = ac.get_rl_dist_xnodes_EXTENDED();;
				
				boolean ZERO_OFF = false;
			
					try {
					
						
						for(int i=0; i< rl_dist.length;i++)
						{
							if( (!ZERO_OFF) || (ZERO_OFF && (rl_dist[i] !=0) ))
							{
							
								out.write(i+"\t"+rl_dist[i]+"\t");
								//out.write(i+"\t"+rlx[i][3]+"\t"+rlx[i][1]+"\t");
								out.flush();
								out.write( System.getProperty(("line.separator")));
								out.flush();
							}
						}
						
						out.flush();
						System.out.println("FINISHED WRITING");
						
					} catch (IOException e) {
						System.out.println("cant write");
						e.printStackTrace();
					}
			
					
				}
				//Catch Exception in the Loop!!!
				catch(Exception e)
				{
					System.out.println("EXCEPTION !!! ");
					e.printStackTrace();
					System.exit(0);
					
				}
				
				
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
