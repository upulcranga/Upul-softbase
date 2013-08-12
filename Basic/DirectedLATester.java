


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

public class DirectedLATester {
	
	public static void main(String[] args) {

	
		DirectedLATester sim = new DirectedLATester();
		sim.simulate();
	}
	
	
	public void simulate()
	{
		System.out.println("LA sim started");
		Random random = new Random();
		NetworkGenerator gen = new NetworkGenerator();
		Network network;

		
		BufferedWriter out2 = null;
		BufferedReader in2 = null;
	
	
	
	network = gen.generateFromFile();
	
	
	
		double gamma;
		int Np; 
		
		BufferedWriter out = null;
	
		
		try {
		
			//out = new BufferedWriter(new FileWriter(".\\..\\Output\\RHO-E\\ER out difference dist 2.xls"));
			out = new BufferedWriter(new FileWriter(".\\..\\Output\\KENNETH\\net2010 t.xls"));
			
		}
			
		 catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}

		
			
				try
				{
					System.out.println("number of nodes: " + network.getSize());
					System.out.println("number of links: " + network.getNoOfLinks());
				
		
					/* Measure r for the grown network */
					AssortativenessCalculator ac = new AssortativenessCalculator(network);
					double calculated_r = Math.round(1000 * ac.calculate_r_based_on_ejk()) / 1000.0;
					
					System.out.println("CALCULATED UNDIRECTED r: " + calculated_r);
					
					
					// calculate r using the 'theoritical' corrrelation method
					double calculated_th_r = Math.round(1000 * ac.calculate_theoritical_directed_r()) / 1000.0;
					
					
					System.out.println("CALCULATED theory R: " + calculated_th_r);
					
					//calculate r using the 'theoritical' corrrelation method
					//calculated_th_r = Math.round(1000 * ac.calculate_r_based_on_ejk()) / 1000.0;
					
					
					//System.out.println("CALCULATED UNDIRECTED theory R: " + calculated_th_r);
					
					calculated_th_r = Math.round(100000 * ac.calculate_theoritical_directed_r_OO()) / 100000.0;
					
					
					System.out.println("CALCULATED theory R out out: " + calculated_th_r);
					
					calculated_th_r = Math.round(1000 * ac.calculate_theoritical_directed_r_InIn()) / 1000.0;
					
					
					System.out.println("CALCULATED theory R in in: " + calculated_th_r);
			
			
			
					//Print Np and <k> of networks
					//double maxD =  network.getCutoff();
					//System.out.println("CutOff: " + maxD);
					
					
					double GAMMA =  network.getGamma();
					System.out.println("GAMMA=" + GAMMA);
					
				
					//double  mi =  network.calculate_MI_straight();
					//System.out.println("Mutual Information: " + mi);
					
					
					//mi =  network.calculate_MI_DirectedNet();
					//System.out.println("Mutual Information Directed: " + mi);
					
					//mi =  network.calculate_MI_Out_DirectedNet();
					//System.out.println("Mutual Information out-out: " + mi);
					
					// mi =  network.calculate_MI_In_DirectedNet();
					//System.out.println("Mutual Information In-In: " + mi);
					
					out.write(""+calculated_r);
					out.flush();
					out.write( System.getProperty(("line.separator")));
					out.flush();
					
					
					out.write(""+ac.calculate_theoritical_directed_r());
					out.flush();
					out.write( System.getProperty(("line.separator")));
					out.flush();
					
					out.write(""+ac.calculate_theoritical_directed_r_OO());
					out.flush();
					out.write( System.getProperty(("line.separator")));
					out.flush();
					
					out.write(""+ac.calculate_theoritical_directed_r_InIn());
					out.flush();
					out.write( System.getProperty(("line.separator")));
					out.flush();
					
					System.exit(0);
					
				
			
					
				//double[] rl_dist = ac.get_rho_e_out_dist();
				double[] rl_dist = ac.get_rho_e_out_dist();
			
				boolean ZERO_OFF = false;
				
				double posnodecount = 0;
				double negnodecount = 0;
				double posrhosum = 0;
				double negrhosum = 0;
				
		

					try {
					
				
						
						
						for(int i=0; i< rl_dist.length;i++)
						{
							if(rl_dist[i]>=0)
							{
								posnodecount = posnodecount + 1;
								posrhosum =  posrhosum + rl_dist[i];
								
							}
							else
							{
								negnodecount = negnodecount + 1;
								negrhosum =  negrhosum  + rl_dist[i];
							}
							if( (!ZERO_OFF) || (ZERO_OFF && (rl_dist[i] !=0) ))
							{
							
								out.write(i+"\t"+rl_dist[i]+"\t");
								out.flush();
								out.write( System.getProperty(("line.separator")));
								out.flush();
							}
						}
						
						out.write(" positive \t"+posnodecount+"\t" + (posrhosum /posnodecount ));
						out.write( System.getProperty(("line.separator")));
						out.flush();
						
						out.write(" negative \t"+negnodecount+"\t" + (negrhosum /negnodecount ));
						out.write( System.getProperty(("line.separator")));
						out.flush();
						
						out.write(" r \t"+(posrhosum+ negrhosum) );
						out.write( System.getProperty(("line.separator")));
						out.flush();
						   
						out.flush();
						System.out.println("FINISHED WRITING");
						
					} 
					catch (IOException e)
					{
						System.out.println("cant write");
						e.printStackTrace();
					}
					
					System.exit(0);
	
					double[][] rl_distX = ac.get_rho_e_in_dist_xlinks_EXTENDED();
					
		

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
