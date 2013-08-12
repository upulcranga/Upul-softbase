package Distributions;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

public abstract class QDistribution {
	
	protected double[] qDist;
	protected double[] pDist;
	protected int length; /* Number of elements in the Distribution. Distribution must start from zero  */
	protected Random random;
	
	/**
	 * This function generates a random degree
	 * NOTE: Random remaining degree SHOULD NOT BE directly generated from the q distribution!
	 * @return a degre randomly generated from this distribution
	 */
	public int getRandomDegree()
	{
		//TODO: Do this in cleverer, faster way!
		int pd = 0;
		double ra = random.nextDouble();
		double total = 0;
		for(int i=0;i< length;i++)
		{
			total = total + pDist[i];
			if(total >= ra)
			{
				pd = i;
				break;
			}
			pd = length-1;
			
		}
		
		return pd;
	}
	
	
	
	
	
	
	public double getDegreeProbability(int degree)
	{
		return pDist[degree];
	}
	public double getRemainingDegreeProbability(int q)
	{
		if((q >= 0) && (q < length))
			return qDist[q];
		else
			return 0;
	}
	
	/**
	 * Calculate and return the average degree of this network based on this remaining degree distribution
	 */
	public double getAverageDegree()
	{

		double avg = 0;
		for(int i=0; i< pDist.length;i++)
		{
			avg = avg + (i*pDist[i]);
		}
		return avg;
	}
	
	public double getAverageRemainingDegree()
	{
		return getAverageDegree()-1;
	}
	
	/**
	 * Converts from remaining degree distribution to degree distribution
	 * @param qd
	 * @return
	 */
	public double[] getDegreeDist(double[] qd)
	{
		double sum = 0;
		for(int i=0; i< qd.length;i++)
		{
			sum = sum + ( qd[i] / (i+1) );
		}
		
		double[] pDist = new double[qd.length+1];
		for(int i=0; i< qd.length;i++)
		{
			pDist[i+1] = qd[i] / (sum*(i+1));
		}
		return pDist;
	}
	
	/**
	 * Converts from degree distribution to remaining degree distribition
	 * @param pd
	 * @return
	 */
	public double[] getRemainingDegreeDist(double[] pd)
	{
		double sum = 0;
		for(int i=0; i< pd.length;i++)
		{
			sum = sum + ( i*pd[i]  );
		}
		
		double[] qDist = new double[pd.length-1];
		for(int i=1; i< pd.length;i++)
		{
			qDist[i-1] = (i*pd[i]) /sum;
		}
		return qDist;
	}
	
	public int getLength()
	{
		return length;
	}
	
	/**
	 * Returns the maximum probability value in this distribution
	 * @return MAX value in Q dist
	 */
	public double getPeakQValue()
	{
		double max = 0;
		for(int i=0; i< qDist.length;i++)
		{
			max = Math.max(qDist[i],max);
		}
		return max;
	}
	
	/**
	 * Calculate the Entropy of this Q distribution, H(Q)
	 * Answer in Nits
	 */
	public double calculateHq() {
		double h = 0;
		for (int i = 0; i < qDist.length; i++) {
			if (qDist[i] != 0)
			{
				h = h - (qDist[i] * Math.log(qDist[i]));
			}
		}
		return h;
	}
	/**
	 * Get the Zeta number for the given index gamma, and number of elements n
	 * 
	 */
	public double zeta(double gamma, int n)
	{
		double zeta = 0;
		for (int i = 1; i <= (n); i++)
		{
			zeta = zeta + Math.pow(i,( -1.0*gamma));
		
		}
		
		return zeta;
	}
	
	/**
	 * This function determines the distribution length that is needed to generate a power
	 * law degree distribution with given gamma and average degree
	 */
	public int calculateDistLength(double gamma, double avg_k)
	{
		if(gamma <0.1){
			System.out.println("Invalid Gamma");
			return 0;
		}
		int n=1; /* We do not allow zero  degree */
		double z1 = 1;
		double z2 = 1;
		double z_ratio = 0;
		//double index = 1;
		while(true)
		{
			z1 = z1 + Math.pow(n, -1*gamma);
			z2 = z2 + Math.pow(n, -1*(gamma-1));
			z_ratio = (z2-1)/(z1-1); /* -1 because we do not want q=0, though its legal */
			if(z_ratio>avg_k)
			{
				return n;
			}
			if(n>100)
			{
				System.out.println("Impossible to get n");
				return n;
			}
			n++;
		}
	}
	
	/**
	 * This function (roughly) calculates the maximum average degree possible
	 * for  power law distribution with a given gamma,
	 * and a maximum allowed distribution length
	 * @return Maximum average degree possible for this distribution
	 */
	public double calculateMaximumAverageDegree(double gamma, int max_dist_length)
	{
		double ak = zeta(gamma-1,1)/zeta(gamma,1);
		
		 ak = (zeta(gamma-1,max_dist_length)-1)/(zeta(gamma,max_dist_length)-1);
		return ak;
	}
	
	public double[] getDist()
	{
		return qDist;
	}
	
	public void print()
	{
		System.out.println("");
		System.out.println("q");
		for (int i = 0; i<qDist.length; i++)
		{
			System.out.print(""+(Math.round(1000*qDist[i]) / 1000.0)+" ");
		}
		System.out.println("");
		System.out.println("");
	}
	
	
	public void printToFile(String Filename)
	{
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(".\\PRE\\"+Filename+".xls"));
			
			for (int i = 0; i<qDist.length; i++)
			{
				out.write(""+(i) +"\t");
				out.write(""+(Math.round(100000*pDist[i]) / 100000.0)+"\t");
				out.write(""+(Math.round(100000*qDist[i]) / 100000.0)+"\t\n");
				out.flush();
			}
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
	}
	

	

}
