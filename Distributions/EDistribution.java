package Distributions;

import java.util.Random;
import java.util.Vector;

import Basic.Link;
import Basic.Network;

public class EDistribution {
	
	private double[][] eDist; 
	private double[][] ijEDist; /* the distribution of i*j*eDist[i][j] */
	private int length; /* Number of elements in the Distribution. Distribution must start from zero. Both dimensions must have the same length  */
	private Random random;
	
	public EDistribution(double[][] arr)
	{
		length = arr.length;
		eDist = arr;
		random = new Random();
	}
	
	public EDistribution(Network network)
	{
		setDistribution(network);
		random = new Random();
	}
	
	
	public void setDistribution(Network network){
		/* Calculate the e
		 *  distribution from the network */
		double[][] newDist = new double [network.getMaxRemainingDegree()+1][network.getMaxRemainingDegree()+1]; /*Note degree will be from zero */
		
		//Go through all links
		for(int i=0;i< network.getAllLinks().size();i++)
		{
			Link thisLink = network.getAllLinks().elementAt(i);
			int qi = thisLink.iRemainingDegree();
			int qj = thisLink.iRemainingDegree();
			newDist[qi][qj] = newDist[qi][qj] + 1;
		}
		
		/* Normalize the distribution */
		for(int i=0;i< newDist.length;i++)
		{
			for(int j=0;j< newDist[i].length;j++)
			{
				newDist[i][j] = newDist[i][j] / (network.getAllLinks().size()); 
			}
		}
		
		/* Set it as the new qistribution */
		eDist = newDist;
		length = newDist.length;
		
	}
	
	public double getProbability(int i, int j)
	{
		if((i >= 0) && (i < length) && (j >= 0) && (j < length))
			return eDist[i][j];
		else
			return 0;
	}
	
	public double[][] getDistribution()
	{
		return eDist;
	}
	
	/* returns i*j*Edist[i][j] Distribution */
	public double[][] getWeightedDistribution()
	{
		double[][] ijEDist = new double[eDist.length][eDist.length];
		for(int i=0;i< eDist.length;i++)
		{
			for(int j=0;j< eDist[i].length;j++)
			{
				ijEDist[i][j] = i*j*eDist[i][j];
			}
		}
		return ijEDist;
	}
	
	/**
	 * Return the scaled distribution, where each row is a probability disribution
	 * @return
	 */
	public double[][] getScaledDistribution()
	{
		double[][] scaledEDist = new double[eDist.length][eDist.length];
		for(int i=0;i< eDist.length;i++)
		{
			double sum = 0;
			for(int j=0;j< eDist[i].length;j++)
			{
				sum = sum + eDist[i][j];
			}
			
			for(int j=0;j< eDist[i].length;j++)
			{
				scaledEDist[i][j] = eDist[i][j] /sum;
			}
		}
		return scaledEDist;
	}
	
	public double Expectation()
	{
		double ex = 0;
		for(int i=0;i< eDist.length;i++)
		{
			for(int j=0;j< eDist.length;j++)
			{
				
				ex = ex + i*j*eDist[i][j];
			}
		}
		
		return ex;
		
	}
	
	public void print()
	{
		for (int i = 0; i<eDist.length; i++)
		{
			for (int j = 0; j<eDist.length; j++)
			{
				System.out.print(""+(Math.round(1000*eDist[i][j]) / 1000.0)+" ");
			}
			System.out.println("");
		}
		System.out.println("");
	}
	
		
}
