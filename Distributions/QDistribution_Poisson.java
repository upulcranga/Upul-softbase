package Distributions;

import java.util.Random;
import java.util.Vector;

import Basic.Link;
import Basic.Network;

/**
 * This class creates a QDistributions
 * COORESPPONDING TO A POISSON DEGREE DISTRIBUTION Pk
 * It DOES NOT create a poisson QDistribtuion
 * Poisson Degree Distribution results from Erdos - Renyi Random Graphs.
 * @author pir011
 *
 */
public class QDistribution_Poisson  extends QDistribution{
	
	private final double VERY_SMALL_PROBABILITY = Math.pow(10, -10); //Any probability less than this is made zero
	/* CONSTRUCTORS */
	
	
	public QDistribution_Poisson(int l, double lambda)
	{
			length = l;
			pDist = new double[l+1];
			double sum = 0;
			pDist[0] = 0; //No 'hanging' nodes
			if(l <= 100) //Poisson approximation is fine
			{
				for(int i=1;i<(l+1);i++) 
				{
					
					double x = (Math.exp(-1*lambda)*Math.pow(lambda, i-1)) / factorial(i-1);
					if(x < VERY_SMALL_PROBABILITY)
					{
						x = 0; //Avoid meaningless calculations
					}
					 sum = sum +  x;
					 if(sum >1.0)
					 {
						 sum = sum - pDist[i-1];
						 break;
					 }
					 pDist[i] = x;
				
				}
			}
			else
			{
				//Normal distribution
				for(int i=1;i<(l+1);i++) 
				{
					double sigma = 4.0;
					double x = (1.0 / (sigma*Math.sqrt(2*Math.PI)))*(Math.exp(-1*(1.0/(sigma*sigma))*(0.5*(i-lambda)*(i-lambda))));
					if(x < VERY_SMALL_PROBABILITY)
					{
						x = 0; //Avoid meaningless calculations
					}
					 sum = sum +  x;
					 if(sum >1.0)
					 {
						 sum = sum - pDist[i-1];
						 break;
					 }
					 pDist[i] = x;
				
				}
			}
			
			
			if(sum <= 1.0)
			{
				 pDist[l] = 1.0 - sum;
			}
			else
			{
				System.out.println("DEBUG CREATION OF EXPONENTIAL DISTRIBUTION");
			}
			
			qDist = getRemainingDegreeDist(pDist);
			random = new Random();
		
	}
	
	public double factorial(double x)
	{
		if(x <=1)
		{
			return 1;
		}
		else
		{
			return x*factorial(x-1);
		}
	}
	
	
	
	

}
