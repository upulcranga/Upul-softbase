package Distributions;

import java.util.Random;
import java.util.Vector;

import Basic.Link;
import Basic.Network;

/**
 * This class creates a QDistribution
 * COORESPPONDING TO A EXPONENTIAL DEGREE DISTRIBUTION Pk
 * It DOES NOT create a Exponential QDistribtuion
 * See newman for the instances where this kind of degree distribution happens in real world
 * @author pir011
 *
 */
public class QDistribution_Exponential extends QDistribution{
	
	private final double VERY_SMALL_PROBABILITY = Math.pow(10, -10); //Any probability less than this is made zero
	/* CONSTRUCTORS */
	
	/**
	 * NOTE: l should be infinity, ideally, for a proper probability distribution that somes up to ONE.
	 * If not, at least Lambda must be much smaller than l (Np) 
	 */
	public QDistribution_Exponential(int l, double lambda)
	{
			length = l;
			pDist = new double[l+1];
			double sum = 0;
			pDist[0] = 0; //No 'hanging' nodes
				
			for (int i = 1; i < (l + 1); i++)
			{
				
				double x = (1 - Math.exp(-1.0 / lambda))*Math.exp((-1.0*(i-1)) / lambda);
				if (x < VERY_SMALL_PROBABILITY) {
					x = 0; // Avoid meaningless calculations
				}
				sum = sum + x;
				if (sum > 1.0) {
					sum = sum - pDist[i - 1];
					break;
				}
				pDist[i] = x;
	
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
