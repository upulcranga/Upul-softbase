package Distributions;

import java.util.Random;


	
	public class QDistribution_Out extends QDistribution{
		
		private double gammaValue = 1000;
		
		/* CONSTRUCTORS */ 
		public QDistribution_Out(double[] arr)
		{
			length = arr.length;
			qDist = arr;
			pDist = getDegreeDist(qDist);
			random = new Random();
		}
		
		public double getExpectation()
		{

			double avg = 0;
			for(int i=0; i< qDist.length;i++)
			{
				avg = avg + (i*qDist[i]);
			}
			return avg;
		}
		
		public double getVariance()
		{
			double v = 0;
			double avg = getExpectation();
			for(int i=0; i< qDist.length;i++)
			{
				v = v + (i*i*qDist[i]);
			}
			v = v - (avg*avg);
			return v;
		}
	}

