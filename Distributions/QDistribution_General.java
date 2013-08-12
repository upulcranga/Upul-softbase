package Distributions;

import java.util.Random;


	
	public class QDistribution_General  extends QDistribution{
		
		private double gammaValue = 1000;
		
		/* CONSTRUCTORS */ 
		public QDistribution_General(double[] arr)
		{
			length = arr.length;
			qDist = arr;
			pDist = getDegreeDist(qDist);
			random = new Random();
		}
	}

