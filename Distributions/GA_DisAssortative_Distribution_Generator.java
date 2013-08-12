package Distributions;
import GA.*;

public class GA_DisAssortative_Distribution_Generator {
	
	QDistribution qDist;
	public GA_DisAssortative_Distribution_Generator(QDistribution q)
	{
		qDist = q;
	}
	
	/**
	 * This method generates the EDistribution with the most possible dis assortativeness
	 * (smallest value of r_min)
	 * for the given Remaining Degree Distribution q(k)
	 * @return
	 */
	public EDistribution generate()
	{
		
		GARunner ru = new GARunner(qDist);
		Genotype best = ru.runGA(100);
		int[] arr = best.getCodeSet();
		//Let us calculate the eDist from this x.
//		Make a double array of x
		double[] x = new double[arr.length];
		for(int i=0;i<x.length;i++)
	    {
			 x[i] = arr[i] / 10000.0; //TODO: Get rid of this magic number
	    }
		
		double[] q = qDist.getDist();
		//Use Newman Formula
		double[][] eArr = new double[q.length][q.length];
		for (int j = 0; j < eArr.length; j++) {
			for (int k = 0; k < eArr.length; k++)
			{
				eArr[j][k] = (q[j] * x[k]) + (q[k] * x[j]) - (x[j] * x[k]);
				//eArr[j][k] = (q[q.length-1-j] * x[q.length-1-k]) + (q[q.length-1-k] * x[q.length-1-j]) - (x[q.length-1-j] * x[q.length-1-k]);	
			}

		}
		EDistribution ed = new EDistribution(eArr);
		//Calculate The  Assortativeness for this EDist
		EDistributionGenerator egg = new EDistributionGenerator();
		double r = egg.calculate_r(ed,qDist);
		System.out.println("best rmin "+r);
		return ed;
	}
	
	
	

}
