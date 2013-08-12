package Distributions;

import java.util.Random;


/**
 * This is a Utility class that has all E Distribution Generating Functions
 * This SHOULD BE THE ONLY PLACE where any method that generates an EDistribution should be implemented
 * Any changes to generating these also SHOULD BE MADE ONLY HERE!
 * @author pir011
 *
 */
public class EDistributionGenerator
{
	/*
	 * Default Constructor
	 */
	
	private final double VERY_SMALL_PROBABILITY = Math.pow(10, -10); //Any probability less than this is made zero
	public EDistributionGenerator()
	{
		
	}
	
	public  EDistribution generateEDistribution_class_A(QDistribution qDist, double r,
			int totslidingSteps) {
		Random random = new Random();
		EDistribution EDist1 = generateAssortativeEDistribution(qDist);
		
		//EDistribution EDistM1 = generateNewmanDisassortativeEDistribution(qDist); 
		EDistribution EDistM1  = generateDisassortativeEDistribution(qDist);
		double rm = calculate_r(EDistM1, qDist);
		System.out.println("rm"+rm);
		EDistribution EDist = generateAssortativeEDistribution(qDist);
		EDistribution EDist0 = generateNonassortativeEDistribution(qDist);
		int totalOneSideSlidingSteps = totslidingSteps / 2;
		int slidingSteps;
		if (r >= 0) {

			slidingSteps = (int) Math
					.round(((double) totalOneSideSlidingSteps * (1.0 - r)));
		} else {
			slidingSteps = (int) Math
					.round(((double) totalOneSideSlidingSteps * (1.0 + r)));
		}

		/*
		 * slide down by 'slidingSteps' number of steps from the r=1
		 * eDistribution
		 */

		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		if (r >= 0) /* Sliding from r=1 to r=0 */
		{
			eArr = EDist.getDistribution();
			for (int i = 0; i < qDist.getLength(); i++) {

				double stepSize = (EDist1.getProbability(i, i) - EDist0
						.getProbability(i, i))
						/ (double) totalOneSideSlidingSteps;
				eArr[i][i] = eArr[i][i] - (stepSize * slidingSteps);

				for (int j = 0; j < eArr[i].length; j++) {

					if (i == j) {
						continue;
					}

					eArr[i][j] = eArr[i][j]
							+ (((double) slidingSteps / (double) totalOneSideSlidingSteps) * EDist0
									.getProbability(i, j));

				}
			}
		} else /* Sliding from r=0 to r=-1 */
		{

			double[][] eArrM1 = EDistM1.getDistribution();
			double[][] eArr0 = EDist0.getDistribution();
			
			for (int i = 0; i < qDist.getLength(); i++) 
			{
				for (int j = 0; j < qDist.getLength(); j++) 
				{
					eArr[i][j] = eArrM1[i][j] - ((double) slidingSteps / (double) totalOneSideSlidingSteps)*(( eArrM1[i][j]- eArr0[i][j]));
				
				}
				
			}

		

		
		}
		
		//Clear up and make small numbers zero
		for (int i = 0; i < qDist.getLength(); i++) 
		{
			for (int j = 0; j < qDist.getLength(); j++) 
			{
				if(eArr[i][j] < VERY_SMALL_PROBABILITY)
				{
					eArr[i][j] = 0;
				}
			
			}
		}

		EDistribution newEDist = new EDistribution(eArr);

		return newEDist;
	}
	
	public  EDistribution generateEDistribution_class_B(QDistribution qDist, double r,
			int totslidingSteps) {
		Random random = new Random();
		EDistribution EDist1 = generateAssortativeEDistribution(qDist);
		EDistribution EDist = generateAssortativeEDistribution(qDist);
		
		EDistribution EDistM1 = generateDisassortativeEDistribution(qDist);
		//EDistribution EDistM1 = generateNewmanDisassortativeEDistribution(qDist);
		
	
	
		int totalOneSideSlidingSteps = totslidingSteps / 2;
		int slidingSteps;
		if (r >= 0) {

			slidingSteps = (int) Math
					.round(((double) totalOneSideSlidingSteps * (1.0 - r)));
		} else {
			slidingSteps = totalOneSideSlidingSteps  + (int) Math
					.round(((double) totalOneSideSlidingSteps * (Math.abs(r))));
		}

		/*
		 * slide down by 'slidingSteps' number of steps from the r=1
		 * eDistribution
		 */

		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		

			double[][] eArrM1 = EDistM1.getDistribution();
			double[][] eArr1 = EDist1.getDistribution();
			
			for (int i = 0; i < qDist.getLength(); i++) 
			{
				for (int j = 0; j < qDist.getLength(); j++) 
				{
					eArr[i][j] = eArr1[i][j] - ((double) slidingSteps / (double) totslidingSteps)*(( eArr1[i][j]- eArrM1[i][j]));
				
				}
				
			}
			
			
		EDistribution newEDist = new EDistribution(eArr);

		return newEDist;
	}
	
	
	
	public EDistribution generateAssortativeEDistribution(QDistribution qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		for (int y = 0; y < qDist.getLength(); y++) {
			
			eArr[y][y] = qDist.getDist()[y];
			
		}

		EDistribution ed = new EDistribution(eArr);
		return ed;
	}
	
	public EDistribution generateNonassortativeEDistribution(QDistribution qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		for (int j = 0; j < qDist.getLength(); j++) {
			for (int k = 0; k < qDist.getLength(); k++) {
				eArr[j][k] = qDist.getDist()[j]* qDist.getDist()[k];
			}
		}

		EDistribution ed = new EDistribution(eArr);
		return ed;
	}
	
	/**
	 * This is a Crude Heuristic Used to get a dis assortative E Distribution
	 * @param qDist
	 * @return
	 */
	public static EDistribution generateDisassortativeEDistribution(QDistribution qDist) {
		
		
		//GA_DisAssortative_Distribution_Generator gaGen = new GA_DisAssortative_Distribution_Generator(qDist);
		
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
	
		int lastIndex = qDist.getLength() - 1;
		double Nq = qDist.getLength();
		

		for (int y = 0; y < qDist.getLength(); y++) 
		{
			if (qDist.getRemainingDegreeProbability(y) <= qDist
					.getRemainingDegreeProbability(lastIndex - y)) {
				/*
				 * when both indices are the 'opposite', there is peak
				 * probability
				 */
				eArr[y][lastIndex - y] = qDist
						.getRemainingDegreeProbability(y);
				/*
				 * Other elements are zero, but no need to explicitly code that
				 * as arrays are automatically initialized to zero
				 */
			} else /* Need to maintain symmetry */
			{
				eArr[y][lastIndex - y] = qDist.getRemainingDegreeProbability(lastIndex - y);
				eArr[y][y] =  qDist.getRemainingDegreeProbability(y) - qDist.getRemainingDegreeProbability(lastIndex - y);	

			}
		}

		EDistribution ed = new EDistribution(eArr);
		boolean ch = checkEDistValidity(ed, qDist);
		if(!ch){
			System.out.println("invalid Er=-1 which is:");
			ed.print();
		}
		else
		{
			//System.out.println("Valid Er=-1 which is:");
			//ed.print();
		}
		//return gaGen.generate();
		return ed;
	}
	
	

	public static boolean checkQDistValidity(QDistribution_PowerLaw q) {
		double[] qArr = q.getDist();

		double sum = 0;
		for (int j = 0; j < qArr.length; j++) {
			sum = sum + qArr[j];
		}
		if ( Math.abs(sum - 1) >= 0.001) { /*Giving room for double imprecision */
			System.out.println("q Sum false");
			return false;
		}
		return true;
	}


	public static boolean checkEDistValidity(EDistribution ed, QDistribution q) {
		double[][] eArr = ed.getDistribution();
		int counter = 0;
		boolean check = true;
		double sum = 0;
		for (int j = 0; j < eArr.length; j++) {
			double rowsum = 0;
			for (int k = 0; k < eArr.length; k++) {
				if (eArr[j][k] < 0) {

					//System.out.println(" < 0");
					return false;

				}

				if (eArr[j][k] > 1) {
					System.out.println(" > 1");
					return false;
				}

				if (eArr[k][j] != eArr[j][k]) {
					System.out.println("non symmetry");
					return false;
				}

				rowsum = rowsum + eArr[j][k];
				sum = sum + eArr[j][k];

			}

			if ( Math.abs(rowsum - q.getRemainingDegreeProbability(j)) > 0.001) { /*Giving room for double imprecision */
				System.out.println("Rowsum false at"+j);
				return false;
			}

		}
		
		if ( Math.abs(sum - 1) > 0.001) { /*Giving room for double imprecision */
			System.out.println("Sum false");
			return false;
		}

		return check;
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
	
	public static double calculate_r(EDistribution eDist, QDistribution qDist) {
		
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
	
	/**
	 * Dissassortativeness as defined by Newman is used in this method
	 * Crude Approximation - works ONLY FOR Np = 7 OR Np = 4!!!!!!
	 * @param qDist
	 * @return
	 */
	public EDistribution generateNewmanDisassortativeEDistribution(
			QDistribution_PowerLaw qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];

		double[] xrand2 = new double[qDist.getLength()];
		double minR = 0;
		double[] xma;
		double res = 0.005;
		boolean validE = false;
		if (qDist.getLength() == 7) {
			//System.out.println("length = seven");
			for (double i1 = 0; i1 < 1.0; i1 = i1 + res) {
				xrand2[0] = i1;
				if ((xrand2[0] > 1 + (0.5 * res))) {

					break;
				}
				for (double i2 = 0; i2 < 1.0; i2 = i2 + res) {
					xrand2[1] = i2;
					if ((xrand2[0] + xrand2[1]) > 1 + (0.5 * res)) {

						break;
					}
					for (double i3 = 0; i3 < 1.0; i3 = i3 + res) {
						xrand2[2] = i3;
						if ((xrand2[0] + xrand2[1] + xrand2[2]) > 1 + (0.5 * res)) {

							break;
						}
						for (double i4 = 0; i4 < 1.0; i4 = i4 + res) {
							xrand2[3] = i4;
							if ((xrand2[0] + xrand2[1] + xrand2[2] + xrand2[3]) > 1 + (0.5 * res)) {

								break;
							}
							for (double i5 = 0; i5 < 1.0; i5 = i5 + res) {
								xrand2[4] = i5;
								if ((xrand2[0] + xrand2[1] + xrand2[2]
										+ xrand2[3] + xrand2[4]) > 1 + (0.5 * res)) {

									break;
								}
								for (double i6 = 0; i6 < 1.0; i6 = i6 + res) {
									xrand2[5] = i6;
									if ((xrand2[0] + xrand2[1] + xrand2[2]
											+ xrand2[3] + xrand2[4] + xrand2[5]) > 1 + (0.5 * res)) {

										break;
									}
									{
										for (double i7 = 0; i7 < 1.0; i7 = i7
												+ res) {
											xrand2[6] = i7;
											{

												if ((xrand2[0] + xrand2[1]
														+ xrand2[2] + xrand2[3]
														+ xrand2[4] + xrand2[5] + xrand2[6]) < 1 - (0.5 * res)) {

													continue;
												} else if ((xrand2[0]
														+ xrand2[1] + xrand2[2]
														+ xrand2[3] + xrand2[4]
														+ xrand2[5] + xrand2[6]) > 1 + (0.5 * res)) {

													break;
												} else {
													xma = xrand2.clone();
													EDistribution Ed = generateNewmanEdist(
															qDist.getDist(),
															xrand2);
													boolean b = checkEDistValidity(
															Ed, qDist);
													if (b) {

														double mq = weightedmean(qDist
																.getDist());
														double mx = weightedmean(xrand2);
														double sq2 = variance(qDist
																.getDist());
														// calculate rd - NewmAN FORMULA
														double rd = -((mq - mx) * (mq - mx))
																/ sq2;

														//Standard formula
														rd = calculate_r(Ed, qDist);
														if (rd < minR) {
															validE = true;
															minR = rd;
															eArr = Ed.getDistribution();
															
														}

													}
													double checker = 0;
												}
											}
										}
									}

								}
							}
						}
					}
				}
			}
		}
		
		else if (qDist.getLength() == 4) {
			//System.out.println("length = four");
			for (double i1 = 0; i1 < 1.0; i1 = i1 + res) {
				xrand2[0] = i1;
				if ((xrand2[0] > 1 + (0.5 * res))) {

					break;
				}
				for (double i2 = 0; i2 < 1.0; i2 = i2 + res) {
					xrand2[1] = i2;
					if ((xrand2[0] + xrand2[1]) > 1 + (0.5 * res)) {

						break;
					}
					for (double i3 = 0; i3 < 1.0; i3 = i3 + res) {
						xrand2[2] = i3;
						if ((xrand2[0] + xrand2[1] + xrand2[2]) > 1 + (0.5 * res)) {

							break;
						}
						
									{
										for (double i4 = 0; i4 < 1.0; i4 = i4
												+ res) {
											xrand2[3] = i4;
											{

												if ((xrand2[0] + xrand2[1]
														+ xrand2[2] + xrand2[3]
														) < 1 - (0.5 * res)) {

													continue;
												} else if ((xrand2[0]
														+ xrand2[1] + xrand2[2]
														+ xrand2[3] ) > 1 + (0.5 * res)) {

													break;
												} else {
													xma = xrand2.clone();
													EDistribution Ed = generateNewmanEdist(
															qDist.getDist(),
															xrand2);
													boolean b = checkEDistValidity(
															Ed, qDist);
													if (b) {

														double mq = weightedmean(qDist
																.getDist());
														double mx = weightedmean(xrand2);
														double sq2 = variance(qDist
																.getDist());
														// calculate rd - NewmAN FORMULA
														double rd = -((mq - mx) * (mq - mx))
																/ sq2;

														//Standard formula
														rd = calculate_r(Ed, qDist);
														if (rd < minR) {
															validE = true;
															minR = rd;
															eArr = Ed.getDistribution();
															
														}

													}
													double checker = 0;
												}
											}
									
							}
						}
					}
				}
			}
		}
		if (!validE) {
			System.out.println("ZERO EJK");
		}
		return new EDistribution(eArr);
	}

	public  EDistribution generateNewmanEdist(double[] q, double[] x) {
		double[][] eArr = new double[q.length][q.length];
		for (int j = 0; j < eArr.length; j++) {
			for (int k = 0; k < eArr.length; k++) {
				eArr[j][k] = (q[j] * x[k]) + (q[k] * x[j]) - (x[j] * x[k]);
				
			}

		}
		EDistribution ed = new EDistribution(eArr);
		return ed;
	}
	
	

}
