package Basic;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.Vector;

import oldClassesArchive.AssortativenessCalculator;

import Distributions.EDistribution;
import Distributions.QDistribution_PowerLaw;

public class LocalrSimulator {
	
	public static void main(String[] args) {

	
		LocalrSimulator sim = new LocalrSimulator();
		sim.simulate();
	}
	
	
	public void simulate()
	{
		Random random = new Random();
		Network network;

		int SLIDINGSTEPS = 20;

		/* Initialize Parameters */
		double ra = 0;
		int N = 1000; //Used to be 7000
		
		double gamma =1;
		int Np = 10; //This has to be relatively small compared to the size of network
		//So that the probabilities can be calculated properly
		
		

		/* Get the Q distribution for the desired average degree */
		QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(Np, gamma, true);



		/* Get Edistribution for desired r and Q distribution */
		EDistribution eDist = generateEDistribution(qDist, ra,
				SLIDINGSTEPS);

		
		BufferedWriter out = null;
		
		try {
		
			out = new BufferedWriter(new FileWriter(".\\LocalR\\Newdef\\localr_ distBin;r="+ra+".xls"));
		}
			
		 catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}

		for(int y=0; y < 1; y++){ 
			System.out.println("y="+y);
			
			
			/*
			 * Build up a network of N nodes the APA way, for the given
			 * EDistribution and Q distribution
			 */
			
			//System.out.println("BUILDING NETWORK");
			int avgTimes = 0;
			while(avgTimes < 1)
			{
				
				try
				{
					network = buildAPANetwork(qDist, eDist, N);
					//Network clonedNetwork = network.clone();
					//System.out.println("NETWORK BUILT");
					
				
				
		
					/* Measure r for the grown network */
					AssortativenessCalculator ac = new AssortativenessCalculator(network);
					double calculated_r = Math.round(1000 * ac.calculate_r()) / 1000.0;
					
					
					System.out.println("CALCULATED R: " + calculated_r);
					double calculated_sum_r = Math.round(1000 * ac.calculate_r_as_Sum_of_Rl()) / 1000.0;
					double calculated_avg_k = Math.round(1000 * network.getAverageDegree()) / 1000.0;
					
					System.out.println("CALCULATED SUM R: " + calculated_sum_r);
					
					
					// calculate r using the 'theoritical' corrrelation method
					double calculated_th_r = Math.round(1000 * ac.calculate_r_based_on_ejk()) / 1000.0;
					
					
					System.out.println("CALCULATED theory R: " + calculated_th_r);
					//System.out.println("Links " + ((network.getAllLinks().size())));
					//System.out.println("Nodes " + ((network.getAllNodes().size())));
					
				boolean ZERO_OFF  =false;
				double[] rl_dist = ac.get_rl_dist_bin();
				//double[] rl_dist = ac.get_rl_dist_xnodes();
				//double[] rl_dist = ac.get_ejkl_dist_xnodes();
					try {
					
						/*
						out.write(""+N+"\t"+Np+"\t"+ra+"\t"+calculated_r+"\t"+gamma);
						out.flush();
						out.write( System.getProperty(("line.separator")));
						out.flush();*/
						
						for(int i=0; i< rl_dist.length;i++)
						{
							if( (!ZERO_OFF) || (ZERO_OFF && (rl_dist[i] !=0) ))
							{
								out.write(i+"\t"+rl_dist[i]+"\t");
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
					avgTimes++;
					
				}
				//Catch Exception in the Loop!!!
				catch(Exception e)
				{
					System.out.println("EXCEPTION AT Y="+y+" trial="+avgTimes);
					continue; //TODO: Debug this Exception
				}
				
				//Print Network for visualization
				network.createTextFile("SimNet gamma=1 r=0 np=10 N= 1000");
				
				
			}
		}
	}
	
	
	/**
	 * This method implements the Assortative Preferential Attachment Method
	 * 
	 * @param qDist
	 * @param eDist
	 * @param networkSize
	 * @return
	 */
	public Network buildAPANetwork(QDistribution_PowerLaw qDist, EDistribution eDist,
			int networkSize) {
		Random random = new Random();

		/* start with m0 nodes. 1:2 ratio between m0 and N is optimal */
		int m0 = (int) (0.5 * networkSize);

		Network network = new Network();

		/* Add half of the nodes with no links */
		/* Gennerate otehr half of the nodes, which will be joining one by one */
		int counter = 0;
		Vector<Node> joiningNodes = new Vector();
		for (int y = 0; y < qDist.getLength(); y++) {

			int nodesWithThisq = (int) Math.round((m0 * qDist
					.getDegreeProbability(y + 1)));
			if (y == qDist.getLength() - 1)// *If the last one, get the balance
											// of nodes */
			{
				nodesWithThisq = m0 - network.getSize();
			}
			for (int k = 1; k <= nodesWithThisq; k++) {
				counter++;
				Node node = new Node(counter, (y + 1));
				network.addNode(node);
				/* Generate another node, with the same degree which will join */
				Node node2 = new Node(m0 + counter, (y + 1));
				joiningNodes.addElement(node2);
			}

		}

		/* Shuffle the set of nodes that are going to join */
		joiningNodes = network.shuffle(joiningNodes);

		/* Add this half of nodes in the 'APA' way */

		int addedLinksCount = 0;
		boolean linkLimit;

		for (int i = 0; i < (m0); i++) /*
										 * Note that m0 = N/2 gives the optimum
										 * conditions
										 */
		{
			linkLimit = false;
			addedLinksCount = 0;

			/* Add new node */
			
			Node node = joiningNodes.elementAt(i); 
			network.addNode(node);

			/* ADD LINK */
			/* scale the eDist for this instance of adding */
			double[] scaledDist = getScaledEDist(node.getAssignedRemainingDegree(),
					eDist, qDist);

			/*
			 * Go through all existing nodes in the network in Random order, and
			 * assign probabilities
			 */
			network.shuffle();
			for (int j = 0; j < network.getSize(); j++)
			// while(true)
			{
				Node networkNode = network.getAllNodes().elementAt(j);
				// Node networkNode = network.getAllNodes().elementAt(
				// (Math.abs(random.nextInt())%network.getSize()) ); //TODO:
				// Check for duplicate links
				double nodeProb = scaledDist[networkNode.getAssignedRemainingDegree()];
				// double nodeProb =
				// scaledDist2[node.getRemainingDegree()][networkNode.getRemainingDegree()];
				if (networkNode.getNumberOfLinks() >= networkNode.getAssignedDegree()) {
					/* all 'Stumps' have been exhausted */
					nodeProb = 0;
				}

				/*
				 * if(networkNode.isLinked(node)) { // This network node is
				 * already linked to joining node nodeProb = 0; }
				 */

				if (networkNode.getID() == node.getID()) {
					/* We dont want to connect the joining ndoe to itself! */
					nodeProb = 0;
				}

				/*
				 * Toss a coin to see whether this node should connect to the
				 * new node
				 */
				if (random.nextDouble() < nodeProb) {
					try {
						/* Add Link */
						network.addLink(node, networkNode);
						addedLinksCount++;
						if (addedLinksCount >= node.getAssignedDegree()) {
							/* Should not add any more links */
							linkLimit = true;
							addedLinksCount = 0;
						}
					} catch (Exception e) {
						// Eat it!
						System.out.println("Problem adding links");
						linkLimit = true;
						addedLinksCount = 0;
					}
				}

				if (linkLimit) /*
								 * No need to add any more links to this joining
								 * node
								 */
				{
					break;
				}

			}

		
		}
		
		System.out.println("Network size before cleanup "+network.getSize());
		//Clean network from any 'hanging' nodes i.e nodes with zero degree
		for (int i = 0; i < network.getSize(); i++)
		{
			Node thisNode = network.getAllNodes().elementAt(i);
			if(thisNode.getNumberOfLinks() <= 0) //Not linked
			{
				try {
					boolean b = network.removeNode((int) thisNode.getID());
					if(b){
						//System.out.println("Linkless Node REMOVED");
					}
					else
					{
						System.out.println("Linkless Node COULD NOT BE REMOVED");
					}
				} catch (Exception e) {
				
					e.printStackTrace();
				}
				
				
			}
		}
		
		System.out.println("Network size after cleanup "+network.getSize());
		// network.printArray(network.getLinkDistribution(),"link dist");
		return network;

	}
	
	/**
	 * 'Scale' a Edistribution
	 */
	private double[] getScaledEDist(int joiningNodeQ, EDistribution eDist,
			QDistribution_PowerLaw qDist) {
		if (joiningNodeQ < 0) {
			System.out.println("A Node Tries to join with Unsuitable remaining Degree: Error");
		}
		double[] scaledDist = new double[eDist.getDistribution().length];

		/* Copy the values from main E distribution */
		double probSum = 0;
		for (int i = 0; i < scaledDist.length; i++) {
			
				scaledDist[i] = (eDist.getProbability(joiningNodeQ, i) / qDist
						.getDegreeProbability(i + 1));
			

			probSum = probSum + scaledDist[i];
		}

		/*
		 * Scale the copied Distribution so that summation of its probabilities
		 * become 1
		 */
		for (int i = 0; i < scaledDist.length; i++) {
			scaledDist[i] = scaledDist[i] / probSum;
		}

		return scaledDist;
	}

	public static EDistribution generateEDistribution(QDistribution_PowerLaw qDist, double r,
			int totslidingSteps) {
		Random random = new Random();
		EDistribution EDist1 = generateAssortativeEDistribution(qDist);
		
		EDistribution EDistM1  = generateDisassortativeEDistribution2(qDist);
		
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

		EDistribution newEDist = new EDistribution(eArr);

		return newEDist;
	}
	
	
	public static EDistribution generateEDistribution2(QDistribution_PowerLaw qDist, double r,
			int totslidingSteps) {
		Random random = new Random();
		EDistribution EDist1 = generateAssortativeEDistribution(qDist);
		EDistribution EDistM1 = generateDisassortativeEDistribution2(qDist);
		EDistribution EDist = generateAssortativeEDistribution(qDist);
	
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
			
			if(r==0){
			
				/*
				for (int i = 0; i < qDist.getLength(); i++) 
				{
					for (int j = 0; j < qDist.getLength(); j++) 
					{
						if(i==j)
						{
							eArr[i][j] = 0.006;
						}
						else if(i==(qDist.getLength()- j-1))
						{
							eArr[i][j] = 0.244;
						}
						else
						{
							eArr[i][j] = 0;
						}
					
					}
					
				}
				*/
			
			}

	
		EDistribution newEDist = new EDistribution(eArr);

		return newEDist;
	}
	
	
	
	public static EDistribution generateAssortativeEDistribution(QDistribution_PowerLaw qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		for (int y = 0; y < qDist.getLength(); y++) {
			
			eArr[y][y] = qDist.getDist()[y];
			
		}

		EDistribution ed = new EDistribution(eArr);
		return ed;
	}
	
	public static EDistribution generateNonassortativeEDistribution(QDistribution_PowerLaw qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		for (int j = 0; j < qDist.getLength(); j++) {
			for (int k = 0; k < qDist.getLength(); k++) {
				eArr[j][k] = qDist.getDist()[j]* qDist.getDist()[k];
			}
		}

		EDistribution ed = new EDistribution(eArr);
		return ed;
	}
	
	public static EDistribution generateDisassortativeEDistribution2(QDistribution_PowerLaw qDist) {
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


	public static boolean checkEDistValidity(EDistribution ed, QDistribution_PowerLaw q) {
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
				System.out.println("Rowsum false");
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
