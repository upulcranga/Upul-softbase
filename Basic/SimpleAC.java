package Basic;

/**
 * This CLASS is used ONLY to calculate LOCAL ASSORTATIVENESS DISTRIBUTIONS
 * for UNDIRECTED NETWORKS.
 * 
 * It cannot be used to calculate GLOBAL ASSORTATIVENESS.
 * It cannot be used with DIRECTED NETWORKS to calculate anything.
 */
import Basic.Link;
import Basic.Node;
import Basic.Network;
import Distributions.EDistribution;
import Distributions.QDistribution_In;
import Distributions.QDistribution_Out;
import Distributions.QDistribution_PowerLaw;

public class SimpleAC {
	
	private Network network;
	private QDistribution_PowerLaw qDist;
	private QDistribution_PowerLaw qDist2;

	private QDistribution_In qDistIn;
	private QDistribution_Out qDistOut;
	private QDistribution_In qDistInDash;
	private QDistribution_Out qDistOutDash;
	private EDistribution eDist;
	boolean removeLinkless  =false;
	
	public SimpleAC(Network net)
	{
		network = net;
		System.out.println("Network size before cleanup cal "+network.getSize());
		//Clean network from any 'hanging' nodes i.e nodes with zero degree
		for (int i = 0; i < network.getSize(); i++)
		{
			Node thisNode = network.getAllNodes().elementAt(i);
			if(thisNode.getNumberOfLinks() <= 0) //Not linked
			{
				if(removeLinkless)
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
		}
	
		
		System.out.println("Network size after cleanup cal "+network.getSize());
		qDist = new QDistribution_PowerLaw(net);
		//qDist2 = new QDistribution_PowerLaw(qDist.getDegreeDist(qDist.getDist()), true);
		//qDistIn = new QDistribution_In(net.getInDegreeDistribution());
		//qDistOut = new QDistribution_Out(net.getOutDegreeDistribution());
		//qDistInDash = new QDistribution_In(net.getInDegreeDistribution_dash());
		//qDistOutDash = new QDistribution_Out(net.getOutDegreeDistribution_dash());
		qDistIn = new QDistribution_In(net.getDegreeDistribution());  //TODO: QUICK FIX- GET RID OF THESE!!!!
		qDistOut = new QDistribution_Out(net.getDegreeDistribution()); //TODO: QUICK FIX- GET RID OF THESE!!!!
		qDistInDash = new QDistribution_In(net.getDegreeDistribution()); //TODO: QUICK FIX- GET RID OF THESE!!!!
		qDistOutDash = new QDistribution_Out(net.getDegreeDistribution()); //TODO: QUICK FIX- GET RID OF THESE!!!!
		//qDistIn = new QDistribution_In(net.getInRemainingDegreeDistribution());
		//qDistOut = new QDistribution_Out(net.getOutRemainingDegreeDistribution());
		eDist = new EDistribution(net);
	}
	
	/* The logic of this function is contained in
	 * "Assortative Mixing in Networks"
	 * MEJ newman
	 * Formula 4
	 */
	public double calculate_r()
	{
		if(qDist.getLength()<=2)
		{
			return 1; /* we do not allow q=0. Therefore there is only one kind of remaining degree. Perfect Matching */
		}
		double r=0;
		
		double sum1 = 0;
		double sum2 = 0;
		double sum3 = 0;
		
		double  M = network.getAllLinks().size(); /* Using the Notation in Newman for network Size */
		if(M==0)
		{
			return 0;
		}
		
		for(int i=0;i<M;i++)
		{
			/* are the j,k in Newman remaining degrees on each end, or degrees as the paper says? */
			Link link = network.getAllLinks().elementAt(i);
			//sum1 = sum1 + (link.iRemainingDegree()*link.jRemainingDegree());
			//sum2 = sum2 + ( 0.5*(link.iRemainingDegree() + link.jRemainingDegree()) );
			//sum3 = sum3 + ( 0.5*(/**/ (link.iRemainingDegree()*link.iRemainingDegree()) + (link.jRemainingDegree()*link.jRemainingDegree())/**/ )  );
			
			sum1 = sum1 + (link.iDegree()*link.jDegree());
			sum2 = sum2 + ( 0.5*(link.iDegree() + link.jDegree()) );
			sum3 = sum3 + ( 0.5*(/**/ (link.iDegree()*link.iDegree()) + (link.jDegree()*link.jDegree())/**/ )  );
		}
		
		if( (sum3/M) - (Math.pow((sum2/M),2))==0 ) /* Denominator Zero */
		{
			//System.out.println("Uncalculable r condition"); 
			/*TODO: Research the configurations that may case this. 
			The only known configuration now is two 2---2 links.
			For this, obviously  r=1 */
			if(M==2) /* ASSUMING THE ABOVE CONDITION */
			{
				r=1.0;
			}
			else
			{
				//System.out.println("Uncalculable r condition: UNKNOWN CONFIGURATION");
				r = 1000;
			}
		}
		else
		{
			r = ( (sum1/M) - Math.pow((sum2/M),2)  ) / ( (sum3/M) - Math.pow((sum2/M),2) );
		}
		
		return r;
	}
	
	public double calculate_averageND(Node node)
	{
		if(node.getNumberOfLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		
		
		
		int sumDegrees = 0;
		for(int y=0; y<node.getLinks().size();y++)
		{
			if(node.getLinks().elementAt(y).iNodeID() != node.getID())
			{
				//This is a neghbour
				sumDegrees = sumDegrees + node.getLinks().elementAt(y).iNode().getNumberOfLinks();
			}
			
			if(node.getLinks().elementAt(y).jNodeID() != node.getID())
			{
				//This is a neghbour
				sumDegrees = sumDegrees + node.getLinks().elementAt(y).jNode().getNumberOfLinks();
			}
		}
		
		double avgD = ((double)sumDegrees )/ ((double)node.getNumberOfLinks());
		
	
		
		return avgD;
	}
	
	
	public double calculate_localR(Node node)
	{
		if(node.getNumberOfLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		double N = network.getSize();
		double Nlinks = network.getNoOfLinks();
		double mew_q = qDist.getExpectation();
		double var_q = qDist.getVariance();
		//Now calculate the average remaining degree of all neighbours!
		
		if(var_q == 0)
		{
			// Only one type of degree!
			return (1.0 / N);
		}
		
		int sumRemainingDegrees = 0;
		for(int y=0; y<node.getLinks().size();y++)
		{
			if(node.getLinks().elementAt(y).iNodeID() != node.getID())
			{
				//This is a neghbour
				sumRemainingDegrees = sumRemainingDegrees + node.getLinks().elementAt(y).iRealRemainingDegree();
			}
			
			if(node.getLinks().elementAt(y).jNodeID() != node.getID())
			{
				//This is a neighbour
				sumRemainingDegrees = sumRemainingDegrees + node.getLinks().elementAt(y).jRealRemainingDegree();
			}
		}
		
		double avgRD = ((double)sumRemainingDegrees )/ ((double)node.getNumberOfLinks());
		
	
		double j = node.getNumberOfLinks()-1;
		double ejkl = (j*(j+1)*avgRD)/(2*Nlinks);
		//double mew_ql = ((j+1)/(2*Nlinks))*mew_q*mew_q;
//		NEW DEFINITION FOR BETA
		double mew_ql = ((j+1)/(2*Nlinks))*j*mew_q;
		 //mew_ql = 0;
		double rl = (ejkl -mew_ql )/ var_q;
		// rl = ejkl;
		return rl;
	}
	
	
	public double calculate_localEjk(Node node)
	{
		if(node.getNumberOfLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		double N = network.getSize();
		double Nlinks = network.getNoOfLinks();
		double mew_q = qDist.getExpectation();
		double var_q = qDist.getVariance();
		//Now calculate the average remaining degree of all neighbours!
		
		
		
		int sumRemainingDegrees = 0;
		for(int y=0; y<node.getLinks().size();y++)
		{
			if(node.getLinks().elementAt(y).iNodeID() != node.getID())
			{
				//This is a neghbour
				sumRemainingDegrees = sumRemainingDegrees + node.getLinks().elementAt(y).iRealRemainingDegree();
			}
			
			if(node.getLinks().elementAt(y).jNodeID() != node.getID())
			{
				//This is a neighbour
				sumRemainingDegrees = sumRemainingDegrees + node.getLinks().elementAt(y).jRealRemainingDegree();
			}
		}
		
		double avgRD = ((double)sumRemainingDegrees )/ ((double)node.getNumberOfLinks());
		
	
		double j = node.getNumberOfLinks()-1;
		double ejkl = (j*(j+1)*avgRD)/(2*Nlinks);
		
		return ejkl;
	}
	
	public double calculate_r_as_Sum_of_Rl()
	{
		double r = 0;
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			r = r + calculate_localR(thisNode);
		}
		return r;
	}
	
	
	
	/**
	public double calculate_r_based_on_ejk_2()
	{
		double[][] ejk = network.getLinkDistribution2();
		double ejk_sum = 0;
		for(int y=0; y<ejk.length;y++)
		{
			for(int z=0; z<ejk.length;z++)
			{
				ejk_sum = ejk_sum + y*z*ejk[y][z];
			}
		}
			
		
		double mew_q = qDist2.getExpectation();
		double var_q = qDist2.getVariance();
		
		if(var_q==0)
		{
			if(mew_q >= 0)
				return 1.0;
			else
				return 0;
		}
		
		return ((ejk_sum - (mew_q*mew_q)) / var_q);
		
	
	}
	*/
	
	
	
	
	/**
	 * Calculate proper rl distribution
	 * Number of bins are 1/100 th of network size
	 * @return
	 */
	public double[] get_rl_dist_bin()
	{
		int bins = (int)(network.getSize() / 10.0);
		double rlmax = -1;
		double rlmin = 1;
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			double thisrl = calculate_localR(thisNode);
			if(thisrl> rlmax)
			{
				rlmax = thisrl;
			}
			if(thisrl< rlmin)
			{
				rlmin = thisrl;
			}
			
		}
		double binSize = rlmax / (double)bins;
		double[] rl = new double[bins];
		double times[] = new double[bins];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			double thisrl = calculate_localR(thisNode);
			int index = (int)(bins*(thisrl - rlmin) / (rlmax - rlmin));
			if(index>= bins)
			{
				index = bins-1; //End correction
			}
			rl[index] = rl[index]+1;
			
		}
		
		
		
		return rl;
	}
	
	public double[] get_rl_dist_xnodes()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[network.getAllNodes().size()];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			
			rl[y] = calculate_localR(thisNode); //Accumulated r for that degree
			/** THIS WAS ADDED TO TEST AN ALTERNATIVE DEFINITION
			if(thisNode.getNumberOfLinks()>0)
			{
		
				rl[y] = rl[y]/((double)thisNode.getNumberOfLinks()); //Degree
			}
			else
			{
				rl[y] = 0;
			}
			*/
			
		}
		
		
		
		return rl;
	}
	
	
	public double[][] get_rl_dist_xnodes_EXTENDED()
	{
		int Np = network.getDegreeDistribution().length;
		double[][] rl = new double[network.getAllNodes().size()][4];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[ind][0] = ind;
			rl[ind][1] = calculate_localR(thisNode); 
			rl[ind][2] = thisNode.getNumberOfLinks();
			rl[ind][3] = thisNode.getNumberOfLinks();
			
		}
		
		
		
		return rl;
	}
	
	
	
	public double[] get_avgND_dist_xnodes()
	{
		int Np = network.getDegreeDistribution().length;
		double[] nd = new double[network.getAllNodes().size()];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			
			nd[y] = calculate_averageND(thisNode); //Accumulated r for that degree
			
			
		}
		
		
		
		return nd;
	}


	public double[] get_ejkl_dist_xnodes()
	{
		int Np = network.getDegreeDistribution().length;
		double[] dist = new double[network.getAllNodes().size()];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			
			dist[y] = calculate_localEjk(thisNode); //Accumulated r for that degree
			//dist[y] = thisNode.getID(); //TODO: Hack - delete this!
		}
		
		
		
		return dist;
	}
	
	public double calculate_localRd(Node node)
	{
		return calculate_localR1(node) + calculate_localR2(node);
	}
	
	
	public double calculate_localR1(Node node)
	{
		if(node.getNumberOfOutLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		double N = network.getSize();
		double Nlinks = network.getNoOfLinks();
		double mew_qIn = qDistIn.getExpectation();
		double var_qIn = qDistIn.getVariance();
		
		double mew_qOut = qDistOut.getExpectation();
		double var_qOut = qDistOut.getVariance();
		
		//Now calculate the average remaining degree of all neighbours!
		
		if(var_qIn == 0 || var_qOut==0)
		{
			// Only one type of degree!
			return (1.0 / N);
		}
		
		int sumNeighbourInDegrees = 0;
		for(int y=0; y<node.getLinks().size();y++)
		{
			if(node.getLinks().elementAt(y).iNodeID() != node.getID())
			{
				//This is a neghbour from whom WE GET AN EDGE
				//sumNeighbourInDegrees = sumNeighbourInDegrees + node.getLinks().elementAt(y).iNode().getNumberOfInLinks();
			}
			
			if(node.getLinks().elementAt(y).jNodeID() != node.getID())
			{
				//This is a neighbour to whom we have an edge
				sumNeighbourInDegrees = sumNeighbourInDegrees + node.getLinks().elementAt(y).jNode().getNumberOfInLinks();
			}
		}
		
		double avrNIND = ((double)sumNeighbourInDegrees )/ ((double)node.getNumberOfOutLinks());
		
	
		double j = node.getNumberOfOutLinks();
		double ejkl = (j*(j)*avrNIND)/(1.0*Nlinks);
		double mew_ql = ((j)/(1.0*Nlinks))*mew_qIn*mew_qOut;
		double rl = (ejkl -mew_ql )/ (Math.sqrt(var_qIn)*Math.sqrt(var_qOut));
		
//		new rho 1
		double kv_out = node.getNumberOfOutLinks();
		double firstpart = (kv_out*kv_out)/(2.0*Nlinks);
		double secondpart = avrNIND - mew_qIn;
		rl = (firstpart*secondpart)/(Math.sqrt(var_qIn)*Math.sqrt(var_qOut));
		return rl;
	}
	
	
	public double calculate_localR2(Node node)
	{
		if(node.getNumberOfInLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		double N = network.getSize();
		double Nlinks = network.getNoOfLinks();
		double mew_qIn = qDistIn.getExpectation();
		double var_qIn = qDistIn.getVariance();
		
		double mew_qOut = qDistOut.getExpectation();
		double var_qOut = qDistOut.getVariance();
		
		//Now calculate the average remaining degree of all neighbours!
		
		if(var_qIn == 0 || var_qOut==0)
		{
			// Only one type of degree!
			return (1.0 / N);
		}
		
		int sumNeighbourOutDegrees = 0;
		for(int y=0; y<node.getLinks().size();y++)
		{
			if(node.getLinks().elementAt(y).iNodeID() != node.getID())
			{
				//This is a neghbour from whom WE GET AN EDGE
				sumNeighbourOutDegrees = sumNeighbourOutDegrees + node.getLinks().elementAt(y).iNode().getNumberOfOutLinks();
			}
			
			if(node.getLinks().elementAt(y).jNodeID() != node.getID())
			{
				//This is a neighbour to whom we have an edge
				//sumNeighbourOutDegrees = sumNeighbourOutDegrees + node.getLinks().elementAt(y).jNode().getNumberOfOutLinks();
			}
		}
		
		double avrNOUTD = ((double)sumNeighbourOutDegrees )/ ((double)node.getNumberOfInLinks());
		
	
		double j = node.getNumberOfInLinks();
		double ejkl = (j*(j)*avrNOUTD)/(1.0*Nlinks);
		double mew_ql = ((j)/(1.0*Nlinks))*mew_qIn*mew_qOut;
		double rl = (ejkl -mew_ql )/ (Math.sqrt(var_qIn)*Math.sqrt(var_qOut));
		
		//new rho 2
		double kv_in = node.getNumberOfInLinks();
		double firstpart = (kv_in*kv_in)/(2.0*Nlinks);
		double secondpart = avrNOUTD - mew_qOut;
		rl = (firstpart*secondpart)/ (Math.sqrt(var_qIn)*Math.sqrt(var_qOut));
	
		
		return rl;
	}
	
	public double calculate_localRout(Node node)
	{
		if(node.getNumberOfOutLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		double N = network.getSize();
		double Nlinks = network.getNoOfLinks();
		double mew_qOutDash = qDistOutDash.getExpectation();
		double var_qOutDash = qDistOutDash.getVariance();
		
		double mew_qOut = qDistOut.getExpectation();
		double var_qOut = qDistOut.getVariance();
		

		if(var_qOut < 1.0E-9){
			var_qOut = 0;
		}
		
		if(var_qOutDash < 1.0E-9){
			var_qOutDash = 0;
		}
		//Now calculate the average remaining degree of all neighbours!
		
		if(var_qOutDash == 0 || var_qOut==0)
		{
			System.out.println("AT LEAST ONE VARIANCE IS ZERO IN QOUT");
			// Only one type of in degree!
			if(var_qOut == 0 && var_qOutDash==0) //both are peaky
			{
				if(qDistOut.getDist()[0]== qDistOutDash.getDist()[0]) // peak at same point, completely assortative
				{
					return (node.getNumberOfOutLinks() / Nlinks);
				}
				else
				{
					//peak at different points, completely disassortative
					return -(node.getNumberOfOutLinks() / Nlinks);
				}
			}
			else //Only one of them is peaky, other is not:
				//ASSUME THAT all kind of connections are possible
			{
				System.out.println("spiky Vs spread q distributions");
				return 0;
			}
		}
		
		int sourceSumNeighbourOutDegrees = 0;
		int targetSumNeighbourOutDegrees = 0;
		
		for(int y=0; y<node.getLinks().size();y++)
		{
			if(node.getLinks().elementAt(y).iNodeID() != node.getID())
			{
				//This is a neghbour from whom WE GET AN EDGE
				sourceSumNeighbourOutDegrees = sourceSumNeighbourOutDegrees + node.getLinks().elementAt(y).iNode().getNumberOfOutLinks();
			}
			
			if(node.getLinks().elementAt(y).jNodeID() != node.getID())
			{
				//This is a neighbour to whom we have an edge
				targetSumNeighbourOutDegrees = targetSumNeighbourOutDegrees + node.getLinks().elementAt(y).jNode().getNumberOfOutLinks();
			}
		}
		
		double avrNOUTD_source = ((double)sourceSumNeighbourOutDegrees )/ ((double)node.getNumberOfInLinks());
		double avrNOUTD_target = ((double)targetSumNeighbourOutDegrees )/ ((double)node.getNumberOfOutLinks());
		
		if(node.getNumberOfInLinks()<=0)
		{
			avrNOUTD_source = 0;
		}
	
		if(node.getNumberOfOutLinks()<=0)
		{
			avrNOUTD_target = 0;
		}
	
		double j = node.getNumberOfInLinks();
		//double ejkl = (j*(j)*avrIND)/(1.0*Nlinks);
		//double mew_ql = ((j)/(1.0*Nlinks))*mew_qIn*mew_qInDash;
		//double rl = (ejkl -mew_ql )/ (Math.sqrt(var_qIn)*Math.sqrt(var_qInDash));
		double kv_in = node.getNumberOfInLinks();
		double kv_out = node.getNumberOfOutLinks();
		
		
		double firstTerm = kv_out/ (2*Nlinks*(Math.sqrt(var_qOut)*Math.sqrt(var_qOutDash)));
		double secondTerm = kv_out*(avrNOUTD_target- mew_qOutDash) + kv_in*(avrNOUTD_source- mew_qOut);
		double rl = firstTerm*secondTerm;
		return rl;

	
		//double j = node.getNumberOfOutLinks();
		//double ejkl = (j*(j)*avrNOUTD)/(1.0*Nlinks);
		//double mew_ql = ((j)/(1.0*Nlinks))*mew_qOutDash*mew_qOut;
		//double rl = (ejkl -mew_ql )/ (Math.sqrt(var_qOutDash)*Math.sqrt(var_qOut));
		//return rl;
	}
	
	public double calculate_localRin(Node node)
	{
		if(node.getNumberOfInLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		double N = network.getSize();
		double Nlinks = network.getNoOfLinks();
		double mew_qIn = qDistIn.getExpectation();
		double var_qIn = qDistIn.getVariance();
		
		double mew_qInDash = qDistInDash.getExpectation();
		double var_qInDash = qDistInDash.getVariance();
		

		
		if(var_qIn < 1.0E-9){
			var_qIn = 0;
		}
		
		if(var_qInDash < 1.0E-9){
			var_qInDash = 0;
		}
		
		if(var_qIn == 0 || var_qInDash==0)
		{
			System.out.println("AT LEAST ONE VARIANCE IS ZERO IN QIN");
			// Only one type of in degree!
			if(var_qIn == 0 && var_qInDash==0) //both are peaky
			{
				if(qDistIn.getDist()[0]== qDistInDash.getDist()[0]) // peak at same point, completely assortative
				{
					return (node.getNumberOfInLinks() / Nlinks);
				}
				else
				{
					//peak at different points, completely disassortative
					return -(node.getNumberOfInLinks() / Nlinks);
				}
			}
			else //Only one of them is peaky, other is not:
				//ASSUME THAT all kind of connections are possible
			{
				System.out.println("spiky Vs spread q distributions");
				return 0;
			}
		}
		
		int sourceSumNeighbourInDegrees = 0;
		int targetSumNeighbourInDegrees = 0;
		for(int y=0; y<node.getLinks().size();y++)
		{
			if(node.getLinks().elementAt(y).iNodeID() != node.getID())
			{
				//This is a neghbour from whom WE GET AN EDGE
				sourceSumNeighbourInDegrees = sourceSumNeighbourInDegrees + node.getLinks().elementAt(y).iNode().getNumberOfInLinks();
			}
			
			if(node.getLinks().elementAt(y).jNodeID() != node.getID())
			{
				//This is a neighbour to whom we have an edge
				targetSumNeighbourInDegrees = targetSumNeighbourInDegrees + node.getLinks().elementAt(y).jNode().getNumberOfInLinks();
			}
		}
		
		double avrIND = ((double)sourceSumNeighbourInDegrees )/ ((double)node.getNumberOfInLinks());
		if(node.getNumberOfInLinks()<=0)
		{
			avrIND = 0;
		}
		double avrINDTarget = ((double)targetSumNeighbourInDegrees )/ ((double)node.getNumberOfOutLinks());
		if(node.getNumberOfOutLinks()<=0)
		{
			avrINDTarget = 0;
		}
	
		double j = node.getNumberOfInLinks();
		double ejkl = (j*(j)*avrIND)/(1.0*Nlinks);
		double mew_ql = ((j)/(1.0*Nlinks))*mew_qIn*mew_qInDash;
		double rl = (ejkl -mew_ql )/ (Math.sqrt(var_qIn)*Math.sqrt(var_qInDash));
		double kv_in = node.getNumberOfInLinks();
		double kv_out = node.getNumberOfOutLinks();
		double firstTerm = kv_in/ (2*Nlinks*(Math.sqrt(var_qIn)*Math.sqrt(var_qInDash)));
		double secondTerm = kv_in*(avrINDTarget- mew_qInDash) + kv_out*(avrIND- mew_qIn);
		rl = firstTerm*secondTerm;
		return rl;
	}
	
	
	public double[] get_rl_dist()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfLinks()] = times[thisNode.getNumberOfLinks()]+1;
			rl[thisNode.getNumberOfLinks()] = rl[thisNode.getNumberOfLinks()]+calculate_localR(thisNode); //Accumulated r for that degree
			
		}
		
		for(int y=0; y<times.length;y++)
		{
			if(times[y] != 0)
			{
				rl[y] = rl[y] / times[y];
			}
			else
			{
				rl[y] = 0;
			}
			
		}
		
		return rl;
	}
	
	/**
	 * Directed Local assortativeness
	 * raw-out
	 * @return
	 */
	public double[] get_rl1_dist()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfOutLinks()] = times[thisNode.getNumberOfOutLinks()]+1;
			rl[thisNode.getNumberOfOutLinks()] = rl[thisNode.getNumberOfOutLinks()]+calculate_localR1(thisNode); //Accumulated r for that degree
			
		}
		
		for(int y=0; y<times.length;y++)
		{
			if(times[y] != 0)
			{
				rl[y] = rl[y] / times[y];
			}
			else
			{
				rl[y] = 0;
			}
			
		}
		
		return rl;
	}
	
	/**
	 * Directed Local assortativeness
	 * raw-in
	 * @return
	 */
	public double[] get_rl2_dist()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfInLinks()] = times[thisNode.getNumberOfInLinks()]+1;
			rl[thisNode.getNumberOfInLinks()] = rl[thisNode.getNumberOfInLinks()]+calculate_localR2(thisNode); //Accumulated r for that degree
			
		}
		
		for(int y=0; y<times.length;y++)
		{
			if(times[y] != 0)
			{
				rl[y] = rl[y] / times[y];
			}
			else
			{
				rl[y] = 0;
			}
			
		}
		
		return rl;
	}
	

	
	public double[] get_rld_dist()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfInLinks()] = times[thisNode.getNumberOfInLinks()]+1;
			rl[thisNode.getNumberOfInLinks()] = rl[thisNode.getNumberOfInLinks()]+calculate_localRd(thisNode); //Accumulated r for that degree
			
		}
		
		for(int y=0; y<times.length;y++)
		{
			if(times[y] != 0)
			{
				rl[y] = rl[y] / times[y];
			}
			else
			{
				rl[y] = 0;
			}
			
		}
		
		return rl;
	}
	/**
	 * Directed Local assortativeness
	 * raw-out
	 * @return
	 */
	public double[] get_rlOUT_dist()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfOutLinks()] = times[thisNode.getNumberOfOutLinks()]+1;
			rl[thisNode.getNumberOfOutLinks()] = rl[thisNode.getNumberOfOutLinks()]+calculate_localRout(thisNode); //Accumulated r for that degree
			
		}
		
		for(int y=0; y<times.length;y++)
		{
			if(times[y] != 0)
			{
				rl[y] = rl[y] / times[y];
			}
			else
			{
				rl[y] = 0;
			}
			
		}
		
		return rl;
	}
	
	public double[] get_rlIN_dist()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfInLinks()] = times[thisNode.getNumberOfInLinks()]+1;
			rl[thisNode.getNumberOfInLinks()] = rl[thisNode.getNumberOfInLinks()]+calculate_localRin(thisNode); //Accumulated r for that degree
			
		}
		
		for(int y=0; y<times.length;y++)
		{
			if(times[y] != 0)
			{
				rl[y] = rl[y] / times[y];
			}
			else
			{
				rl[y] = 0;
			}
			
		}
		
		return rl;
	}
	
	

}
