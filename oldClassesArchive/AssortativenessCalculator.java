package oldClassesArchive;

import Basic.Link;
import Basic.Node;
import Basic.Network;
import Distributions.EDistribution;
import Distributions.QDistribution_In;
import Distributions.QDistribution_Out;
import Distributions.QDistribution_PowerLaw;

public class AssortativenessCalculator {
	
	private Network network;
	private QDistribution_PowerLaw qDist;
	private QDistribution_PowerLaw qDist2;
	private QDistribution_In qDistIn;
	private QDistribution_Out qDistOut;
	private QDistribution_In qDistInDash;
	private QDistribution_Out qDistOutDash;
	private EDistribution eDist;
	boolean removeLinkless  =false;
	
	public AssortativenessCalculator(Network net)
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
		qDistIn = new QDistribution_In(net.getInDegreeDistribution());
		qDistOut = new QDistribution_Out(net.getOutDegreeDistribution());
		qDistInDash = new QDistribution_In(net.getInDegreeDistribution_dash());
		qDistOutDash = new QDistribution_Out(net.getOutDegreeDistribution_dash());
		//qDistIn = new QDistribution_In(net.getInRemainingDegreeDistribution());
		//qDistOut = new QDistribution_Out(net.getOutRemainingDegreeDistribution());
		eDist = new EDistribution(net);
	}
	/* The logic of this function is contained in
	 * "Assortative Mixing in Networks"
	 * MEJ newman
	 * Formula 4
	 * However, here it is calculated for a generic quantity of the node,
	 * given by node 'reading'.
	 */
	public double calculate_nondegree_r()
	{
	
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
	
			Link link = network.getAllLinks().elementAt(i);
			Node iNode = link.iNode();
			Node jNode = link.jNode();
			
			sum1 = sum1 + (iNode.getReading()*jNode.getReading());
			sum2 = sum2 + ( 0.5*(iNode.getReading()+ jNode.getReading()) );
			sum3 = sum3 + ( 0.5*(/**/ (iNode.getReading()*iNode.getReading()) + (jNode.getReading()*jNode.getReading())/**/ )  );
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
				System.out.println("Uncalculable r condition: UNKNOWN CONFIGURATION");
				r = 1.0;
			}
		}
		else
		{
			r = ( (sum1/M) - Math.pow((sum2/M),2)  ) / ( (sum3/M) - Math.pow((sum2/M),2) );
		}
		
		return r;
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
	
	
	public double calculate_localR_sensor(Node node)
	{
		if(node.getNumberOfLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		double N = network.getSize();
		double Nlinks = network.getNoOfLinks();
		double mew_q = network.getqExpectation_sensor();
		double var_q =network.getqVariance_sensor();
		
		//Now calculate the average neighbour reading of all neighbours!
		
		if(var_q == 0)
		{
			// Only one type of degree!
			return (1.0 / N);
		}
		
		double sumNeighbourReading = 0;
		for(int y=0; y<node.getLinks().size();y++)
		{
			if(node.getLinks().elementAt(y).iNodeID() != node.getID())
			{
				//This is a neghbour
				sumNeighbourReading = sumNeighbourReading + node.getLinks().elementAt(y).iNode().getReading();
			}
			
			if(node.getLinks().elementAt(y).jNodeID() != node.getID())
			{
				//This is a neighbour
				sumNeighbourReading = sumNeighbourReading + node.getLinks().elementAt(y).jNode().getReading();
			}
			
			if(node.getLinks().elementAt(y).iNodeID() == node.getLinks().elementAt(y).jNodeID())
			{
				System.out.println("self link: shouldnt happen!");
				//This is a self-link. Add anyway.
				sumNeighbourReading = sumNeighbourReading + node.getLinks().elementAt(y).iNode().getReading();
			}
		}
		
		double avgNR = ((double)sumNeighbourReading )/ ((double)node.getNumberOfLinks());
		
	
		double d = node.getNumberOfLinks();
		double reading = node.getReading();
		
		double firstpart = reading*d/(2*Nlinks*var_q);
		double rl = firstpart*(avgNR -mew_q );
		return rl;
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
		//double mew_ql = ((j+1)/(2*Nlinks))*mew_q*mew_q; //TODO: Old def -  should not be active!!!
		//NEW DEFINITION FOR BETA
		double mew_ql = ((j+1)/(2*Nlinks))*j*mew_q;
		double rl = (ejkl -mew_ql )/ var_q;
		
		return rl;
	}
	
	/**
	 * Latest -  'new' definition for local assortativity
	 * called 'Node assortativity'
	 * Derived mathematically, from linked assortativity
	 * @param node
	 * @return
	 */
	public double calculate_nodeR_undirected(Node node)
	{
		if(node.getNumberOfLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		
		double N = network.getSize();
		double M = network.getNoOfLinks();
		double mew_q = qDist.getExpectation();
		double var_q = qDist.getVariance();
		double dv = node.getNumberOfLinks(); //Note: this is degree, not number of neighbours, which might be less
		//Now calculate the average remaining degree of all neighbours!
		double Mnode = node.getLinks().size();
		
		if(var_q == 0)
		{
			// Only one type of degree!
			return (1.0 / N);
		}
		
		double sum = 0;
		for(int y=0; y<Mnode;y++)
		{
			Link l = node.getLinks().elementAt(y);
			long  xn = l.iNodeID();
			long yn = l.jNodeID();
			
			 if(xn==yn)
			 {
				 System.out.println("self link in network");
			 }
			sum =  sum + ((l.iNodeLinks()-1 - mew_q)*(l.jNodeLinks() -1- mew_q)); //sinec our calculations are in terms of 'remaining degrees'
			 
		}
		

		double rl1 = sum;
		double rl =  rl1/(2.0*M*var_q);
		
		if(Math.abs(rl) >1)
		{
			//shouldn't be: for debugging
			rl =  rl*1;
		}
		
		return rl;
	}
	
	
	/**
	 * Latest -  'new' 'traversing' definition for local assortativity
	 * Also called 'Node assortativity'
	 * Derived by traversal and allocation, as opposed to the analytical derivation
	 * @param node
	 * @return
	 */
	public double calculate_deltaR_undirected(Node node, double r, double S)
	{
		if(node.getNumberOfLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		
		double N = network.getSize();
	
		double lambda = (r+1)/N;
		
		//lambda=0;
		//S = 1;
		
		double delta = node.getAverageNeighbourDifference() /  S;
		
		double delta_hat = lambda -  delta;
		
		return delta_hat;
		
		
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
		
		
		double firstTerm = kv_out/ (2.0*Nlinks*(Math.sqrt(var_qOut)*Math.sqrt(var_qOutDash)));
		double secondTerm = kv_out*(avrNOUTD_target- mew_qOutDash) + kv_in*(avrNOUTD_source- mew_qOut);
		double rl = firstTerm*secondTerm;
		//return rl;

	
		//double j = node.getNumberOfOutLinks();
		//double ejkl = (j*(j)*avrNOUTD)/(1.0*Nlinks);
		//double mew_ql = ((j)/(1.0*Nlinks))*mew_qOutDash*mew_qOut;
		//double rl = (ejkl -mew_ql )/ (Math.sqrt(var_qOutDash)*Math.sqrt(var_qOut));
		//return rl;
		
		
		double alpha1 = kv_out*kv_out*avrNOUTD_target / (1.0*Nlinks);
		double beta = (kv_out / ( 2.0*Nlinks))* (kv_out*mew_qOutDash + kv_in*mew_qOut);
		//rl = (alpha1 - beta)/(Math.sqrt(var_qOut)*Math.sqrt(var_qOutDash));
		return rl;
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
		double firstTerm = kv_in/ (2.0*Nlinks*(Math.sqrt(var_qIn)*Math.sqrt(var_qInDash)));
		double secondTerm = kv_in*(avrIND- mew_qInDash) + kv_out*(avrINDTarget- mew_qIn);
		rl = firstTerm*secondTerm;
		
		double alpha1 = kv_in*kv_in*avrIND / (1.0*Nlinks);
		double beta = (kv_in / ( 2.0*Nlinks))* (kv_in*mew_qInDash + kv_out*mew_qIn);
		//rl = (alpha1 - beta)/(Math.sqrt(var_qIn)*Math.sqrt(var_qInDash));
		return rl;
	}
	
	
	public double calculate_local_r_out_e(Link link)
	{
		
		double N = network.getSize();
		double E = network.getNoOfLinks();
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
			// Only one type of out degree!
			if(var_qOut == 0 && var_qOutDash==0) //both are peaky
			{
				if(qDistOut.getDist()[0]== qDistOutDash.getDist()[0]) // peak at same point, completely assortative
				{
					return ( 1.0/ E);
				}
				else
				{
					//peak at different points, completely disassortative
					return -(1.0 / E);
				}
			}
			else //Only one of them is peaky, other is not:
				//ASSUME THAT all kind of connections are possible
			{
				System.out.println("spiky Vs spread q distributions");
				return 0;
			}
		}
		
		Node sourcenode  = link.iNode();
		Node targetnode = link.jNode();
		
		double rho = (sourcenode.getNumberOfOutLinks() - mew_qOut)*(targetnode.getNumberOfOutLinks() - mew_qOutDash);
		rho =  rho /  (Math.sqrt(var_qOutDash)*Math.sqrt(var_qOut)*E);
		return rho;
	}
	
	public double calculate_local_r_in_e(Link link)
	{
		
		double N = network.getSize();
		double E = network.getNoOfLinks();
		double mew_qInDash = qDistInDash.getExpectation();
		double var_qInDash = qDistInDash.getVariance();
		
		double mew_qIn= qDistIn.getExpectation();
		double var_qIn = qDistIn.getVariance();
		

		if(var_qIn < 1.0E-9){
			var_qIn = 0;
		}
		
		if(var_qInDash < 1.0E-9){
			var_qInDash = 0;
		}
		//Now calculate the average remaining degree of all neighbours!
		
		if(var_qInDash == 0 || var_qIn==0)
		{
			System.out.println("AT LEAST ONE VARIANCE IS ZERO IN QIn");
			// Only one type of In degree!
			if(var_qIn == 0 && var_qInDash==0) //both are peaky
			{
				if(qDistIn.getDist()[0]== qDistInDash.getDist()[0]) // peak at same point, completely assortative
				{
					return ( 1.0/ E);
				}
				else
				{
					//peak at different points, completely disassortative
					return -(1.0 / E);
				}
			}
			else //Only one of them is peaky, other is not:
				//ASSUME THAT all kind of connections are possible
			{
				System.out.println("spiky Vs spread q distributions");
				return 0;
			}
		}
		
		Node sourcenode  = link.iNode();
		Node targetnode = link.jNode();
		
		double rho = (sourcenode.getNumberOfInLinks() - mew_qInDash)*(targetnode.getNumberOfInLinks() - mew_qIn);
		rho =  rho /  (Math.sqrt(var_qInDash)*Math.sqrt(var_qIn)*E);
		return rho;
	}
	
	/**
	 * This method is implemented to calculate edge-based congruity
	 * @date 14.04.2012
	 * @param link
	 * @return congruity of link
	 */
	public double calculate_congruity_in_e_undirected(Link link)
	{
		
		double N = network.getSize();
		double E = network.getNoOfLinks();
		QDistribution_PowerLaw sensor_qDist = new QDistribution_PowerLaw(network.getqkDistribution_sensor());
		double mew_q = sensor_qDist.getExpectation();
		double var_q = sensor_qDist.getVariance();
		

		if(var_q < 1.0E-9){
			var_q = 0;
		}
		
	
		//Now calculate the average remaining degree of all neighbours!
		
		if(var_q==0)
		{
			System.out.println(" VARIANCE IS ZERO IN Q sensor");
		
			return ( 1.0/ E);
				
		}
		
		Node sourcenode  = link.iNode();
		Node targetnode = link.jNode();
		
		double rho = (sourcenode.getReading() - mew_q)*(targetnode.getReading() - mew_q);
		rho =  rho /  ((var_q)*E);
		return rho;
	}
	
	public double get_mew_q()
	{
		return qDist.getExpectation();
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
	
	public double calculate_r_as_Sum_of_Rl1()
	{
		double r = 0;
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			r = r + calculate_localR1(thisNode);
		}
		return r;
	}
	
	public double calculate_r_as_Sum_of_Rl2()
	{
		double r = 0;
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			r = r + calculate_localR2(thisNode);
		}
		return r;
	}
	
	public double calculate_rOUT_as_Sum_of_RLOUT()
	{
		double r = 0;
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			r = r + calculate_localRout(thisNode);
		}
		return r;
	}
	
	public double calculate_rIN_as_Sum_of_RLIN()
	{
		double r = 0;
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			r = r + calculate_localRin(thisNode);
		}
		return r;
	}
	
	public double calculate_r_based_on_ejk()
	{
		double[][] ejk = network.getLinkDistribution();
		double ejk_sum = 0;
		for(int y=0; y<ejk.length;y++)
		{
			for(int z=0; z<ejk.length;z++)
			{
				ejk_sum = ejk_sum + y*z*ejk[y][z];
			}
		}
			
		
		double mew_q = qDist.getExpectation();
		double var_q = qDist.getVariance();
		
		if(var_q==0)
		{
			if(mew_q >= 0)
				return 1.0;
			else
				return 0;
		}
		
		return ((ejk_sum - (mew_q*mew_q)) / var_q);
		
	
	}
	
	//UNFINISHED
	public double calculate_r_based_on_links()
	{
		double r = 0;
		double mew = 0;
		double sigma = 0;
		double M = network.getAllLinks().size();
		for(int y=0; y< M; y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			double x =  thisLink.iNodeLinks();
			double z =  thisLink.jNodeLinks();
			mew = mew + x +z;
			//sigma = sigma + ((x*x)/(2*M)) - ((x*x)/(4*M*M)) + ((y*y)/(2*M)) - ((y*y)/(4*M*M));
			sigma = sigma + (x*x) + (z*z);
		}
		
		
	
			
		mew = mew / (2*M);
		
		sigma =  sigma / (2*M);
		
		sigma  =  sigma -  (mew*mew);
		
		for(int y=0; y< M; y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			r = r +((thisLink.iNodeLinks() -mew)*( thisLink.jNodeLinks()-mew));
		}

		r = r /(M*sigma); //Note well that what is calculated as 'sigma' is, in fact, variance.
		return r;
		
	
	}
	
	public double calculate_theoritical_r_sensornet()
	{
		double[][] ejk = network.getPjkDistribution();
		double ejk_sum = 0;
		for(int y=0; y<ejk.length;y++)
		{
			for(int z=0; z<ejk.length;z++)
			{
				ejk_sum = ejk_sum + y*z*ejk[y][z];
			}
		}
			
		QDistribution_PowerLaw sensor_qDist = new QDistribution_PowerLaw(network.getqkDistribution_sensor());
		double mew_q = sensor_qDist.getExpectation();
		double var_q = sensor_qDist.getVariance();
		
		if(var_q==0)
		{
			if(mew_q >= 0)
				return 1.0;
			else
				return 0;
		}
		
		return ((ejk_sum - (mew_q*mew_q)) / var_q);
		
	
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
	 * This method calculated r for a directed networks based on theoritical
	 * formulae, ejk and in, out degree distributions
	 * @return
	 */
	public double calculate_theoritical_directed_r()
	{
		double[][] ejk = network.getDirectedLinkDistribution();
		double ejk_sum = 0;
		for(int y=0; y<ejk.length;y++)
		{
			for(int z=0; z<ejk.length;z++)
			{
				ejk_sum = ejk_sum + y*z*ejk[y][z];
			}
		}
			
		
		double mew_qIn = qDistIn.getExpectation();
		double var_qIn = qDistIn.getVariance();
		
		double mew_qOut = qDistOut.getExpectation();
		double var_qOut = qDistOut.getVariance();
		
		if(var_qIn==0 ||var_qOut ==0 )
		{
			if(mew_qIn >= 0 && mew_qOut>=0)
				return 1.0;
			else
				return 0;
		}
		
		return ((ejk_sum - (mew_qIn*mew_qOut)) / (Math.sqrt(var_qIn)*Math.sqrt(var_qOut)));
		
	
	}
	
	/**
	 * The following two methods calculate assortativeness in directed networks 
	 * according to alternative definitions to that of Newman
	 * @return
	 */
	public double calculate_theoritical_directed_r_OO()
	{
		double[][] ejk = network.getDirectedLinkOutOutDistribution();
		double ejk_sum = 0;
		for(int y=0; y<ejk.length;y++)
		{
			for(int z=0; z<ejk.length;z++)
			{
				ejk_sum = ejk_sum + y*z*ejk[y][z];
			}
		}
			
		
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
		
		if(var_qOut == 0 || var_qOutDash==0)
		{
			System.out.println("AT LEAST ONE VARIANCE IS ZERO IN QOUT");
			if(var_qOut == 0 && var_qOutDash==0) //both are peaky
			{
				if(qDistOut.getDist()[0]== qDistOutDash.getDist()[0]) // peak at same point, completely assortative
				{
					return 1.0;
				}
				else
				{
					//peak at different points, completely disassortative
					return -1.0;
				}
			}
			else //Only one of them is peaky, other is not:
				//ASSUME THAT all kind of connections are possible
			{
				System.out.println("spiky Vs spread q out distributions");
				return 0;
			}
		}
		double f = ejk_sum - (mew_qOut*mew_qOutDash);
		
		double s=Math.sqrt(var_qOut)*Math.sqrt(var_qOutDash);
		 return (f/s);
		 
		 
		
		
	
	}
	
	
	
	public double calculate_theoritical_directed_r_InIn()
	{
		double[][] ejk = network.getDirectedLinkInInDistribution();
		double ejk_sum = 0;
		for(int y=0; y<ejk.length;y++)
		{
			for(int z=0; z<ejk.length;z++)
			{
				ejk_sum = ejk_sum + y*z*ejk[y][z];
			}
		}
			
		
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
					return 1.0;
				}
				else
				{
					//peak at different points, completely disassortative
					return -1.0;
				}
			}
			else //Only one of them is peaky, other is not:
				//ASSUME THAT all kind of connections are possible
			{
				System.out.println("spiky Vs spread q in distributions");
				return 0;
			}
		}
		
		
		double f = ejk_sum - (mew_qIn*mew_qInDash);
		
		double s = Math.sqrt(var_qIn)*Math.sqrt(var_qInDash);
		
		return f / s;
		
		//return ((ejk_sum - (mew_qIn*mew_qInDash)) / (Math.sqrt(var_qIn))*Math.sqrt(var_qInDash));
		
	
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
	
	public double[] get_rl_dist_sensor()
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
	
	//To calculate the number of positive links based on assortativity
	public int getPositiveLinks()
	{
	
		double[] rl = get_rl_e_dist_sensor();
		int pos=0;
		
		for(int y=0; y<network.getAllLinks().size();y++)
		{
			if(rl[y]>0)
				pos++;
		}
		
		return pos;
		
		
		
		
	}
	
	public double[][] get_rl_e_dist_sensor_MULTI()
	{
		int Np = network.getDegreeDistribution().length;
		//double[] rl = new double[2*Np];
		double[][] rl = new double[network.getNoOfLinks()][3];
		double times[] = new double[2*Np];
		for(int y=0; y<network.getAllLinks().size();y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			times[thisLink.iNode().getNumberOfLinks()+thisLink.jNode().getNumberOfLinks()] = times[thisLink.iNode().getNumberOfLinks()+thisLink.jNode().getNumberOfLinks()]+1;
			//rl[thisLink.iNode().getNumberOfLinks()+thisLink.jNode().getNumberOfLinks()] = rl[thisLink.iNode().getNumberOfLinks()+thisLink.jNode().getNumberOfLinks()]+calculate_congruity_in_e_undirected(thisLink); //Accumulated r for that degree
			rl[y][0] = calculate_congruity_in_e_undirected(thisLink);
			//rl[y][1] = thisLink.iNode().getReading() + thisLink.jNode().getReading();
			//rl[y][2] = thisLink.iNode().getPageRank() + thisLink.jNode().getPageRank();
			
			rl[y][1] = Math.abs(thisLink.iNode().getReading() - thisLink.jNode().getReading());
			rl[y][2] = Math.abs(thisLink.iNode().getPageRank() - thisLink.jNode().getPageRank());
		}
		
		for(int y=0; y<times.length;y++)
		{
			if(times[y] != 0)
			{
				//rl[y] = rl[y] / times[y];
			}
			else
			{
			//	rl[y] = 0;
			}
			
		}
		
		return rl;
	}
	
	public double[] get_rl_e_dist_sensor()
	{
		int Np = network.getDegreeDistribution().length;
		//double[] rl = new double[2*Np];
		double[] rl = new double[network.getNoOfLinks()];
		double times[] = new double[2*Np];
		for(int y=0; y<network.getAllLinks().size();y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			times[thisLink.iNode().getNumberOfLinks()+thisLink.jNode().getNumberOfLinks()] = times[thisLink.iNode().getNumberOfLinks()+thisLink.jNode().getNumberOfLinks()]+1;
			//rl[thisLink.iNode().getNumberOfLinks()+thisLink.jNode().getNumberOfLinks()] = rl[thisLink.iNode().getNumberOfLinks()+thisLink.jNode().getNumberOfLinks()]+calculate_congruity_in_e_undirected(thisLink); //Accumulated r for that degree
			rl[y] = rl[y] + calculate_congruity_in_e_undirected(thisLink);
		}
		
		for(int y=0; y<times.length;y++)
		{
			if(times[y] != 0)
			{
				//rl[y] = rl[y] / times[y];
			}
			else
			{
			//	rl[y] = 0;
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
	
	
	
	/**
	 * This is the new method which calculates the node assortativity distribution
	 * @return
	 */
	public double[] get_rnode_dist_undirected()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfLinks()] = times[thisNode.getNumberOfLinks()]+1;
			rl[thisNode.getNumberOfLinks()] = rl[thisNode.getNumberOfLinks()]+calculate_nodeR_undirected(thisNode); //Accumulated r for that degree
			
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
	 * This is the newest method which calculates the delta assortativity distribution
	 * @return
	 */
	public double[] get_rdelta_dist_undirected()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		//double r =  this.calculate_r_based_on_ejk();
		double r =  this.calculate_r_based_on_links();
		double S = network.sumOfDDAverage();
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfLinks()] = times[thisNode.getNumberOfLinks()]+1;
			rl[thisNode.getNumberOfLinks()] = rl[thisNode.getNumberOfLinks()]+calculate_deltaR_undirected(thisNode,r,S); //Accumulated r for that degree
			
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
	 * This is the new method which calculates the node assortativity distribution
	 * with single node IDs at the X axis
	 * @return
	 */
	public double[] get_rnode_dist_undirected_xnode()
	{
		double sum = 0;
		int Np = network.getSize();
		double[] rl = new double[Np];

		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			rl[y] = calculate_nodeR_undirected(thisNode); 
			sum = sum + rl[y];
		}
		
		
		System.out.println("sum of LA "+sum);
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
	
	
	/**
	 * This method is used to calculate edge based local assortativity
	 * Since nodes have no 'degree', the distribution is simply a random 'edge ID' vs edge's local assortativity
	 * The standard deviation etc of this distribution could be calculated.
	 * Also, assortativity of subnetworks could be calculated.
	 */
	public double[] get_rho_e_out_dist_xlinks()
	{

		double[] rl = new double[network.getAllLinks().size()];
	
		for(int y=0; y<network.getAllLinks().size();y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			
			rl[y] = calculate_local_r_out_e(thisLink); 
			
		}
		
		
		
		return rl;
	}
	
	public double[] get_rho_e_in_dist_xlinks()
	{

		double[] rl = new double[network.getAllLinks().size()];
	
		for(int y=0; y<network.getAllLinks().size();y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			
			rl[y] = calculate_local_r_in_e(thisLink); 
			
		}
		
		
		return rl;
	}
	
	
	/**
	 * This method is used to calculate edge based local assortativity
	 * Since nodes have no 'degree', the distribution is simply a random 'edge ID' vs edge's local assortativity
	 * The standard deviation etc of this distribution could be calculated.
	 * Also, assortativity of subnetworks could be calculated.
	 */
	public double[] get_rho_e_out_dist()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
	
		for(int y=0; y<network.getAllLinks().size();y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			int tot = Math.abs(thisLink.iNode().getNumberOfOutLinks() - thisLink.jNode().getNumberOfOutLinks());
			times[tot] = times[tot]+1;
			rl[tot] = rl[tot] +  calculate_local_r_out_e(thisLink); 
			
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
	
	
	public double[] get_rho_e_in_dist()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
	
		for(int y=0; y<network.getAllLinks().size();y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			int tot = thisLink.iNode().getNumberOfInLinks() + thisLink.jNode().getNumberOfInLinks();
			times[tot] = times[tot]+1;
			rl[tot] = rl[tot] +  calculate_local_r_in_e(thisLink); 
			
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
	
	
	public double[] get_rl_sensor_dist()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfLinks()] = times[thisNode.getNumberOfLinks()]+1;
			//times[(int)thisNode.getReading()] = times[(int)thisNode.getReading()]+1;
			rl[thisNode.getNumberOfLinks()] = rl[thisNode.getNumberOfLinks()]+calculate_localR_sensor(thisNode); //Accumulated r for that degree
			//rl[(int)thisNode.getReading()] = rl[(int)thisNode.getReading()]+calculate_localR_sensor(thisNode);
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
	
	public double[] get_rl_sensor_dist_xnodes()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[network.getSize()];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			//times[thisNode.getNumberOfLinks()] = times[thisNode.getNumberOfLinks()]+1;
			//times[(int)thisNode.getReading()] = times[(int)thisNode.getReading()]+1;
			rl[(int)thisNode.getID()] = calculate_localR_sensor(thisNode); //Accumulated r for that degree
			//rl[(int)thisNode.getReading()] = rl[(int)thisNode.getReading()]+calculate_localR_sensor(thisNode);
		}
		
		for(int y=0; y<times.length;y++)
		{
			if(times[y] != 0)
			{
				//rl[y] = rl[y] / times[y];
			}
			else
			{
				//rl[y] = 0;
			}
			
		}
		
		return rl;
	}
	
	
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
	
	public double[] get_rl1_dist_xnodes()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[network.getAllNodes().size()];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[ind] = calculate_localR1(thisNode); 
			
		}
		
		
		
		return rl;
	}
	
	public double[] get_rl2_dist_xnodes()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[network.getAllNodes().size()];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			rl[ind] = calculate_localR2(thisNode); 
			
			
		}
		
		
		
		return rl;
	}
	
	public double[] get_rlOUT_dist_xnodes()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[network.getAllNodes().size()];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[ind] = calculate_localRout(thisNode); 
			
		}
		
		
		
		return rl;
	}
	
	public double[] get_rlIN_dist_xnodes()
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[network.getAllNodes().size()];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[ind] = calculate_localRin(thisNode); 
			
		}
		
		
		
		return rl;
	}
	
	public double[][] get_rlIN_dist_xnodes_EXTENDED()
	{
		int Np = network.getDegreeDistribution().length;
		double[][] rl = new double[network.getAllNodes().size()][4];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[ind][0] = ind;
			rl[ind][1] = calculate_localRin(thisNode); 
			rl[ind][2] = thisNode.getNumberOfLinks();
			rl[ind][3] = thisNode.getNumberOfInLinks();
			
		}
		
		
		
		return rl;
	}
	
	/**
	 * This is the new method which calculates the node assortativity distribution
	 * with single node IDs at the X axis
	 * EXTENDED method to get more information
	 * @return
	 */
	public double[][] get_rnode_dist_undirected_xnode_EXTENDED()
	{
		double sum = 0;
		int Np = network.getSize();
		double[][] rl = new double[Np][4];

		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[y][0] = ind;
			rl[y][1] = calculate_nodeR_undirected(thisNode); 
			rl[y][2] = thisNode.getNumberOfLinks();
			rl[y][3] = rl[y][1] ; //Repeat LA -  for easiness of graph  plotting in Excel
			sum = sum + rl[y][1];
		}
		
		
		System.out.println("sum of LA "+sum);
		return rl;
	}
	

	/**
	 * This is the NEWEST method which calculates the node assortativity distribution  by traversals
	 * (calculates the delta_hats)
	 * with single node IDs at the X axis
	 * EXTENDED method to get more information
	 * @return
	 */
	public double[][] get_rdelta_dist_undirected_xnode_EXTENDED()
	{
		double sum = 0;
		int Np = network.getSize();
		double[][] rl = new double[Np][4];
		//double r = this.calculate_r_based_on_ejk();
		double r = this.calculate_r_based_on_links();
		double S = network.sumOfDDAverage();

		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[y][0] = ind;
			rl[y][1] = calculate_deltaR_undirected(thisNode, r, S); 
			rl[y][2] = thisNode.getNumberOfLinks();
			rl[y][3] = rl[y][1] ; //Repeat LA -  for easiness of graph  plotting in Excel
			sum = sum + rl[y][1];
		}
		
		
		System.out.println("sum of delta-hat NA "+sum);
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
	
	public double[][] get_rlOUT_dist_xnodes_EXTENDED()
	{
		int Np = network.getDegreeDistribution().length;
		double[][] rl = new double[network.getAllNodes().size()][4];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[ind][0] = ind;
			rl[ind][1] = calculate_localRout(thisNode); 
			rl[ind][2] = thisNode.getNumberOfLinks();
			rl[ind][3] = thisNode.getNumberOfOutLinks();
			
		}
		
		
		
		return rl;
	}
	
	public double[][] get_rl1_dist_xnodes_EXTENDED()
	{
		int Np = network.getDegreeDistribution().length;
		double[][] rl = new double[network.getAllNodes().size()][4];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[ind][0] = ind;
			rl[ind][1] = calculate_localR1(thisNode); 
			rl[ind][2] = thisNode.getNumberOfLinks();
			rl[ind][3] = thisNode.getNumberOfOutLinks();
			
		}
		
		
		
		return rl;
	}
	
	
	
	public double[][] get_rho_e_out_dist_xlinks_EXTENDED()
	{
		int Np = network.getDegreeDistribution().length;
		double[][] rl = new double[network.getAllLinks().size()][4];
		
		for(int y=0; y<network.getAllLinks().size();y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			int ind = y;
			
			rl[ind][0] = ind;
			rl[ind][1] = calculate_local_r_out_e(thisLink); 
			rl[ind][2] = thisLink.iNode().getNumberOfLinks();
			rl[ind][3] = thisLink.jNode().getNumberOfLinks();
			rl[ind][2] = thisLink.iNode().getID();
			rl[ind][3] = thisLink.jNode().getID();
			
		}
		
		
		
		return rl;
	}
	
	public double[][] get_rho_e_in_dist_xlinks_EXTENDED()
	{
		int Np = network.getDegreeDistribution().length;
		double[][] rl = new double[network.getAllLinks().size()][4];
		
		for(int y=0; y<network.getAllLinks().size();y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			int ind = y;
			
			rl[ind][0] = ind;
			rl[ind][1] = calculate_local_r_in_e(thisLink); 
			rl[ind][2] = thisLink.iNode().getNumberOfLinks();
			rl[ind][3] = thisLink.jNode().getNumberOfLinks();
			rl[ind][2] = thisLink.iNode().getID();
			rl[ind][3] = thisLink.jNode().getID();
			
		}
		
		
		
		return rl;
	}
	
	public double[][] get_rlD_dist_xnodes_EXTENDED()
	{
		int Np = network.getDegreeDistribution().length;
		double[][] rl = new double[network.getAllNodes().size()][4];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[ind][0] = ind;
			rl[ind][1] = calculate_localR1(thisNode) + calculate_localR1(thisNode); 
			rl[ind][2] = thisNode.getNumberOfLinks();
			rl[ind][3] = thisNode.getNumberOfOutLinks();
			
		}
		
		
		
		return rl;
	}
	
	public double[][] get_rl2_dist_xnodes_EXTENDED()
	{
		int Np = network.getDegreeDistribution().length;
		double[][] rl = new double[network.getAllNodes().size()][4];
		
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			int ind = (int)thisNode.getID();
			
			rl[ind][0] = ind;
			rl[ind][1] = calculate_localR2(thisNode); 
			rl[ind][2] = thisNode.getNumberOfLinks();
			rl[ind][3] = thisNode.getNumberOfInLinks();
			
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


}
