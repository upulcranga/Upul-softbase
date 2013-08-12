package oldClassesArchive;

import Basic.Link;
import Basic.Node;
import Basic.Network;
import Distributions.EDistribution;
import Distributions.QDistribution_In;
import Distributions.QDistribution_Out;
import Distributions.QDistribution_PowerLaw;

/**
 * This is a class to calculate assortativity and
 * node assortativity without involving distributions that can be memory hungry
 * @author piraveenan
 * @date 14.07.2012
 */
public class AssortativenessCalculatorLight {
	
	private Network network;

	boolean removeLinkless  =false;
	
	public AssortativenessCalculatorLight(Network net)
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
		
		if(this.get_variance_q()==0)
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
		double mew_q = this.get_mew_q();
		double var_q = this.get_variance_q();
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
			sum =  sum + (((double)(l.iNodeLinks()-0.0) - mew_q)*((double)(l.jNodeLinks() -0.0)- mew_q)); //sinec our calculations are in terms of 'remaining degrees'
			 
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
	
	
	
	
	
	
	
	
	public double get_variance_q()
	{
		double sigma = 0;
		double mew = 0;

		double M = network.getAllLinks().size();
		for(int y=0; y< M; y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			double x =  thisLink.iNodeLinks();
			double z =  thisLink.jNodeLinks();
			mew = mew + x +z;
			sigma = sigma + (x*x) + (z*z);
		}
		
		
	
			
		mew = mew / (2*M);
		
		sigma =  sigma / (2*M);
		
		sigma  =  sigma -  (mew*mew);
		
		return sigma;
	}
	
	public double get_mew_q()
	{
		
		double mew = 0;

		double M = network.getAllLinks().size();
		for(int y=0; y< M; y++)
		{
			Link thisLink = network.getAllLinks().elementAt(y);
			double x =  thisLink.iNodeLinks();
			double z =  thisLink.jNodeLinks();
			mew = mew + x +z;
		
		}
		
		
	
			
		mew = mew / (2*M);
		
	
		
		return mew;
	}
	

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




}
