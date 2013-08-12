package Basic;

import java.util.Random;
import java.util.Vector;

public class SIRNetwork extends Network {


/**
 * Note well that the synchronized, and not the transient, states are taken into account.
 * @return
 */
public int countInfectedNodes()
{
	int count = 0;
	for(int i=0; i<getSize(); i++)
	{
		if(      ((SIRNode)(getAllNodes().elementAt(i))).   getSynchronizedSIRState().equals("I")     )
		{
			count++;
		}
	}
	
	
	return count;
}

public int countRecoveredNonImmunizedNodes()
{
	int count = 0;
	for(int i=0; i<getSize(); i++)
	{
		if(      ((SIRNode)(getAllNodes().elementAt(i))).   getSynchronizedSIRState().equals("R")     )
		{
			count++;
		}
	}
	
	
	return count;
}

public int countImmunizedNodes()
{
	int count = 0;
	for(int i=0; i<getSize(); i++)
	{
		if(      ((SIRNode)(getAllNodes().elementAt(i))).   getSynchronizedSIRState().equals("V")     )
		{
			count++;
		}
	}
	
	
	return count;
}

public int countSusceptibleNodes()
{
	int count = 0;
	for(int i=0; i<getSize(); i++)
	{
		if(      ((SIRNode)(getAllNodes().elementAt(i))).   getSynchronizedSIRState().equals("S")     )
		{
			count++;
		}
	}
	
	
	return count;
}

public void SynchronizeAllSIRStates()
{
	for(int i=0; i<getSize(); i++)
	{
		 SIRNode thisNode = ((SIRNode)(getAllNodes().elementAt(i)));
		 thisNode.synchState = thisNode.getSIRState(); //Not ideal implementation. TODO: need to make safer.
	}
}


/**
 * This method should be invoked only if recovery means vaccination
 * All nodes with stare "R" are given state "V"
 *
 * needs to be invoked after synchronization has been done
 */
public void confirmVaccination()
{
	for(int i=0; i<getSize(); i++)
	{
		 SIRNode thisNode = ((SIRNode)(getAllNodes().elementAt(i)));
		 if(thisNode.getSynchronizedSIRState().equals("R"))
		 {
			 thisNode.setSIRState("V");
			 thisNode.synchState =  "V";
		 }
			 
	}
}

/**
 * Return a set of all susceptible nodes, based on synchronized states
 * @return
 */
public Vector getAllSusceptibleNodes()
{
	Vector<SIRNode> susNodes =  new Vector<SIRNode>();
	
	for(int i=0; i<getSize(); i++)
	{
		SIRNode n = (SIRNode)(getAllNodes().elementAt(i));
		if(      n.   getSynchronizedSIRState().equals("S")     )
		{
			susNodes.addElement(n);
		}
	}
	
	return susNodes;
}


/**
 * This method is used to 'force' a POSITIVE end of epidemic, by making all nodes which are in 'I' state to have 'R' state
 *
 */
public void forceRecoveryofAllNodes()
{
	for(int i=0; i<this.getAllNodes().size(); i++)
	{
		SIRNode n = (SIRNode)(getAllNodes().elementAt(i));
		if(      n.   getSynchronizedSIRState().equals("I")     )
		{
			n.setSIRState("R");
		}
	}
	
	this.SynchronizeAllSIRStates();
}

/**
 * To make sure that the network is re-initialized for another simulation
 * of primary or secondary infection spread
 */
public void reInititializeAllNodes()
{
	for(int i=0; i<this.getAllNodes().size(); i++)
	{
		SIRNode n = (SIRNode)(getAllNodes().elementAt(i));
		
		n.setSIRState("S");
		
	}
	
	this.SynchronizeAllSIRStates();
}

/**
 * To make sure that the network is re-initialized for another simulation secondary infection spread
 * Networks that were vaccinated on primary simulation need to bepreserved
 */
public void reInititializeNonVaccinatedNodes()
{
	this.SynchronizeAllSIRStates();
	for(int i=0; i<this.getAllNodes().size(); i++)
	{
		SIRNode n = (SIRNode)(getAllNodes().elementAt(i));
		if(n.synchState.equals("V"))
		{
			//Do nothing
		}
		else
		{
			n.setSIRState("S");
		}
	}
	
	this.SynchronizeAllSIRStates();

}

/**
 * Grow this network according to BA preferential attachment, by n nodes and m links
 * The network is grown while retaining the states of existing nodes
 * All new nodes are given the state S
 * @param n increase size by this many nodes
 * @param m increase size by this many links
 * 
 * @return return the increased network
 */
public SIRNetwork growNetwork_PA(int n, int m)
{
	Random rand = new Random();

	double n0 =  this.getSize();
	double m0 =  this.getNoOfLinks(); //initial node and link sizes

	double p = m0/n0; //A ratio which indicates how many links per joining node should be created
	for(int i = 0; i < n; i++) //assuming the original network has IDs starting from 0
	{
		SIRNode node = new SIRNode((int)(n0+i));

		double probs[] =  new double[this.getSize()];
		
		for(int j = 0; j <this.getSize();j++)
		{
			Node destNode = this.getAllNodes().elementAt(j);
			double prob = p;
			if(this.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
			{
				 prob = p*((double)destNode.getNumberOfLinks()) / (double)this.getNoOfLinks();
				
			}
			else //Destination node has no links
			{
				prob = p;
				System.out.println("Destination has no links: Should MOT HAPPEN in growth scenario");
			}
			probs[j] = prob;
		}
		
		for(int j = 0; j <this.getSize();j++)
		{
			Node destNode = this.getAllNodes().elementAt(j);
			//Throw dice
			if(rand.nextDouble() <=  probs[j])
			{
				try 
				{
					this.addLink(node, destNode);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
			
		this.addNode(node);
	
	}

	return this;
}

/**
 * Grow this network according to ER random network model, by n nodes and m links
 * The network is grown while retaining the states of existing nodes
 * All new nodes are given the state S
 * @param n increase size by this many nodes
 * @param m increase size by this many links
 * 
 * @return return the increased network
 */
public SIRNetwork growNetwork_ER(int n, int m)
{
	Random rand = new Random();

	double n0 =  this.getSize();
	double m0 =  this.getNoOfLinks(); //initial node and link sizes

	for(int  i=0; i<n; i++) //assuming original network has ID starting with 0
	{
		SIRNode node = new SIRNode((int)(n0+i),(int)(n0+n) ); 
		this.addNode(node);
	}
	


	for(int  i=0; i<m0; i++)
	{
		int idx =  (Math.abs(rand.nextInt()) % this.getSize());
		int idy = (Math.abs(rand.nextInt()) % this.getSize());
		
		try {
			this.addLink(idx, idy);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}
	return this;
}

/**
* Grow this network according to SMALL WORLD  network model,proposed by Strgattz, by n nodes and m links
* The network is grown while retaining the states of existing nodes
 *@param K Number of nearest neighbours connected in the orginical ring -  Needs to be EVEN NUMBER
	 * @param p rewiring probability
	 * @return small world network
	 * Note well: Number of Links m = 0.5*N*K
	 * WARNING: returns a different network, yet messes up the original network also. Must be assigned to the object of itself
	 * Warning -  completely rewires - only the states of the original vaccinated nodes are maintained
*/
public SIRNetwork growNetwork_SW(int n, int K, double p)
{
	Random rand = new Random();

	double n0 =  this.getSize();
	double m0 =  this.getNoOfLinks(); //initial node and link sizes
	int N =  (int)(n0+n);
	
	int halfK = (int)(0.5*K);
	int M = halfK*N;
	int[][] links =  new int[N][N];

	for(int  i=0; i<n; i++) //assuming original network has ID starting with 0
	{
		SIRNode node = new SIRNode((int)(n0+i),(int)(n0+n) ); 
		this.addNode(node);
	}
	
	SIRNetwork newNet =  new SIRNetwork();
	
	for(int  i=0; i<this.getSize(); i++)
	{
		SIRNode node =  (SIRNode)(this.nodes.elementAt(i));
		SIRNode newNode =  node.clone();
		newNet.addNode(newNode);
	}

	
//	creat ring
	int linkCount = 0;
	for(int  i=0; i< N; i++)
	{
		int idx = i;
		for(int  j=1; j<= (halfK); j++)
		{
			int idy = (N+i -j)%N;
			int idz = (i + j)%N;
			links[idx][idy]=1;
			links[idx][idz]=1;
			linkCount = linkCount + 2;
		}
	}
	
	//Adjacency matrix has been created
	//travel in loop halfK times, stochastically destroy links clockwise, and add random links
	for(int  i=1; i<= halfK; i++)
	{
		for(int  j=0; j< N; j++)
		{
			if(rand.nextDouble() <=p)
			{
				links[j][(j+i)%N]=0;
				int randLink = (Math.abs(rand.nextInt())%N);
				links[j][randLink]=1;
			}
		}
	}
	
	//Now create links according to adjacency matrix
	
	for(int  x=0; x< N; x++)
	{
		for(int  y=(x+1); y< N; y++) //We do not need double links or self links
		{
			if(links[x][y]	>0)
			{
				try {
					newNet.addLink(x,y);
				} catch (Exception e) {
					System.out.println("SMALL WORLD NET: COULD NOT CREATE LINK BETWEEN "+x+ " and "+y);
					e.printStackTrace();
				}
			}
		}
	}
	
	System.out.println("SMALL WORLD NET CREATED WITH "+newNet.getSize()+" NODES AND "+newNet.getNoOfLinks()+ " LINKS");

	return newNet;

}


public SIRNetwork clone()
{
	 SIRNetwork copy = new SIRNetwork() ;
	 copy.networkSize = networkSize;
	 // Need to copy nodes, not give the same set to the cloned network
	 Vector<Node> n = new Vector<Node>(nodes.size());
	 {
		 for(int y=0; y<nodes.size();y++)
		{
			n.addElement(((Node)nodes.elementAt(y)).clone());
		}
	 }
	 copy.nodes = n;
	 copy.maxDegree = maxDegree;
	 // Rebuild the links with the cloned nodes */
	 Vector<Link> l = new Vector<Link>(this.links.size());
	 for(int y=0; y<links.size();y++)
	 {
		 Link thisLink = links.elementAt(y);
		 int id1 = (int)((Link)thisLink).iNodeID();
		 int id2 = (int)((Link)thisLink).jNodeID();
		 
		 try {
			Node n1 = copy.getNode(id1);
			Node n2 = copy.getNode(id2);
			 try
			 {
				copy.addLink((Node)n1, (Node)n2);
			 } catch (Exception e) {
				e.printStackTrace();
			 }
		} catch (Exception e1) {
			System.out.println("NODE WITH GIVEN ID does not exist in copy network");
			e1.printStackTrace();
		}
		
		 
		
		 //l.addElement(new DirectedLink(n.elementAt(id1-1), n.elementAt(id2-1)));
	 }
	 //copy.links = l;
	 copy.random = new Random();
	 return copy;
}


}