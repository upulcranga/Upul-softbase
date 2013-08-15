/**
 * @author Piraveenan
 * Key methods: addNode, addLink
 */

package Basic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.Random;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Stack;

import oldClassesArchive.AssortativenessCalculator;



import Distributions.EDistribution;
import Distributions.QDistribution_PowerLaw;

public class Network {

	protected Vector<Node> nodes;
	protected Vector<Link> links;


	protected int networkSize;
	protected int maxDegree;
	protected Random random;

	protected final int INF = Integer.MAX_VALUE;

	public Network()
	{
		networkSize = 0;
		nodes = new Vector<Node>();
		maxDegree = 1000;
		links = new Vector<Link>();
		random = new Random(System.currentTimeMillis());
	}

	public Network clone()
	{
		Network copy = new Network() ;
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

	/**
	 * This method casts the network into an SIR network, which is a subclass of networks
	 * All the nodes are 'converted' into SIRNodes. SIRNode is a subclass of Node.
	 * note that Nodes are not deep-cloned. That is, node attributes are not preserved, only the topology and node IDs are preserved,
	 */
	public SIRNetwork cloneToSIRNet()
	{
		SIRNetwork copy = new SIRNetwork() ;
		copy.networkSize = networkSize;


		Vector<Node> n = new Vector<Node>(nodes.size());
		{
			for(int y=0; y<nodes.size();y++)
			{
				Node thisNode = nodes.elementAt(y);
				n.addElement(thisNode.cloneToSIRNode());
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
				SIRNode n1 = (SIRNode)copy.getNode(id1);
				SIRNode n2 = (SIRNode)copy.getNode(id2);
				try
				{
					copy.addLink((SIRNode)n1, (SIRNode)n2);
				} catch (Exception e) {
					e.printStackTrace();
				}
			} catch (Exception e1) {
				System.out.println("NODE WITH GIVEN ID does not exist in copy network");
				e1.printStackTrace();
			}



		}

		copy.random = new Random(System.currentTimeMillis());
		return copy;
	}

	/* Note that this method only adds a node, does not link it to any other node.
	 * Therefore the Node is still isolated.
	 * It needs to be explicitly linked.
	 */
	public void addNode(Node node)
	{
		nodes.addElement(node);
		if(node.getAssignedDegree() >= maxDegree)
		{
			maxDegree = node.getAssignedDegree();
		}
		networkSize++;
	}

	public void selfLinkCheck()
	{
		for(int j=0;j< links.size();j++)
		{
			Link thisLink = links.elementAt(j);
			if(thisLink.iNodeID() == thisLink.jNodeID())
			{
				System.out.println("SELF LINK " + thisLink.iNodeID() + " "+ thisLink.jNodeID());
			}
		}
	}

	public int realMaxDegree()
	{
		int max = 0;
		for(int i=0; i< nodes.size(); i++)
		{
			int d =  nodes.elementAt(i).getNumberOfLinks();
			if(d>max)
			{
				max = d;
			}
		}

		return max;
	}

	public boolean removeNode(int nodeID) throws Exception
	{

		for(int i=0;i< networkSize;i++) 
		{
			Node node = nodes.elementAt(i);
			if(node.getID()== nodeID)
			{
				//Remove links
				for(int j=0;j< links.size();j++)
				{
					Link thisLink = links.elementAt(j);
					if(thisLink.iNodeID()==nodeID || thisLink.jNodeID()==nodeID)
					{
						//System.out.println("removed link when node destroyed "+thisLink.iNodeID()+ " "+thisLink.jNodeID());
						links.remove(j);
						j--; //new code - CHECK!!!
					}
				}

				//Remove Node
				nodes.remove(i);
				//System.out.println("NODE("+nodeID+") ACTUALLY REMOVED");
				networkSize--;
				return true;
			}
		}

		//If control comes here, node with the given ID not found
		System.out.println("NODE WITH GIVEN ID("+ nodeID+") NOT FOUND");
		return false;



	}



	public void removeAllLinks()
	{
		links.clear();
	}

	public void removeMarkedLinks()
	{
		int i = 0;
		while(true)
		{
			if(i >= links.size())
			{
				break;
			}
			Link  link = links.elementAt(i);
			if(link.isMarked())
			{
				links.remove(i);
				i--;
			}
			i++;
		}
		System.out.println("all marked links removed");

	}

	/*
	 * get all nodes with degree one; or with degrees less than a CUTOFF
	 */
	public Vector<Node> oneLinkNodes()
	{
		Vector<Node> olNodes =  new Vector();

		double CUTOFF = 2.0 + (0.01*nodes.size());
		for(int i=0; i< nodes.size(); i++)
		{
			Node thisNode =  nodes.elementAt(i);
			if(thisNode.getNumberOfLinks()<=CUTOFF)
			{
				olNodes.addElement(thisNode);
			}
		}

		return  olNodes;


	}

	/*
	 * This method deletes a marked Link and replaces it with two links.
	 * The two new links are given to the nodes that shared the deleted links,
	 * and these new nodes are connected now with two other Nodes based on Preferential attachment.
	 */

	public void replaceMarkedLinks()
	{
		Vector<Node> oln = oneLinkNodes();
		int olSize = oln.size();
		Random rand  =  new Random(System.currentTimeMillis());
		int i = 0;
		while(true)
		{
			if(i >= links.size())
			{
				break;
			}
			Link  link = links.elementAt(i);
			if(link.isMarked())
			{
				links.remove(i);
				i--;

				//create two new Links, with nodes selected preferentially based on their degrees
				int linkIDMax = getNoOfLinks();

				Node n1 = link.iNode();
				Node n2 = link.jNode();

				Node sNode; //Lets select the Node that is now linkless. If both Nodes have links, the second one is selected.
				if(n1.getNumberOfLinks()==0)
				{
					//System.out.println("LONELY NODE -  TO BE LINKED AGAIN");
					sNode = n1;
				}
				else if(n2.getNumberOfLinks()==0)
				{
					//System.out.println("LONELY NODE -  TO BE LINKED AGAIN");
					sNode = n2;
				}
				else
				{
					if(n1.getNumberOfLinks()> n2.getNumberOfLinks()) //select hub
					{
						sNode = n1;
					}
					else
					{
						sNode = n2;
					}
				}


				long peerID = (getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); //Select preferentially
				long peer2ID = (getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID();

				//long peerID = Math.abs((rand.nextInt()))%getSize(); //select randomly
				//long peer2ID = Math.abs((rand.nextInt()))%getSize(); 

				/* Select anti - preferentially */
				if(olSize>0)
				{
					long ind = 0;
					for(int y = nodes.size()-1; y >=0; y--){ //traverse the sorted nodes from peripherals to hubs
						double probl = (double)y / (0.5* nodes.size()*(nodes.size()+1)); //Higher indexed nodes -  i.e peripherals -  have better probability
						probl = (double)y /  nodes.size();
						// N/2(N+1) formula at the bottom
						if(rand.nextDouble() < probl) //throw coin
						{
							ind = ((Node)nodes.elementAt((int)y)).getID(); //Noet that the vector index is DIFFERENT from node ID
							break;
						}
						if(y==0)
						{
							ind = ((Node)nodes.elementAt((int)y)).getID(); //y = 0 now
							System.out.println("HUB CHOSEN");
						}
					}

					peerID = ind;

					//Do again
					ind = 0;
					for(int y = nodes.size()-1; y >=0; y--){ //traverse the sorted nodes from peripherals to hubs
						double probl = (double)y / (0.5* nodes.size()*(nodes.size()+1)); //Higher indexed nodes -  i.e peripherals -  have better probability
						probl = (double)y /  nodes.size();
						// N/2(N+1) formula at the bottom
						if(rand.nextDouble() < probl) //throw coin
						{
							ind = ((Node)nodes.elementAt((int)y)).getID();
							break;
						}
						if(y==0)
						{
							ind = ((Node)nodes.elementAt((int)y)).getID(); //y = 0 now
							System.out.println("HUB CHOSEN");
						}
					}


					peer2ID = ind;

					//long index = Math.abs((rand.nextInt()))%olSize; //select from single-degrees Nodes
					//peerID = oln.elementAt((int)index).getID();
					//index = Math.abs((rand.nextInt()))%olSize; //select from single-degrees Nodes
					//peer2ID = oln.elementAt((int)index).getID();
				}
				else
				{
					System.out.println("No one degreed Node exists");
				}


				try {

					Node peerNode = getNode((int)peerID);
					Node peer2Node = getNode((int)peer2ID);

					//sNode.addLink(peerNode); //No need to add links individually; the networks add link function takes care of that!
					//peerNode.addLink(sNode);
					//System.out.println("REPLACED WITH "+ peerNode.getNumberOfLinks()+ " "+sNode.getNumberOfLinks());
					addLink((int)peerID, (int)sNode.getID());
					i++;
					addLink((int)peer2ID, (int)sNode.getID());
					i++;

				} catch (Exception e) {
					System.out.println("problem REPLACING link");
					e.printStackTrace();
				}
			}
			i++;
		}
		//System.out.println("all marked links removed");

	}

	public Node getNode(int nodeID) throws Exception
	{
		Node node;
		if(nodeID>networkSize)
		{
			double z = 0;
			//throw new Exception("Node with Given ID Does not exist");
		}

		/*Search for the node with the given ID */

		for(int i=0;i< networkSize;i++) /* Note well: Node ID starts from one not zero */
		{
			node = nodes.elementAt(i);
			if(node.getID()== nodeID)
			{
				return node;
			}
		}

		/* If node has not been found during search, exit with exception */
		throw new Exception("Node with Given ID "+ nodeID+"  Does not exist");

	}

	/** Finds a given Nodes position in the nodes vector
	 * 
	 * 
	 * @param nodeID
	 * @return
	 */

	public int getPosition(int nodeID) 
	{
		Node node;
		if(nodeID>networkSize)
		{
			double z = 100000;
			//throw new Exception("Node with Given ID Does not exist");
		}

		/*Search for the node with the given ID */

		for(int i=0;i< networkSize;i++) /* Note well: Node ID starts from one not zero */
		{
			node = nodes.elementAt(i);
			if(node.getID()== nodeID)
			{
				return i;
			}
		}

		return 100000; //A BIG NUMBER



	}


	public void addLink(int sourceID, int destID) throws Exception
	{
		Node sourceNode = getNode(sourceID);
		Node destNode = getNode(destID);
		addLink(sourceNode, destNode);
	}

	public void addLink(Node sourceNode, Node destNode) throws Exception
	{
		sourceNode.addLink(destNode, true);
		if(sourceNode.getID()!= destNode.getID()) //Check for self links
		{
			destNode.addLink(sourceNode, false);
		}
		if(sourceNode.getAssignedDegree() >= maxDegree)
		{
			maxDegree = sourceNode.getAssignedDegree();
		}
		if(destNode.getAssignedDegree() >= maxDegree)
		{
			maxDegree = destNode.getAssignedDegree();
		}

		links.addElement(new Link(sourceNode, destNode));

	}

	public int getMaxDegree()
	{
		return maxDegree;
	}

	public int getMaxRemainingDegree()
	{
		return maxDegree-1;
	}


	public int getSize()
	{
		return networkSize;
	}

	public Vector<Node> getAllNodes()
	{
		return nodes;
	}

	public int getNoOfLinks()
	{
		return links.size();
	}
	public Vector<Link> getAllLinks()
	{
		return links;
	}

	/*
	 * Shuffle Nodes randomly
	 */
	public void shuffle()
	{
		for(int i=0;i< nodes.size();i++) 
		{
			int rnd = ( ((int)Math.abs(random.nextInt())) %nodes.size() );
			Node temp = nodes.elementAt(i); /* save the element at position i */
			nodes.setElementAt(nodes.elementAt(rnd),i); /* put the ndoe at random position in position i */
			nodes.setElementAt(temp,rnd); /* Put the noed that was in potision into random position */	
		}	
	}

	/* Move this function somewhere else */
	public Vector shuffle(Vector v)
	{

		for(int i=0;i< v.size();i++) 
		{
			int rnd = ( Math.abs(random.nextInt())%v.size() );
			Object temp = v.elementAt(i); /* save the element at position i */
			v.setElementAt(v.elementAt(rnd),i); /* put the ndoe at random position in position i */
			v.setElementAt(temp,rnd); /* Put the noed that was in potision into random position */	
		}

		return v;
	}

	/*
	 * This method returns a sorted list of links based on the sum of their degrees on each end.
	 * SMALLEST FIRST
	 */
	public void sortLinksBasedOnCombinedDegree()
	{
		System.out.println("SORT SIZE: " + getNoOfLinks());
		Vector<Link> sortedLinks = new Vector<Link>();

		Link[] allLinks = new Link[links.size()];
		links.copyInto(allLinks); //The links have been copied into the array

		Arrays.sort(allLinks); //Now the array has been sorted
		for(int i = 0; i< allLinks.length; i++)
		{
			sortedLinks.addElement(allLinks[i]);
			//System.out.println("Link with rd "+ allLinks[i].iNodeLinks()+ " "+ allLinks[i].jNodeLinks());
		}

		links =  sortedLinks; //now the actual link list is sorted  - No Copy

	}

	/**
	 * Reassigns the ID of a node and updates all links accordingly
	 * This method assumes that in the time of checking, the old ID of the node is not reassigend to any link
	 * @param newID
	 * @param node
	 */
	public void reassignNodeID(long newID, long oldID, Node node)
	{
		for(int i = 0; i< links.size(); i++)
		{
			Link thisLink =  links.elementAt(i);
			if(thisLink.iNodeID()== oldID)
			{
				System.out.println("Node " +newID + " remembered as "+ oldID);
			}

			if(thisLink.jNodeID()== oldID)
			{
				System.out.println("Node " +newID + " remembered as "+ oldID);
			}
		}

		for(int i = 0; i< nodes.size(); i++)
		{
			Node thisNode =  nodes.elementAt(i);
			for(int j = 0; j< thisNode.getLinks().size(); j++)
			{
				Link thisLink =  links.elementAt(j);
				if(thisLink.iNodeID()== oldID)
				{
					System.out.println("Node " +newID + " remembered as "+ oldID);
				}

				if(thisLink.jNodeID()== oldID)
				{
					System.out.println("Node " +newID + " remembered as "+ oldID);
				}
			}
		}

	}

	/**
	 * Sort Nodes based on degree
	 *
	 */
	public void sortNodes()
	{

		Vector<Node> sortedNodes = new Vector<Node>();

		Node[] allNodes = new Node[nodes.size()];
		nodes.copyInto(allNodes); 

		Arrays.sort(allNodes); //Now the array has been sorted
		for(int i = 0; i< allNodes.length; i++)
		{
			sortedNodes.addElement(allNodes[i]);
			//System.out.println("Link with rd "+ allLinks[i].iNodeLinks()+ " "+ allLinks[i].jNodeLinks());
		}

		nodes=  sortedNodes; //now the actual link list is sorted  - No Copy

	}

	public boolean isMember(Vector<Node> nodeSet, int nodeID)
	{
		boolean exists = false;
		for(int i=1;i<= nodeSet.size();i++) /* Note well: Node ID starts from one not zero */
		{
			if(nodeSet.elementAt(i).getID() ==nodeID)
			{
				exists = true;
			}
		}

		return exists;
	}

	/**
	 * Degree distribution
	 * @return
	 */
	public double[]  getDegreeDistribution()
	{
		double[] dist = new double[maxDegree+1];
		for(int i=0;i<nodes.size();i++) 
		{
			//dist[nodes.elementAt(i).getDegree()]= (1.0/nodes.size())+dist[nodes.elementAt(i).getDegree()];
			dist[nodes.elementAt(i).getNumberOfLinks()]= (1.0/nodes.size())+dist[nodes.elementAt(i).getNumberOfLinks()];
		}

		return dist;
	}

	public int getHighestDegree()
	{
		//double[] dist = getDegreeDistribution();
		int Np = 0;
		/*
		for(int i=0;i<dist.length;i++) 
		{
			if(dist[i] > 0 )
			{
				Np = i;
			}
		}*/

		for(int i=0;i<links.size();i++) 
		{
			if(links.elementAt(i).iNodeLinks() > Np )
			{
				Np = links.elementAt(i).iNodeLinks();
			}
		}

		return Np;
	}

	public int getHighestBCNodeID() 
	{

		double Np = 0;
		long id = (int) nodes.elementAt(0).getID();

		for(int i=0;i<nodes.size();i++) 
		{
			if(nodes.elementAt(i).getBCCentrality() >= Np )
			{
				Np = nodes.elementAt(i).getBCCentrality();
				id = nodes.elementAt(i).getID();
			}


		}


		return (int)id;
	}

	public int getHighestCCNodeID() 
	{

		double Np = 0;
		long id = (int) nodes.elementAt(0).getID();

		for(int i=0;i<nodes.size();i++) 
		{
			if(nodes.elementAt(i).getCCCentrality() >= Np )
			{
				Np = nodes.elementAt(i).getCCCentrality();
				id = nodes.elementAt(i).getID();
			}


		}


		return (int)id;
	}

	public int getHighestDegreeNodeID() throws Exception
	{

		int Np = 0;
		int id = getSize()+1;

		for(int i=0;i<links.size();i++) 
		{
			if(links.elementAt(i).iNodeLinks() >= Np )
			{
				Np = links.elementAt(i).iNodeLinks();
				id = (int) links.elementAt(i).iNode().getID();
			}

			if(links.elementAt(i).jNodeLinks() >= Np )
			{
				Np = links.elementAt(i).jNodeLinks();
				id = (int) links.elementAt(i).jNode().getID();
			}
		}

		if(id==(getSize() + 1))
		{
			//No node is found with maximum id
			//Cannot be
			System.out.println("No maximum exists");
			if(getSize()==0)
			{
				System.out.println("Empty network");
				throw new Exception("Empty network");

			}
			else
			{
				System.out.println("Complete fracturization of graph");
				id = (int) nodes.elementAt(0).getID();//just pick one!!
			}
		}
		return id;
	}

	/**
	 * In Degree distribution
	 * @return
	 */
	public double[]  getInRemainingDegreeDistribution()
	{
		double[] dist = new double[maxDegree];/* Length of remaining degree distribution is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{
			//the in degree should be calculated for j nodes, since they have the incoming links
			dist[links.elementAt(i).jNode().getNumberOfInLinks()-1]= (1.0/(double)(1.0*links.size()))+dist[links.elementAt(i).jNode().getNumberOfInLinks()-1];
		}

		return dist;
	}

	/**
	 * Degree distribution
	 * @return
	 */
	public double[]  getInDegreeDistribution()
	{
		double[] dist = new double[maxDegree];/* Length of remaining degree distribution is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{
			//the in degree should be calculated for j nodes, since they have the incoming links
			dist[links.elementAt(i).jNode().getNumberOfInLinks()]= (1.0/(double)(1.0*links.size()))+dist[links.elementAt(i).jNode().getNumberOfInLinks()];
		}

		return dist;
	}

	/**
	 * in Degree distribution OUT OF nodes
	 * @return
	 */
	public double[]  getInDegreeDistribution_dash()
	{
		double[] dist = new double[maxDegree];/* Length of remaining degree distribution is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{

			dist[links.elementAt(i).iNode().getNumberOfInLinks()]= (1.0/(double)(1.0*links.size()))+dist[links.elementAt(i).iNode().getNumberOfInLinks()];
		}

		return dist;
	}



	/**
	 * out Degree distribution (actually qout(k-1)
	 * @return
	 */
	public double[]  getOutDegreeDistribution()
	{
		double[] dist = new double[maxDegree];/* Length of remaining degree distribution is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{

			dist[links.elementAt(i).iNode().getNumberOfOutLinks()]= (1.0/(double)(1.0*links.size()))+dist[links.elementAt(i).iNode().getNumberOfOutLinks()];
		}

		return dist;
	}

	/**
	 * out Degree distribution INTO nodes
	 * @return
	 */
	public double[]  getOutDegreeDistribution_dash()
	{
		double[] dist = new double[maxDegree];/* Length of remaining degree distribution is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{

			dist[links.elementAt(i).jNode().getNumberOfOutLinks()]= (1.0/(double)(1.0*links.size()))+dist[links.elementAt(i).jNode().getNumberOfOutLinks()];
		}

		return dist;
	}

	/**
	 * out Degree distribution
	 * @return
	 */
	public double[]  getOutRemainingDegreeDistribution()
	{
		double[] dist = new double[maxDegree];/* Length of remaining degree distribution is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{
			//the in degree should be calculated for j nodes, since they have the incoming links
			dist[links.elementAt(i).iNode().getNumberOfOutLinks()-1]= (1.0/(double)(1.0*links.size()))+dist[links.elementAt(i).iNode().getNumberOfOutLinks()-1];
		}

		return dist;
	}


	/**
	 * Degree distribution
	 * @return
	 */
	public double[]  getRemainingDegreeDistribution()
	{
		double[] dist = new double[maxDegree];/* Length of remaining degree distribution is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{
			//dist[links.elementAt(i).iRemainingDegree()]= (1.0/(2*links.size()))+dist[links.elementAt(i).iRemainingDegree()];
			//dist[links.elementAt(i).jRemainingDegree()]= (1.0/(2*links.size()))+dist[links.elementAt(i).jRemainingDegree()];
			dist[links.elementAt(i).iNodeLinks()-1]= (1.0/(double)(2.0*links.size()))+dist[links.elementAt(i).iNodeLinks()-1];
			dist[links.elementAt(i).jNodeLinks()-1]= (1.0/(double)(2.0*links.size()))+dist[links.elementAt(i).jNodeLinks()-1];
		}

		return dist;
	}



	public double[][]  getLinkDistribution()
	{
		double[][] dist = new double[maxDegree][maxDegree];
		/* Length of Q(k), and therefore E(j,k), is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{			
			//dist[links.elementAt(i).iRemainingDegree()][links.elementAt(i).jRemainingDegree()]= (1.0/links.size())+dist[links.elementAt(i).iRemainingDegree()][links.elementAt(i).jRemainingDegree()];		
			/* we need to maintain the Eij = Eji symmetry  - TODO: Check the logic here */
			dist[links.elementAt(i).iNodeLinks()-1][links.elementAt(i).jNodeLinks()-1]= (0.5/links.size())+dist[links.elementAt(i).iNodeLinks()-1][links.elementAt(i).jNodeLinks()-1];	
			dist[links.elementAt(i).jNodeLinks()-1][links.elementAt(i).iNodeLinks()-1]= (0.5/links.size())+dist[links.elementAt(i).jNodeLinks()-1][links.elementAt(i).iNodeLinks()-1];

		}

		return dist;
	}

	/**
	 * For sensor networks
	 * In this case we assume that directedness does not matter
	 * @return
	 */
	public double[][]  getPjkDistribution()
	{
		double[][] dist = new double[maxDegree][maxDegree];
		/* Length of Q(k), and therefore E(j,k), is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{			

			/* we need to maintain the Eij = Eji symmetry  - TODO: Check the logic here */
			dist[(int)links.elementAt(i).iNode().getReading()][(int)links.elementAt(i).jNode().getReading()]= (0.5/links.size())+dist[(int)links.elementAt(i).iNode().getReading()][(int)links.elementAt(i).jNode().getReading()];
			dist[(int)links.elementAt(i).jNode().getReading()][(int)links.elementAt(i).iNode().getReading()]= (0.5/links.size())+dist[(int)links.elementAt(i).jNode().getReading()][(int)links.elementAt(i).iNode().getReading()];

		}

		return dist;
	}

	/**
	 * For sensor networks
	 * In this case we assume that directedness does not matter
	 * @return
	 */
	public double[]  getqkDistribution_sensor()
	{
		int maxReading = 1000;
		double[] dist = new double[maxReading];
		/* Length of Q(k), and therefore E(j,k), is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{			


			dist[(int)links.elementAt(i).iNode().getReading()]= (0.5/links.size())+dist[(int)links.elementAt(i).iNode().getReading()];
			dist[(int)links.elementAt(i).jNode().getReading()]= (0.5/links.size())+dist[(int)links.elementAt(i).jNode().getReading()];
		}

		return dist;
	}

	/**
	 * For sensor networks
	 * In this case we assume that directedness does not matter
	 * @return
	 */
	public double[]  getpkDistribution_sensor()
	{
		int maxReading = 1000;
		double[] dist = new double[maxReading];

		for(int i=0;i<nodes.size();i++) 
		{			


			dist[(int)nodes.elementAt(i).getReading()]= (1.0/nodes.size())+dist[(int)nodes.elementAt(i).getReading()];
		}

		return dist;
	}

	/**
	 * For sensor networks
	 * In this case we assume that directedness does not matter
	 * @return
	 */
	public double  getqExpectation_sensor()
	{
		double mew = 0;
		double[] dist = getqkDistribution_sensor();


		for(int i=0;i<dist.length;i++) 
		{
			mew = mew + (double)i*dist[i];
		}

		//mew =  mew / (double)dist.length;
		return mew;
	}


	/**
	 * For sensor networks
	 * In this case we assume that directedness does not matter
	 * @return
	 */
	public double  getqVariance_sensor()
	{
		double mew1 = 0;
		double mew2 = 0;
		double var = 0;
		double mew = 0;
		double[] dist = getqkDistribution_sensor();

		for(int k=0;k< dist.length; k++)
		{
			//System.out.println(k+"\t"+dist[k]);

		}

		int ha = 0;

		//System.exit(0);

		for(int i=0;i<dist.length;i++) 
		{
			mew1 = mew1 + (double)i*dist[i];
			mew2 = mew2 + (double)i*(double)i*(dist[i]);
		}


		//mew1 =  mew1 / (double)dist.length;
		//mew2 =  mew2 / (double)dist.length;
		double sigma = mew2 - (mew1*mew1);

		mew = getqExpectation_sensor();

		for(int i=0;i<dist.length;i++) 
		{
			var = var + (  ((double)i - mew)*((double)i - mew)*dist[i]   );

		}


		return var;
	}

	/**
	 * For sensor networks
	 * In this case we assume that directedness does not matter
	 * @return
	 */
	public double  getpVariance_sensor()
	{
		double mew1 = 0;
		double mew2 = 0;
		double[] dist = getpkDistribution_sensor();


		for(int i=0;i<dist.length;i++) 
		{
			mew1 = mew1 + (double)i*dist[i];
			mew2 = mew2 + (double)i*(double)i*(dist[i]);
		}



		double sigma = mew2 - (mew1*mew1);

		return sigma;
	}

	/**
	 * For sensor networks
	 * In this case we assume that directedness does not matter
	 * @return
	 */
	public double  getpMean_sensor()
	{
		double mew1 = 0;
		double mew2 = 0;
		double[] dist = getpkDistribution_sensor();


		for(int i=0;i<dist.length;i++) 
		{
			mew1 = mew1 + (double)i*dist[i];

		}





		return mew1;
	}




	public double[][]  getDirectedLinkDistribution()
	{
		double[][] dist = new double[maxDegree][maxDegree];
		/* Length of Q(k), and therefore E(j,k), is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{			
			//No symmetry is enforced here
			dist[links.elementAt(i).jNode().getNumberOfInLinks()][links.elementAt(i).iNode().getNumberOfOutLinks()]= (1.0/links.size())+dist[links.elementAt(i).jNode().getNumberOfInLinks()][links.elementAt(i).iNode().getNumberOfOutLinks()];		

		}

		return dist;
	}

	public double[][]  getDirectedLinkOutOutDistribution()
	{
		double[][] dist = new double[maxDegree][maxDegree];
		/* Length of Q(k), and therefore E(j,k), is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{			
			//No symmetry is enforced here
			dist[links.elementAt(i).jNode().getNumberOfOutLinks()][links.elementAt(i).iNode().getNumberOfOutLinks()]= (1.0/links.size())+dist[links.elementAt(i).jNode().getNumberOfOutLinks()][links.elementAt(i).iNode().getNumberOfOutLinks()];		

		}

		return dist;
	}

	public double[][]  getDirectedLinkInInDistribution()
	{
		double[][] dist = new double[maxDegree][maxDegree];
		/* Length of Q(k), and therefore E(j,k), is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		for(int i=0;i<links.size();i++) 
		{			
			//No symmetry is enforced here
			dist[links.elementAt(i).jNode().getNumberOfInLinks()][links.elementAt(i).iNode().getNumberOfInLinks()]= (1.0/links.size())+dist[links.elementAt(i).jNode().getNumberOfInLinks()][links.elementAt(i).iNode().getNumberOfInLinks()];		

		}

		return dist;
	}




	public double[][] getConditionalDegreeDistribution()
	{
		double[][] dist = new double[maxDegree+1][maxDegree+1];
		double[][] Ejk = getLinkDistribution();
		double[] Qk = getRemainingDegreeDistribution();
		for(int i=1;i<=maxDegree;i++) 
		{
			for(int j=1;j<=maxDegree;j++) 
			{
				if(Qk[i-1]!=0)
					dist[i][j] = (Ejk[i-1][j-1] / Qk[i-1]);
				else
					dist[i][j] = 0;
			}

		}
		return dist;
	}
	public void printArray(double[][] arr, String arrName)
	{
		System.out.print(arrName+" ");
		for(int y=0; y<arr.length;y++)
		{
			for(int j=0;j<arr[y].length;j++) 
			{
				System.out.println("("+y+","+j+") \t"+arr[y][j]+"");
			}
		}
		System.out.println("");
	}

	public void printArraytoFile(String[] arr, String fileName)
	{
		BufferedWriter out = null;


		try {

			out = new BufferedWriter(new FileWriter(fileName+".xls"));

		}

		catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}

		for(int y=0; y<arr.length;y++)
		{

			try {
				out.write(y+"\t"+arr[y]+"\t");
				out.flush();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			try {
				out.write(System.getProperty(("line.separator")));
				out.flush();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

	public void printArraytoFile(double[][] arr, String fileName)
	{
		BufferedWriter out = null;


		try {

			out = new BufferedWriter(new FileWriter(fileName+".xls"));

		}

		catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}

		for(int y=0; y<arr.length;y++)
		{
			for(int j=0;j<arr[y].length;j++) 
			{
				try {
					out.write(arr[y][j]+"\t");
					out.flush();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			try {
				out.write("\n");
				out.flush();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
		System.out.println("");
	}

	public void printArray(double[] arr, String arrName)
	{
		System.out.print(arrName+" ");
		for(int y=0; y<arr.length;y++)
		{
			System.out.println(y+"\t"+arr[y]+"\t");
		}
		System.out.println("");
	}

	//TODO: Now the network can compute its own current q(k) and q(k/k')
	//It should have functions to calculate its H(q), H(q/q') amd I(q) for its current state.

	/**
	 * Calculate the Entropy of this Q distribution, H(Q)
	 * Answer in Bits
	 */
	public double calculateHq() {
		double h = 0;
		double[] qd = getRemainingDegreeDistribution();
		for (int i = 0; i < qd.length; i++) {
			if (qd[i] != 0)
			{
				h = h - (qd[i] * Math.log(qd[i]));
			}
		}
		return (h/Math.log(2)); /* Convert to bits befreo returning */
	}

	/**
	 * Calculate the conditional Entropy of this Q distribution, H(Q)
	 * Answer in Bits
	 */
	public double calculateHqc()
	{
		double in =  0;
		double[][] pi = getConditionalDegreeDistribution();
		double[] qd = getRemainingDegreeDistribution();
		for (int i = 1; i < pi.length; i++) {
			for (int j = 1; j < pi[i].length; j++) {
				if ( (pi[i][j] != 0))
				{
					in = in - (  qd[i-1]*pi[i][j] * Math.log(pi[i][j])  ); //TODO: Check this logic.
				}
			}
		}

		return ( in / Math.log(2)); /* Convert to bits befreo returning */
	}

	/**
	 * Calculate the mutual Information, I(q) = H(q) - H(q/q') of this network.
	 * @return
	 * TODO: Check logic, change to the unit that Sole uses (bits ?)
	 */
	public double calculateMI()
	{
		return ( calculateHq() - calculateHqc() );
	}

	/**
	 * Calculate MI WITHOUT CALCULATING entropies
	 * @return
	 */
	public double calculate_MI_straight() {

		double[][] ejk = getLinkDistribution();
		double[] qk = getRemainingDegreeDistribution();
		double i = 0;

		for (int j = 0; j < qk.length; j++)
		{
			for (int k = 0; k < qk.length; k++)
			{  double addition =  ejk[j][k] / (qk[j]*qk[k]);
			if(addition != 0)
			{
				if( (qk[j] != 0) && (qk[k] != 0) )
				{
					i = i + (ejk[j][k])*Math.log(  ejk[j][k] / (qk[j]*qk[k])  );
				}
			}
			}

		}
		i  =  i /Math.log(2);

		return i;
	}


	/**
	 * Calculate MI assuming the network is directed
	 * @return
	 */
	public double calculate_MI_DirectedNet() {

		double[][] ejk = getDirectedLinkDistribution();
		double[] qkin = getInDegreeDistribution();
		double[] qkout = getOutDegreeDistribution();
		double i = 0;

		for (int j = 0; j < qkin.length; j++)
		{
			for (int k = 0; k < qkout.length; k++)
			{  double addition =  ejk[j][k] / (qkin[j]*qkout[k]);
			if(addition != 0)
			{
				if( (qkin[j] != 0) && (qkout[k] != 0) )
				{
					i = i + (ejk[j][k])*Math.log(  ejk[j][k] / (qkin[j]*qkout[k])  );
				}
			}
			}

		}
		i  =  i /Math.log(2);

		return i;
	}

	public double calculate_MI_In_DirectedNet() {

		double[][] ejk = getDirectedLinkInInDistribution();
		double[] qkin = getInDegreeDistribution();
		double[] qkdashin = getInDegreeDistribution_dash();
		double i = 0;

		for (int j = 0; j < qkin.length; j++)
		{
			for (int k = 0; k < qkdashin.length; k++)
			{  double addition =  ejk[j][k] / (qkin[j]*qkdashin[k]);
			if(addition != 0)
			{
				if( (qkin[j] != 0) && (qkdashin[k] != 0) )
				{
					i = i + (ejk[j][k])*Math.log(  ejk[j][k] / (qkin[j]*qkdashin[k])  );
				}
			}
			}

		}
		i  =  i /Math.log(2);

		return i;
	}

	public double calculate_MI_Out_DirectedNet() {

		double[][] ejk = getDirectedLinkOutOutDistribution();
		double[] qkout = getOutDegreeDistribution();
		double[] qkdashout = getOutDegreeDistribution_dash();
		double i = 0;

		for (int j = 0; j < qkout.length; j++)
		{
			for (int k = 0; k < qkdashout.length; k++)
			{  double addition =  ejk[j][k] / (qkout[j]*qkdashout[k]);
			if(addition != 0)
			{
				if( (qkout[j] != 0) && (qkdashout[k] != 0) )
				{
					i = i + (ejk[j][k])*Math.log(  ejk[j][k] / (qkout[j]*qkdashout[k])  );
				}
			}
			}

		}
		i  =  i /Math.log(2);

		return i;
	}

	/**
	 * Calculate and return the average degree of this network based on this  degree distribution
	 */
	public double getAverageDegree()
	{
		double[] pDist = getDegreeDistribution();
		double avg = 0;
		for(int i=0; i< pDist.length;i++)
		{
			avg = avg + (i*pDist[i]);
		}
		return avg;
	}


	/**
	 * Calculate and return the average BC of this network
	 */
	public double getAverageBC()
	{

		double bc = 0;
		for(int i=0; i< nodes.size();i++)
		{
			bc = bc + nodes.elementAt(i).getBCCentrality();
		}

		bc = bc / (double)this.getSize();
		return bc;
	}


	/**
	 * Calculate and return the average BC of this network
	 */
	public double getAverageCC()
	{

		double cc = 0;
		for(int i=0; i< nodes.size();i++)
		{
			cc = cc + nodes.elementAt(i).getCCCentrality();
		}

		cc = cc / (double)this.getSize();
		return cc;
	}


	/**
	 * Calculate and return the Gamma of this network based on this degree distribution, assuming it is scale free!
	 * this method is wrongly implemented -  see next
	 */
	public double getGamma_wrong()
	{
		double[] pDist = getDegreeDistribution(); ///TODO: Complete This!
		double gamma = 1;
		double min_error = 10000000;
		double fit_gamma = 0;
		for(int i=0; i< 800;i++)
		{
			gamma = 0.5  + 0.01*i;



			//Calculate Error
			double sqError = 0;
			for(int j=1; j< pDist.length;j++)
			{
				double pr = pDist[1]*Math.pow(j, -1.00*gamma); //Predicted Gamma: problem: assumes going through first point.
				sqError = sqError + ((pr - pDist[j])* (pr - pDist[j]));
			}

			//Compare squared Errors
			if(sqError < min_error)
			{
				min_error = sqError;
				fit_gamma = gamma;
			}
		}
		return fit_gamma;
	}

	public double getGamma()
	{
		double[] pDist = getDegreeDistribution(); ///TODO: Complete This!

		//fill out x and y values

		double[] x = new double[pDist.length];
		double[] y= new double[pDist.length];
		x[0] = 0;
		y[0] = 0;
		for(int j=1; j< pDist.length;j++)
		{
			x[j] = Math.log(j);
			if(pDist[j]!=0)
			{
				y[j] = Math.log(pDist[j]);
			}
			else
			{
				y[j] = 0;
			}
		}

		//Get the required sums

		double sum_xx = 0;
		double sum_xy = 0;
		double sum_x = 0;
		double sum_y = 0;

		double n = x.length;
		double fit_gamma = 0;

		for(int i=1; i< x.length;i++)
		{

			sum_xx = sum_xx + x[i]*x[i];
			sum_xy = sum_xy + x[i]*y[i];
			sum_x = sum_x + x[i];
			sum_y = sum_y + y[i];
		}

		double deno = n*sum_xx - sum_x*sum_x;
		double nume = sum_y*sum_xx - sum_x*sum_xy;
		fit_gamma = nume / deno;
		return -fit_gamma;
	}


	/**
	 * Calculate and return the cutoff of degree distribution
	 * of this network
	 */
	public double getCutoff()
	{
		double[] pDist = getDegreeDistribution();
		double Np = 0;
		for(int i=0; i< pDist.length;i++)
		{
			if(pDist[i] > 0.0000000001 ) //small number
				Np = i;
		}
		return Np;
	}


	/**
	 * This Method is included to change a network to standard text format,
	 * which can then be read into tools like Pajek, Cytoscape and Walrus
	 * Each link is printed on a separate line, with Node IDs separated by space
	 * @since 25.01.2008
	 * @param fileName
	 */
	public void createTextFile(String fileName)
	{
		BufferedWriter out = null;

		try {

			out = new BufferedWriter(new FileWriter(fileName+".txt"));
			for(int i=0; i< links.size();i++)
			{
				Link thisLink = links.elementAt(i);


				long xID = thisLink.iNodeID();
				long yID = thisLink.jNodeID();
				out.write(xID+" "+yID+"\n");
				//System.out.println(xID+" "+yID);
				out.flush();
			}
			//System.out.println("Created Txt File " + fileName + ".txt");
		}

		catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}
	}

	public void printToConsole()
	{




		for(int i=0; i< links.size();i++)
		{
			Link thisLink = links.elementAt(i);
			long xID = thisLink.iNodeID();
			long yID = thisLink.jNodeID();

			System.out.println(xID+" "+yID);

		}



	}

	/**
	 * Find out the number of UNIQUE NEIGHBOURS to a node
	 * based on the link list maintained by network
	 * Self links and double links should be omitted
	 * @param nodeID
	 * @return
	 */
	public int getNeighbourNumber(long nodeID)
	{
		int no = 0;

		boolean noDoubleLink =  true;

		for(int i=0; i<this.links.size();i++)
		{
			noDoubleLink =  true;
			Link thisLink = this.links.elementAt(i);
			if(thisLink.iNodeID()== thisLink.jNodeID()) //self link
			{
				// do not count
				continue;
			}

			//Identify the neighbour in this Link
			Node nei1 =  thisLink.iNode();
			if(nei1.getID() == nodeID)
			{
				nei1 = thisLink.jNode(); //Now we have identified the neighbour in this link
				if(nei1.getID() == nodeID) //Just a double check
				{
					System.out.println("SELF LINK: Shouldn't have occured here");
				}
			}

			for(int j=0; j< i; j++) //check already traversed elements for double entry
			{
				Link otherLink = this.links.elementAt(j);
				Node nei2 =  otherLink.iNode();
				if(nei2.getID() == nodeID)
				{
					nei2 = otherLink.jNode(); //Now we have identified the neighbour in this link
					if(nei2.getID() == nodeID) //Just a double check
					{
						System.out.println("SELF LINK: Shouldn't have occured here");
					}
				}

				//Is this neighbour already identified?
				if(nei1.getID() == nei2.getID()) //double link
				{
					noDoubleLink =  false;
					break;
				}
			}

			//Self-links have been checked for ith element
			if(noDoubleLink)
			{
				no++;
			}


		}

		//Now all links which were not self links and double links have been copied into l

		return no;

	}


	/**
	 * Return the list of UNIQUE NEIGHBOURS to a node
	 * based on the link list maintained by network
	 * Self links and double links should be omitted
	 * @param nodeID
	 * @return
	 */
	public Vector<Link> getUniqueNeighbourList(long nodeID)
	{
		int no = 0;
		Vector<Link> l = new Vector<Link>();
		boolean noDoubleLink =  true;

		for(int i=0; i<this.links.size();i++)
		{
			noDoubleLink =  true;
			Link thisLink = this.links.elementAt(i);
			if(thisLink.iNodeID()== thisLink.jNodeID()) //self link
			{
				// do not count
				continue;
			}

			//Identify the neighbour in this Link
			Node nei1 =  thisLink.iNode();
			if(nei1.getID() == nodeID)
			{
				nei1 = thisLink.jNode(); //Now we have identified the neighbour in this link
				if(nei1.getID() == nodeID) //Just a double check
				{
					System.out.println("SELF LINK: Shouldn't have occured here");
				}
			}

			for(int j=0; j< i; j++) //check already traversed elements for double entry
			{
				Link otherLink = this.links.elementAt(j);
				Node nei2 =  otherLink.iNode();
				if(nei2.getID() == nodeID)
				{
					nei2 = otherLink.jNode(); //Now we have identified the neighbour in this link
					if(nei2.getID() == nodeID) //Just a double check
					{
						System.out.println("SELF LINK: Shouldn't have occured here");
					}
				}

				//Is this neighbour already identified?
				if(nei1.getID() == nei2.getID()) //double link
				{
					noDoubleLink =  false;
					break;
				}
			}

			//Self-links have been checked for ith element
			if(noDoubleLink)
			{
				l.addElement(this.links.elementAt(i));
			}


		}

		//Now all links which were not self links and double links have been copied into l

		return l;

	}

	public void mergeNetworksIntoTextFile(String fileName, Network other, long sourceID, long targetID)
	{
		BufferedWriter out = null;

		try {

			out = new BufferedWriter(new FileWriter(fileName+".txt"));
			for(int i=0; i< links.size();i++)
			{
				Link thisLink = links.elementAt(i);
				long xID = thisLink.iNodeID();
				long yID = thisLink.jNodeID();
				out.write(xID+" "+yID+"\n");
				out.flush();
			}

			for(int i=0; i< other.links.size();i++)
			{
				Link thisLink = other.links.elementAt(i);
				long xID = thisLink.iNodeID() + getSize()+1; //id conversion
				long yID = thisLink.jNodeID() + getSize()+1;
				out.write(xID+" "+yID+"\n");
				out.flush();
			}
			sourceID = sourceID  + getSize()+1; //convert target ID to the new ID system
			// Now print the merging link
			out.write(sourceID+" "+targetID+"\n");
			out.flush();


			//System.out.println("Created Txt File " + fileName + ".txt");
		}

		catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}
	}


	public double  getWeightedDegreeSum()
	{

		/* Length of Q(k), and therefore E(j,k), is one less than length of degree
		distribution, which is (maxDegree+1), allowing for zero degree */
		double sum = 0;
		for(int i=0;i<nodes.size();i++) 
		{			
			Node thisNode = nodes.elementAt(i);
			sum = sum + (thisNode.getNumberOfInLinks()*thisNode.getWeight());
		}

		return sum;
	}

	/*
	 * This method calculates and returns all the assortative Links in the network
	 * Assortative Links are those which have the sum of LA of their nodes as a positive value
	 */
	public Vector<Link> assortativeLinks(){
		Vector<Link>  alinks =  new Vector();
		for(int i=0;i<links.size();i++) 
		{
			Link thisLink = links.elementAt(i);
			if(localR(thisLink) > 0 )

			{
				alinks.addElement(thisLink);
			}
		}

		return alinks;

	}

	public double localR(Link link)
	{
		return Math.abs((calculate_localR(link.iNode()) - calculate_localR(link.jNode())));
	}

	public double calculate_localR(Node node)
	{
		if(node.getNumberOfLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		double N = getSize();
		double Nlinks = getNoOfLinks();
		QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(this);
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
		double mew_ql = ((j+1)/(2*Nlinks))*mew_q*mew_q;
		double rl = (ejkl -mew_ql )/ var_q;
		return rl;
	}

	//THE FOLLOWING FUNCTIONS ARE USED IN  ANALYZING RUCH CLUB PHENOMENA

	/**
	 * THis function determines if a node can be considered ' Rich', given a rich club percentage
	 * THis function assumes the ndoes are already sorted according to degree
	 */
	public boolean isRich(Node node, double percent)
	{
		int richNodes = (int)(getSize()*percent*0.01);

		int cutOffDegree =  nodes.elementAt(richNodes).getNumberOfLinks();
		if(node.getNumberOfLinks() >= cutOffDegree)
		{
			return true;
		}
		return false;
	}

	public boolean isRichLink(Link link, double percent)
	{
		if(isRich(link.iNode(), percent) && isRich(link.jNode(), percent) )
		{
			return true;
		}

		return false;
	}

	/**
	 * THis method returns the rich club coefficient, in percenttage,
	 * for a given number of rich club percentage
	 * @param percent
	 * @return
	 */
	public double richClubCoefficient(double percent)
	{
		sortNodes();
		int richNodes = (int)(getSize()*percent*0.01);
		double possibleRichLinks = richNodes*(richNodes-1)*0.5;

		double richLinks = 0;
		for(int i=0; i< links.size(); i++)
		{
			if(isRichLink(links.elementAt(i), percent))
			{
				richLinks++;
			}
		}

		double phi = 100*richLinks / possibleRichLinks;

		return phi;
	}

	/**
	 * This function returns a double array of rich club profile of a network
	 * The profile is as percentage points, from 1% to 100%.
	 */
	public double[][] richClubProfile()
	{
		double arr[][] = new double[100][2];
		sortNodes();
		for(int i=1; i< 100; i++)
		{
			arr[i-1][0] = i;
			arr[i-1][1] = richClubCoefficient(i);
		}

		return arr;
	}

	//FROM HERE ONWARDS, THE FUNCTIONS IMPLEMENTED ARE MODIFIED VERSIONS OF NAMAL SENARATNE'S CODE
	//These functions implement shortest path calculations and betweenness centrality calculations

	/**
	 * Initialize the graph for Shortest Path calculation
	 *
	 * @param   sourceNodeID    the source node used to calculate the shortest paths
	 */
	public void initializeForSP(int sourceNodeID)
	{
		int MAX_SIZE =  getSize();
		int i = 0;   
		while(i < MAX_SIZE && nodes.elementAt(i) != null)
		{
			/* set the distance from source to inifinity */
			nodes.elementAt(i).setDistance(INF);     

			/* clear all the predecessor information, so we can start out fresh */
			nodes.elementAt(i).clearPredecessors();  

			/* none of the nodes are in the queue, so the state is set to NOT_QUEUED */
			nodes.elementAt(i).setNodeState(NodeState.NOT_ENQUEUED);

			i++;
		}

		/* distance of source node is set to zero */
		try {
			getNode(sourceNodeID).setDistance(0);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			System.out.println("Node with given ID "+ sourceNodeID + " DOES NTO EXSIST");
			e.printStackTrace();
		} 

		//System.out.println("finished initialization of ndoes");
	}

	/**
	 * Initialize the graph for Shortest Path calculation
	 *
	 * @param   sourceNodeID    the source node used to calculate the shortest paths
	 */
	public void initializeForIndirectAdjacency(int sourceNodeID)
	{
		int MAX_SIZE =  getSize();
		int i = 0;   
		while(i < MAX_SIZE && nodes.elementAt(i) != null)
		{
			/* set the distance from source to inifinity */
			nodes.elementAt(i).setDistance(INF);     

			/* clear all the predecessor information, so we can start out fresh */
			// nodes.elementAt(i).clearPredecessors();  

			/* none of the nodes are in the queue, so the state is set to NOT_QUEUED */
			// nodes.elementAt(i).setNodeState(NodeState.NOT_ENQUEUED);

			i++;
		}

		/* distance of source node is set to zero */
		try {
			getNode(sourceNodeID).setDistance(0);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			System.out.println("Node with given ID "+ sourceNodeID + " DOES NTO EXSIST");
			e.printStackTrace();
		} 

		//System.out.println("finished initialization of ndoes");
	}


	/**
	 * Calculates the shorted paths from the specified source node to all <br \>
	 * other nodes. The algorithm used is Breadth First Search since the graph is
	 * a unweighted di-graph
	 *
	 * @param   sourceNodeID    the source node
	 * @throws Exception 
	 */
	public void calculateShortestPath(int sourceNodeID) throws Exception
	{
		/* Initialize the graph for shortest path calculation with souceNode */
		initializeForSP(sourceNodeID);

		/* The queue used as the FIFO queue */
		LinkedBlockingQueue<Node> queue = new LinkedBlockingQueue<Node>();      

		/* start by putting the source node to the queue */
		Node sourceNode = getNode(sourceNodeID);
		queue.put(sourceNode);
		sourceNode.setNodeState(NodeState.QUEUED);

		while(queue.size() > 0)     /* iterate until the queue is not empty */
		{
			Node u = queue.poll();  /* retrieve and remove the head of the queue */

			/* retrieve the adjacent nodes of the dequeued node */
			List<Integer> adjacency = u.getAdjacency();
			Iterator<Integer> itadjacency = adjacency.iterator();

			/* for each adjacent node, calculate the distance */
			while(itadjacency.hasNext())
			{
				Node v = getNode(itadjacency.next().intValue());

				/* Check whether the node's distance value is freezed 
                   If the node is added to the queue and dequeued its state
                   is set to DEQUEUED, a node with its state set to DEQUEUED
                   has the shortest distance from the source node
                   So no need of adding it again to the queue */
				if(v.getNodeState() != NodeState.DEQUEUED)
				{
					/* This condition happens when the node is previously added
                       to the queue but not yet dequeued.
                       The distance of the node cannot change, as the BFS
                       expands nodes with shortest distance first. So the calculated
                       distance for the node through the predecessor 'u' could be 
                       either equal or greater than the node's current distance.
                       If it's greater, then just ignore it. If it is equal, the node
                       has another new predecessor */
					if(v.getNodeState() == NodeState.QUEUED && 
							(v.getDistance() == u.getDistance() + 1))
					{
						/* add the newly identified predecessor */
						v.addPredecessorNode((int)u.getID());    
					}

					/* This condition becomes true only when the node is not added
                       to the queue earlier */
					if(v.getNodeState() == NodeState.NOT_ENQUEUED)
					{
						v.setDistance(u.getDistance() + 1);
						v.addPredecessorNode((int)u.getID());
						queue.put(v);
						v.setNodeState(NodeState.QUEUED);
					}
				}
			}
			u.setNodeState(NodeState.DEQUEUED);
		}     

		Node node = getNode(sourceNodeID);
		//System.out.println("shortest paths calculated");
	}   

	/**
	 * Prints the shortest paths from sourceNodeID to rest of the nodes
	 */
	public void printShortestPath(int sourceNodeID)
	{
		int i = 0;
		while(nodes.elementAt(i) != null)
		{
			System.out.println("shortest path lenght for node " + i + " : " + nodes.elementAt(i).getDistance());
			System.out.println("\tShortest paths...");


			List<List> paths1 = getShortestPaths(i, sourceNodeID);
			for(int j = 0; j < paths1.size(); j++)
			{
				Iterator itVertexList = ((List)paths1.get(j)).iterator();
				System.out.print("\t\t");
				while(itVertexList.hasNext())
				{
					System.out.print(((Integer)itVertexList.next()).intValue() + " ");
				}
				System.out.println();
			}

			i++;
		}
	}


	/**
	 * Return a list of shortest paths from source node to destination node<br\>
	 * Depth First Search is used to traverse from destination node to <br \>
	 * source node to calculate all of the shortest paths.
	 *
	 * @param   destNodeID  the destination node's ID
	 * @param   sourceNode  the source node's ID
	 */
	public List<List> getShortestPaths(int destNodeID, int sourceNode)
	{      
		/* the number of predecessors of destination node */
		try {
			int predecessorNo = getNode(destNodeID).getPredecessor().size();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			System.out.println("dest node for the shortest path doesnt exist");
			e.printStackTrace();
		}
		List<List> paths = new ArrayList<List>();

		// retrieving predecessors of the destination node
		Iterator<Integer> itnodes;
		try {
			itnodes = getNode(destNodeID).getPredecessor().iterator();


			while(itnodes.hasNext())
			{
				int pred = itnodes.next().intValue();

				// create a new path for the new predecessor
				List<Integer> path = new ArrayList<Integer>();

				// add the node numbers to identified path
				path.add(Integer.valueOf(destNodeID));
				path.add(Integer.valueOf(pred));

				// Do a Depth First Search to identify all the shortest paths from
				// source node to destination node
				DFS(path, paths, pred, sourceNode);
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return paths;
	}

	/**
	 * Return a  shortest path from source node to destination node<br\>
	 * Depth First Search is used to traverse from destination node to <br \>
	 * source node to calculate a shortest path.
	 * 
	 * This METHOD IS USED TO ONLY SEE IF A PATH EXISTS BETWEEN SOURCE AND DESTINATION
	 *
	 * @param   destNodeID  the destination node's ID
	 * @param   sourceNode  the source node's ID
	 */
	public List<List> getAShortestPath(int destNodeID, int sourceNode)
	{      
		/* the number of predecessors of destination node */
		try {
			int predecessorNo = getNode(destNodeID).getPredecessor().size();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			System.out.println("dest node for the shortest path doesnt exist");
			e.printStackTrace();
		}
		List<List> paths = new ArrayList<List>();

		// retrieving predecessors of the destination node
		Iterator<Integer> itnodes;
		try {
			itnodes = getNode(destNodeID).getPredecessor().iterator();


			while(itnodes.hasNext())
			{
				int pred = itnodes.next().intValue();

				// create a new path for the new predecessor
				List<Integer> path = new ArrayList<Integer>();

				// add the node numbers to identified path
				path.add(Integer.valueOf(destNodeID));
				path.add(Integer.valueOf(pred));

				// Do a Depth First Search to identify all the shortest paths from
				// source node to destination node
				boolean isPath =  false;
				isPath = DFS1(path, paths, pred, sourceNode);
				if(isPath)
				{
					break;
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return paths;
	}


	public void updateAllPairSPArrays_rc(int sourceNodeID, float[] BCDistribution, double weighttotal) throws Exception
	{
		int gSize = getSize();
		int sigmaSTV, sigmaST;

		// represents the number of shortest paths from sourceNode to all the other nodes
		// in the graph
		int[] noOfShortPaths = new int[gSize];

		// the entry (i,j) in the array represents the number shortest paths
		// between sourceNode and i, which is going through j
		int [][] noOfCentralityPaths = new int[gSize][gSize];

		int i = 0;
		while(i < getSize() && nodes.elementAt(i) != null)
		{
			int ithID = (int)nodes.elementAt(i).getID();  
			List<List> paths = getShortestPaths(ithID, sourceNodeID);

			// save the number of shortest paths from source Node to (i)th node
			noOfShortPaths[i] = paths.size();

			List path;
			int node;
			for(int j = 0; j < paths.size(); j++)
			{
				path = paths.get(j);
				for(int k = 0; k < path.size(); k++)
				{
					node = ((Integer)path.get(k)).intValue();
					if(node != sourceNodeID && node != ithID)            /* ignore the source and destination nodes */
					{

						//spCentralityArray[sourceNodeID][i][node]++; 
						noOfCentralityPaths[ithID][node] ++;
					}
				}
			}
			i++;
		}

		// update the BCDistribution array
		for(int v = 0; v < gSize; v++)  // for each v
		{
			double w_s = getNode(sourceNodeID).getRiskLevel();

			double w_v = getNode(v).getRiskLevel();

			double scale = w_s / (weighttotal -  w_v);
			for(int t = 0; t < gSize; t++)
			{
				if((sourceNodeID != t) && (sourceNodeID != v) && (t != v))
				{
					sigmaSTV = noOfCentralityPaths[t][v];
					sigmaST = noOfShortPaths[t];
					float delta = (float)sigmaSTV / (float)sigmaST;
					if(sigmaST != 0)
					{
						BCDistribution[v] += (float)(scale*delta);
					}
				}
			}
		}
	}



	public void updateAllPairSPArrays_rc_st(int sourceNodeID, float[] BCDistribution, double weight_diff_total) throws Exception
	{
		int gSize = getSize();
		int sigmaSTV, sigmaST;

		// represents the number of shortest paths from sourceNode to all the other nodes
		// in the graph
		int[] noOfShortPaths = new int[gSize];

		// the entry (i,j) in the array represents the number shortest paths
		// between sourceNode and i, which is going through j
		int [][] noOfCentralityPaths = new int[gSize][gSize];

		int i = 0;
		while(i < getSize() && nodes.elementAt(i) != null)
		{
			int ithID = (int)nodes.elementAt(i).getID();  
			List<List> paths = getShortestPaths(ithID, sourceNodeID);

			// save the number of shortest paths from source Node to (i)th node
			noOfShortPaths[i] = paths.size();

			List path;
			int node;
			for(int j = 0; j < paths.size(); j++)
			{
				path = paths.get(j);
				for(int k = 0; k < path.size(); k++)
				{
					node = ((Integer)path.get(k)).intValue();
					if(node != sourceNodeID && node != ithID)            /* ignore the source and destination nodes */
					{

						//spCentralityArray[sourceNodeID][i][node]++; 
						noOfCentralityPaths[ithID][node] ++;
					}
				}
			}
			i++;
		}

		// update the BCDistribution array
		for(int v = 0; v < gSize; v++)  // for each v
		{

			//Now we need to calculate the scale factor for s-t risk centrality 

			//Remove nodes with ID  = v from the global sum

			double total_st_removal  = 0;
			double st_r =  0;
			for(int r=0; r< getSize(); r++)
			{
				st_r  = getNode(r).getRiskLevel() - getNode(v).getRiskLevel();
				if(st_r < 0)
				{
					st_r = 0;
				}

				total_st_removal   = total_st_removal   + st_r;

				st_r  = getNode(v).getRiskLevel() - getNode(r).getRiskLevel();
				if(st_r < 0)
				{
					st_r = 0;
				}

				total_st_removal   = total_st_removal   + st_r;

			}

			double subtotal_st_diff = weight_diff_total - total_st_removal;
			double scale = 1.0;



			for(int t = 0; t < gSize; t++)
			{

				//Calculate the scaling factor;
				double local_st_diff   = getNode(sourceNodeID).getRiskLevel() - getNode(t).getRiskLevel();
				if(local_st_diff < 0)
				{
					local_st_diff = 0;
				}

				if(subtotal_st_diff!=0)
				{
					scale = local_st_diff / subtotal_st_diff;
				}
				else
				{
					scale = 0;
				}

				if((sourceNodeID != t) && (sourceNodeID != v) && (t != v))
				{
					sigmaSTV = noOfCentralityPaths[t][v];
					sigmaST = noOfShortPaths[t];
					float delta = (float)sigmaSTV / (float)sigmaST;
					if(sigmaST != 0)
					{
						BCDistribution[v] += (float)(scale*delta);
					}
				}
			}
		}
	}



	public void updateAllPairSPArrays(int sourceNodeID, float[] BCDistribution) throws Exception
	{
		int gSize = getSize();
		int sigmaSTV, sigmaST;

		// represents the number of shortest paths from sourceNode to all the other nodes
		// in the graph
		int[] noOfShortPaths = new int[gSize];

		// the entry (i,j) in the array represents the number shortest paths
		// between sourceNode and i, which is going through j
		int [][] noOfCentralityPaths = new int[gSize][gSize];

		int i = 0;
		while(i < getSize() && nodes.elementAt(i) != null)
		{
			int ithID = (int)nodes.elementAt(i).getID();
			List<List> paths = getShortestPaths(ithID, sourceNodeID);

			// save the number of shortest paths from source Node to (i)th node
			noOfShortPaths[i] = paths.size();

			List path;
			int node;
			for(int j = 0; j < paths.size(); j++)
			{
				path = paths.get(j);
				for(int k = 0; k < path.size(); k++)
				{
					node = ((Integer)path.get(k)).intValue();
					if(node != sourceNodeID && node != ithID)            /* ignore the source and destination nodes */
					{

						//spCentralityArray[sourceNodeID][i][node]++; 
						noOfCentralityPaths[ithID][node] ++;
					}
				}
			}
			i++;
		}

		// update the BCDistribution array
		for(int v = 0; v < gSize; v++)  // for each v
		{
			for(int t = 0; t < gSize; t++)
			{
				if((sourceNodeID != t) && (sourceNodeID != v) && (t != v))
				{
					sigmaSTV = noOfCentralityPaths[t][v];
					sigmaST = noOfShortPaths[t];
					if(sigmaST != 0)
					{
						BCDistribution[v] += (float)(sigmaSTV)/sigmaST;
					}
				}
			}
		}
	}


	public int[][] returnPathArrays(int sourceNodeID, int[][] noOfPathsArray  ) throws Exception
	{

		int gSize = getSize();


		int sigmaSTV, sigmaST;


		int i = 0;
		while(i < getSize() && nodes.elementAt(i) != null)
		{
			int ithID = (int)nodes.elementAt(i).getID();
			int check = (int)nodes.elementAt(i).getNumberOfLinks();
			List<List> paths = getAShortestPath((int)nodes.elementAt(i).getID(), sourceNodeID);

			sigmaST = paths.size();
			if(sigmaST != 0)
			{
				noOfPathsArray[sourceNodeID][ithID] = sigmaST;
			}

			i++;
			continue;


		}

		return noOfPathsArray;


	}

	/**
	 * Perform Depth First Search from the nodeID to sourceNode
	 *
	 * @param   vertexList  the list which stores a collection of vertices
	 *                      in a single shortest path
	 * @param   pathList    the list which stores all the shortest paths from
	 *                      sourceNode to destination (destination node is specified)
	 *                      in the first call to DFS()
	 * @param   nodeID      the node ID to start the current iteration of DFS
	 *                      The initial call to DFS will have the destination node's ID 
	 *                      as the value for this argument
	 * @param   sourceNode  the source node from which the shortest path begins
	 * @throws Exception 
	 */
	private void DFS(List vertexList, List pathList, int nodeID, int sourceNode) throws Exception
	{
		int pred;

		if(nodeID == sourceNode)        /* termination of one traversal */
		{
			pathList.add(vertexList);   /* save the identified shortest path */
			return;
		}
		else
		{            
			// retrieving the list of predecessors for the intermediate node
			// in the shortest path
			Iterator<Integer> itnodes = getNode(nodeID).getPredecessor().iterator();


			List<Integer> newList;
			while(itnodes.hasNext())  /* for each predecessor a new shortest path exists */
			{
				pred = itnodes.next().intValue();

				// creates a new copy of the vertext list found so far
				// this allows the next invocation of the same loop 
				// to use the vertices found before comming to this intermediate node
				newList = copyList(vertexList);

				// add one predecessor to the newly created vertex list
				// e.g. If vertexList := 10 7 5, then after this line
				//      newList becomes := 10 7 5 3, given 3 is the newly identified predecessor
				newList.add(Integer.valueOf(pred));

				// Do DFS through new predecessor
				// So the new DFS will be performed on the path 10 7 5 3, considering the above example
				// If 3 and 9 are the predecessor values of 5, 
				// the next invocation of the following DFS will traverse through 10 7 5 9
				// likewise, this line will execute for each new predecessor value
				DFS(newList, pathList, pred, sourceNode);
			}
		}
	}

	/**
	 * Perform Depth First Search from the nodeID to sourceNode
	 *
	 * @param   vertexList  the list which stores a collection of vertices
	 *                      in a single shortest path
	 * @param   pathList    the list which stores all the shortest paths from
	 *                      sourceNode to destination (destination node is specified)
	 *                      in the first call to DFS()
	 * @param   nodeID      the node ID to start the current iteration of DFS
	 *                      The initial call to DFS will have the destination node's ID 
	 *                      as the value for this argument
	 * @param   sourceNode  the source node from which the shortest path begins
	 * thus method returns a boolean upon finding a shortest path, rather than void.
	 * @throws Exception 
	 */
	private boolean DFS1(List vertexList, List pathList, int nodeID, int sourceNode) throws Exception
	{
		int pred;

		if(nodeID == sourceNode)        /* termination of one traversal */
		{
			pathList.add(vertexList);   /* save the identified shortest path */
			return true;
		}
		else
		{            
			// retrieving the list of predecessors for the intermediate node
			// in the shortest path
			Iterator<Integer> itnodes = getNode(nodeID).getPredecessor().iterator();


			List<Integer> newList;
			while(itnodes.hasNext())  /* for each predecessor a new shortest path exists */
			{
				pred = itnodes.next().intValue();

				// creates a new copy of the vertext list found so far
				// this allows the next invocation of the same loop 
				// to use the vertices found before comming to this intermediate node
				newList = copyList(vertexList);

				// add one predecessor to the newly created vertex list
				// e.g. If vertexList := 10 7 5, then after this line
				//      newList becomes := 10 7 5 3, given 3 is the newly identified predecessor
				newList.add(Integer.valueOf(pred));

				// Do DFS through new predecessor
				// So the new DFS will be performed on the path 10 7 5 3, considering the above example
				// If 3 and 9 are the predecessor values of 5, 
				// the next invocation of the following DFS will traverse through 10 7 5 9
				// likewise, this line will execute for each new predecessor value
				DFS1(newList, pathList, pred, sourceNode);
			}
		}

		return false;
	}

	/**
	 * Creates an exact copy of the array list used to hold vertex list for 
	 * a pariticular shortest path
	 *
	 * @param   vertexList  the ArrayList which needs to be copied
	 *
	 * @return  the cloned arrayList (exact copy)
	 */
	private List<Integer> copyList(List<Integer> vertexList)
	{
		List<Integer> newList = new ArrayList<Integer>();

		Integer val;
		for(int i = 0; i < vertexList.size(); i++)
		{
			val = new Integer(((Integer)vertexList.get(i)).intValue());
			newList.add(val);
		}
		return newList;
	}

	/**
	 * Calculates the betweenness centrality indices for all of the vertices
	 * in the graph
	 *
	 * <p>
	 * This is the implementation of the algorithm proposed in the research paper
	 * "A Faster Algorithm for Betweenness Centrality" by <br\>
	 * Ulrik Brandes, Department of Computer & Information Science, University of Konstanz
	 * This algorithm calculates the Betweenness Centrality indices for an acyclic
	 * graph in O(mn) time where n is the number of verices and m is the number of
	 * edges 
	 * The memory complexity of the algorithm is O(n^2 + nm) for acyclic graphs
	 * @throws Exception 
	 */
	public float[] calcBetweenessCentrality() throws Exception
	{
		int gSize = getSize();

		float[] BCentrality = new float[gSize];
		for(int i = 0; i < gSize; i++)   // for each vertex in the graph
		{
			int ithID = (int)nodes.elementAt(i).getID();
			initializeForSP(i);

			Stack<Node> S = new Stack<Node>();

			// predInSP[w] gives the list of predecessors of w who lies between 
			// the shortest path between sourcenode (denoted by i) and w
			List[] predInSP = new ArrayList[gSize];
			for(int j = 0; j < predInSP.length; j++)
			{
				int jthID = (int)nodes.elementAt(j).getID();
				predInSP[j] = new ArrayList<Node>();
			}

			// the entry sigma[t] represents the number shortest pathes from
			// source (denoted by i) to t
			int[] sigma = new int[gSize];
			sigma[i] = 1;

			LinkedBlockingQueue<Node> Q = new LinkedBlockingQueue<Node>();

			// put the source node to queue to initialize the queue
			Q.put(getNode(i));

			while(Q.size() > 0)
			{
				Node node = (Node)Q.poll();
				S.push(node);

				ArrayList<Integer> adjacency = node.getAdjacency();
				for(Integer nodeID : adjacency)
				{
					int neighbourID = nodeID.intValue();
					Node neighbour = getNode(neighbourID);

					// neighbour found for the first time
					if(neighbour.getDistance() == Integer.MAX_VALUE) 
					{
						Q.put(neighbour);
						neighbour.setDistance(node.getDistance() + 1);
					}

					if(neighbour.getDistance() == (node.getDistance() + 1))
					{
						sigma[neighbourID] += sigma[(int)node.getID()];
						predInSP[neighbourID].add(node);
					}
				}
			}

			//here we can fill the 'indirect-adjacency' matrix

			float[] delta = new float[gSize];
			while(!S.empty())
			{
				Node node = S.pop();
				int nodeID = (int)node.getID();
				List<Node> preds = (List)predInSP[nodeID];
				for(Node pred : preds)    // for each predecessor of 'node'
				{
					int predID = (int)pred.getID();
					delta[predID] += (float)(sigma[predID]*(1 + delta[nodeID]))/sigma[nodeID];
				}

				if(nodeID != i)
				{
					double Nscale = (getSize() - 2)*(getSize() - 1);

					BCentrality[nodeID] += (delta[nodeID]) / (Nscale); // scaling 
				}
			}
		}

		return BCentrality;
	}

	/**
	 * Calculates the CLOSENESS centrality indices for all of the vertices
	 * in the graph
	 *
	 * <p>
	 * This is the MODIFIED implementation of the algorithm proposed in the research paper
	 * "A Faster Algorithm for Betweenness Centrality" by <br\>
	 * Ulrik Brandes, Department of Computer & Information Science, University of Konstanz
	 * This algorithm calculates the Closeness Centrality indices for an acyclic
	 * graph in O(mn) time where n is the number of verices and m is the number of
	 * edges 
	 * The memory complexity of the algorithm is O(n^2 + nm) for acyclic graphs
	 * @throws Exception 
	 */
	public float[] calcClosenessCentrality() throws Exception
	{
		int gSize = getSize();

		float[] CCentrality = new float[gSize]; //Varaiable name not changed
		for(int i = 0; i < gSize; i++)   // for each vertex in the graph
		{
			int ithID = (int)nodes.elementAt(i).getID();
			initializeForSP(i);

			Stack<Node> S = new Stack<Node>();

			// predInSP[w] gives the list of predecessors of w who lies between 
			// the shortest path between sourcenode (denoted by i) and w
			List[] predInSP = new ArrayList[gSize];
			for(int j = 0; j < predInSP.length; j++)
			{
				int jthID = (int)nodes.elementAt(j).getID();
				predInSP[j] = new ArrayList<Node>();
			}

			// the entry sigma[t] represents the LENGTH of shortest pathes from
			// source (denoted by i) to t
			int[] sigma = new int[gSize];
			sigma[i] = 0;

			LinkedBlockingQueue<Node> Q = new LinkedBlockingQueue<Node>();

			// put the source node to queue to initialize the queue
			Q.put(getNode(i));

			while(Q.size() > 0)
			{
				Node node = (Node)Q.poll();
				S.push(node);

				ArrayList<Integer> adjacency = node.getAdjacency();
				for(Integer nodeID : adjacency)
				{
					int neighbourID = nodeID.intValue();
					Node neighbour = getNode(neighbourID);

					// neighbour found for the first time
					if(neighbour.getDistance() == Integer.MAX_VALUE) 
					{
						Q.put(neighbour);
						neighbour.setDistance(node.getDistance() + 1);
					}

					if(neighbour.getDistance() == (node.getDistance() + 1))
					{
						sigma[neighbourID] = neighbour.getDistance();
						predInSP[neighbourID].add(node);
					}
				}
			}

			//Now sigma is populated for all neighbours of node i
			//Summation should give invertion of CC

			float ccsum = 0;
			for(int j = 0; j < gSize; j++)   // for each vertex in the graph
			{
				ccsum = ccsum + sigma[j];
			}

			CCentrality[i] = (float)1.0 / ccsum;

		}

		return CCentrality;
	}


	/**
	 * Calculates the Risk centrality indices for all of the vertices
	 * in the graph
	 *
	 * <p>
	 * This implementation is similar to the implementation of the betweenness Centrality algorithm proposed in the research paper
	 * "A Faster Algorithm for Betweenness Centrality" by <br\>
	 * Ulrik Brandes, Department of Computer & Information Science, University of Konstanz
	 * This algorithm calculates the Risk Centrality indices for an acyclic
	 * graph in O(mn) time where n is the number of verices and m is the number of
	 * edges 
	 * The memory complexity of the algorithm is O(n^2 + nm) for acyclic graphs
	 * @throws Exception 
	 * 
	 * NB: This method currently returns BC values only
	 * TODO: Implement RC values
	 */
	public float[] calcRiskCentrality() throws Exception
	{
		//first calculate total risk-weight, this will be needed for scaling.
		double total_riskweight = 0;
		for(int p=0; p< getSize(); p++)
		{
			total_riskweight = total_riskweight + getNode(p).getRiskLevel();
		}

		//centrality bit
		int gSize = getSize();

		float[] BCentrality = new float[gSize];
		for(int i = 0; i < gSize; i++)   // for each vertex in the graph
		{
			int ithID = (int)nodes.elementAt(i).getID();
			initializeForSP(i);

			Stack<Node> S = new Stack<Node>();

			// predInSP[w] gives the list of predecessors of w who lies between 
			// the shortest path between sourcenode (denoted by i) and w
			List[] predInSP = new ArrayList[gSize];
			for(int j = 0; j < predInSP.length; j++)
			{
				int jthID = (int)nodes.elementAt(j).getID();
				predInSP[j] = new ArrayList<Node>();
			}

			// the entry sigma[t] represents the number shortest pathes from
			// source (denoted by i) to t
			int[] sigma = new int[gSize];
			sigma[i] = 1;

			LinkedBlockingQueue<Node> Q = new LinkedBlockingQueue<Node>();

			// put the source node to queue to initialize the queue
			Q.put(getNode(i));

			while(Q.size() > 0)
			{
				Node node = (Node)Q.poll();
				S.push(node);

				ArrayList<Integer> adjacency = node.getAdjacency();
				for(Integer nodeID : adjacency)
				{
					int neighbourID = nodeID.intValue();
					Node neighbour = getNode(neighbourID);

					// neighbour found for the first time
					if(neighbour.getDistance() == Integer.MAX_VALUE) 
					{
						Q.put(neighbour);
						neighbour.setDistance(node.getDistance() + 1);
					}

					if(neighbour.getDistance() == (node.getDistance() + 1))
					{
						sigma[neighbourID] += sigma[(int)node.getID()];
						predInSP[neighbourID].add(node);
					}
				}
			}

			//here we can fill the 'indirect-adjacency' matrix

			float[] delta = new float[gSize];
			while(!S.empty())
			{
				Node node = S.pop();
				int nodeID = (int)node.getID();
				List<Node> preds = (List)predInSP[nodeID];
				for(Node pred : preds)    // for each predecessor of 'node'
				{
					int predID = (int)pred.getID();
					delta[predID] += (float)(sigma[predID]*(1 + delta[nodeID]))/sigma[nodeID];
				}

				if(nodeID != i)
				{
					//i is the source node, get its infection weight
					double w_i = getNode(i).getRiskLevel(); 
					double w_v = getNode(nodeID).getRiskLevel(); 



					double scale = (0.5*(getSize() - 1)*(getSize() - 2)*w_v) / total_riskweight;
					scale = (w_i / (total_riskweight-w_v));
					if(total_riskweight==w_v)  //the node considered for centralioty is the only node with any risk
					{
						scale = 0;
					}


					// BCentrality[nodeID] += ((scale*delta[nodeID])/ ((getSize() - 1)*(getSize() - 2)) ); //scaled version implemented

					//new def;

					double Nscale = getSize() - 2;
					//Nscale = 37; //TODO: Delete this: only for Alberta network this is valid.

					BCentrality[nodeID] += (scale*delta[nodeID])/ (Nscale); //scaled version implemented
					int forbreakpoint = 0; //just to add breakpoint
				}
			}

			// Here we need to scale for each node based on the node's risk factor
		}

		return BCentrality;
	}

	/**
	 * Calculates the Risk centrality indices for all of the vertices
	 * in the graph
	 *
	 * <p>
	 * This implementation is similar to the implementation of the betweenness Centrality algorithm proposed in the research paper
	 * "A Faster Algorithm for Betweenness Centrality" by <br\>
	 * Ulrik Brandes, Department of Computer & Information Science, University of Konstanz
	 * The memory complexity of the algorithm is O(n^2 + nm) for acyclic graphs
	 * @throws Exception 
	 * 
	 * NB: This method calculates the weights based on both source and target node weights: hence the _st part in method name
	 * The memory complexity of the algorithm is O(n^2 + n^2m) for acyclic graphs
	 * TODO: NOTE WELL!!! "ThIS IMPLEMENTATION IS INCOMPLETE AND WRONG DUE TO QUEING SYSTEM. NEEDS TO BE COMPLETED" <= WRITTEN ON ON 1ST JULY 2011.
	 */
	public float[] calcRiskCentrality_st() throws Exception
	{
		//first calculate total of Tramp(sw - tw) this will be needed for scaling.
		double total_st_diff = 0;
		double st = 0;
		for(int s=0; s< getSize(); s++)
		{
			for(int t=0; t< getSize(); t++)
			{
				st  = getNode(s).getRiskLevel() - getNode(t).getRiskLevel();
				if(st < 0)
				{
					st = 0;
				}

				total_st_diff  = total_st_diff  + st;
			}

		}

		//centrality bit
		int gSize = getSize();

		float[] BCentrality = new float[gSize];
		for(int i = 0; i < gSize; i++)   // for each vertex in the graph
		{


			int ithID = (int)nodes.elementAt(i).getID();
			initializeForSP(i);

			Stack<Node> S = new Stack<Node>();

			// predInSP[w] gives the list of predecessors of w who lies between 
			// the shortest path between sourcenode (denoted by i) and w
			List[] predInSP = new ArrayList[gSize];
			for(int j = 0; j < predInSP.length; j++)
			{
				int jthID = (int)nodes.elementAt(j).getID();
				predInSP[j] = new ArrayList<Node>();
			}

			// the entry sigma[t] represents the number shortest pathes from
			// source (denoted by i) to t
			double[] sigma = new double[gSize];
			sigma[i] = 1;

			LinkedBlockingQueue<Node> Q = new LinkedBlockingQueue<Node>();

			// put the source node to queue to initialize the queue
			Q.put(getNode(i));

			while(Q.size() > 0)
			{
				Node node = (Node)Q.poll();
				S.push(node);

				ArrayList<Integer> adjacency = node.getAdjacency();
				for(Integer nodeID : adjacency)
				{
					int neighbourID = nodeID.intValue();
					Node neighbour = getNode(neighbourID);

					// neighbour found for the first time
					if(neighbour.getDistance() == Integer.MAX_VALUE) 
					{
						Q.put(neighbour);
						neighbour.setDistance(node.getDistance() + 1);
					}

					if(neighbour.getDistance() == (node.getDistance() + 1))
					{
						sigma[neighbourID] += sigma[(int)node.getID()];
						predInSP[neighbourID].add(node);
					}
				}
			}

			//here we can fill the 'indirect-adjacency' matrix

			double [] delta = new double[gSize];
			double [] scaledDelta = new double[gSize]; //for souce -  target weight based scale
			while(!S.empty())
			{
				Node node = S.pop();
				int nodeID = (int)node.getID();

				//node with 'NodeID' is the target node, node with 'i'th ID is the source node. Calculate the total of 'source-target' difference

				double local_st_diff = 0;
				double w_s = getNode(i).getRiskLevel(); 
				double w_t = getNode(nodeID).getRiskLevel(); 

				local_st_diff = w_s - w_t;
				if(local_st_diff< 0)
				{
					local_st_diff = 0;
				}


				List<Node> preds = (List)predInSP[nodeID];
				for(Node pred : preds)    // for each predecessor of 'node'
				{
					int predID = (int)pred.getID();

					//debugging code
					if(predID==5)
					{

						if(w_s > w_t)
						{
							//System.out.println("source " + i+" target " + nodeID+" prevID "+ predID +"  w_s "+ w_s +" w_t " + w_t);

						}
					}

					if(predID==5)
					{


						// System.out.println("source " + i+" target " + nodeID+" counted ");

					}

					//Remove nodes with ID prevID from the global sum

					double total_st_removal  = 0;
					double st_r =  0;
					for(int r=0; r< getSize(); r++)
					{
						st_r  = getNode(r).getRiskLevel() - getNode(predID).getRiskLevel();
						if(st_r < 0)
						{
							st_r = 0;
						}

						total_st_removal   = total_st_removal   + st_r;

						st_r  = getNode(predID).getRiskLevel() - getNode(r).getRiskLevel();
						if(st_r < 0)
						{
							st_r = 0;
						}

						total_st_removal   = total_st_removal   + st_r;

					}

					double subtotal_st_diff = total_st_diff - total_st_removal;

					//Calculate the scaling factor;


					double scale = local_st_diff / subtotal_st_diff;

					if(total_st_diff ==0)
					{
						scale = 0;
						System.out.println("WARNING: WEIGHT denominator zero");
					}




					delta[predID] += (double)((sigma[predID]*(1 + delta[nodeID])))/(double)sigma[nodeID];
					scaledDelta[predID] += (double)(scale*(sigma[predID]*(1 + (delta[nodeID])))/(double)sigma[nodeID]);
					//System.out.println("DELTA " + predID+ " FED " + (double)(scale*(sigma[predID]*(1 + delta[nodeID])))/(double)sigma[nodeID]);
					if( (predID ==5) && (nodeID != i))
					{
						System.out.println("source " + i+" target " + nodeID+"  w_s "+ w_s +" w_t " + w_t+" count now "+ delta[predID]);
						System.out.println("scaledcount  now "+ scaledDelta[predID]);
					}
				}

				if(nodeID != i)
				{    if( (nodeID ==5) ) 
				{                	
					System.out.println("count ADD now "+ delta[nodeID]);
					System.out.println("scaledcount ADD now "+ scaledDelta[nodeID]);
				}

				BCentrality[nodeID] += (scaledDelta[nodeID]); //scaled version implemented, but note that no dividing by n-1 or n-2 needed, they have been 'absorbed' earlier.
				int breakpointwatch = 0;
				}
			}


		}

		return BCentrality;
	}



	/**
	 * Calculates the betweenness centrality indices for all of the vertices
	 * in the graph
	 *
	 * <p>
	 * This is the implementation of the algorithm proposed in the research paper
	 * "A Faster Algorithm for Betweenness Centrality" by <br\>
	 * Ulrik Brandes, Department of Computer & Information Science, University of Konstanz
	 * This algorithm calculates the Betweenness Centrality indices for an acyclic
	 * graph in O(mn) time where n is the number of verices and m is the number of
	 * edges 
	 * The memory complexity of the algorithm is O(n^2 + nm) for acyclic graphs
	 * @throws Exception 
	 */
	public boolean[][] returnIndirectAdjacencyArray(boolean[][] adjArray) throws Exception
	{
		adjArray = new boolean[adjArray.length][adjArray.length]; //clear the array!
		int gSize = getSize();

		LinkedBlockingQueue<Node> Q = new LinkedBlockingQueue<Node>();


		for(int i = 0; i < gSize; i++)   // for each vertex in the graph
		{
			//int ithID = (int)nodes.elementAt(i).getID();
			initializeForIndirectAdjacency(i);




			// the entry sigma[t] represents the number shortest pathes from
			// source (denoted by i) to t
			//int[] sigma = new int[gSize];
			//  sigma[i] = 1;


			Q.clear(); //clear the queue

			// put the source node to queue to initialize the queue
			Q.put(getNode(i));

			while(Q.size() > 0)
			{
				Node node = (Node)Q.poll();


				ArrayList<Integer> adjacency = node.getAdjacency();
				for(Integer nodeID : adjacency)
				{
					int neighbourID = nodeID.intValue();
					Node neighbour = getNode(neighbourID);

					// neighbour found for the first time
					if(neighbour.getDistance() == Integer.MAX_VALUE) 
					{
						Q.put(neighbour);
						neighbour.setDistance(node.getDistance() + 1);
					}

					if(neighbour.getDistance() == (node.getDistance() + 1))
					{
						//sigma[neighbourID] += sigma[(int)node.getID()];
						adjArray[i][neighbourID] = true; //sigma[y]; //since sigma is increased, array value will be non-zero: actual value does nto matter 
						//as long as it is non zero
						// predInSP[neighbourID].add(node);
					}
				}
			}

			//here we can fill the 'indirect-adjacency' matrix
			//for(int y=0; y< sigma.length; y++)
			//{
			// adjArray[i][y] = sigma[y];
			//}

			//  int check = 0 ;
			if(i%400==0)
			{
				System.out.println("ii="+i);
			}

		}

		//return BCentrality;
		return adjArray;
	}
	//Code to calculate / get biggest component


	/**
	 * This method calculates and returns the size of the bigggest
	 * connected component in the network.
	 * Network is assumed as UNDIRECTED 
	 *
	 */
	public int getBiggestComponentSize()
	{
		int biggestComponentSize = 0;
		this.sortNodes();

		int[] indexFrequency = new int[getSize()+1];

		Vector<Node> nodeList =  new Vector<Node>();



		//Copy nodes to a new list
		//Note that nodes are NOT CLONED, so each node is still part of the network
		for(int i = 0; i< nodes.size();i++)
		{
			nodeList.addElement(nodes.elementAt(i));
		}

		int componentIndex = 1;

		Node thisNode = nodeList.elementAt(0);
		thisNode.setComponentIndex(componentIndex);
		indexFrequency[componentIndex]++;
		//Now that index is assigned, delete from the list
		//nodeList.removeElementAt(0);

		//Go through all nodes in the list, and if they are neighbours of this node,
		//assign them to the same component

		assignNeighboursToComponent(thisNode, nodeList, 0, componentIndex,indexFrequency);



		//Now go through the rest
		//if they are already assigned a component, number their neighbours

		//if they are not, assign a new number and number the neighbours with the same number
		for(int j = 1; j< nodeList.size();j++)
		{
			thisNode = nodeList.elementAt(j);
			if(thisNode.isComponentAssigned())
			{
				assignNeighboursToComponent(thisNode, nodeList, j, thisNode.getComponentIndex(),indexFrequency);
			}
			else
			{
				componentIndex++;
				thisNode.setComponentIndex(componentIndex);
				indexFrequency[componentIndex]++;
				assignNeighboursToComponent(thisNode, nodeList, j, thisNode.getComponentIndex(),indexFrequency);
			}
		}

		//The above implementation is WRONG, but with nodes sorted in degree order, gives APPROXIMATELY RIGHT RESULTS

		//assume correct numbers have been assigned
		int maxFrequency = 0;
		int biggestComponentIndex = 0;
		for(int j = 0; j< indexFrequency.length;j++)
		{
			if(indexFrequency[j] > maxFrequency)
			{
				maxFrequency = indexFrequency[j];
				biggestComponentSize = maxFrequency;
				biggestComponentIndex = j;
			}
		}

		//    	de assign all components, now that the needed information has been calculated
		for(int j = 1; j< nodeList.size();j++)
		{
			thisNode = nodeList.elementAt(j);
			thisNode.resetComponentAssignedStatus();
		}
		return biggestComponentSize;


	}
	/**
	 * This method calculates and returns the size of the bigggest
	 * connected component in the network.
	 * Network is assumed as UNDIRECTED 
	 * This method uses the recursive call
	 *
	 */

	/**
    public int getBiggestComponentSize_recursive(int[][] paths)
    {
    	int biggestComponentSize = 0;
    	this.sortNodes();

    	int[] indexFrequency = new int[paths.length+1];

    	Vector<Node> nodeList =  new Vector<Node>();



    	//Copy nodes to a new list
    	//Note that nodes are NOT CLONED, so each node is still part of the network
    	for(int i = 0; i< nodes.size();i++)
    	{
    		nodeList.addElement(nodes.elementAt(i));
    	}

    	int componentIndex = 1;

    	Node thisNode;

    	//Now go through the list and assign recursively


    	//if they are not already assigned, assign a new number and number the neighbours with the same number
    	for(int j = 0; j< nodeList.size();j++)
    	{
    		thisNode = nodeList.elementAt(j);
    		//Do a depth search first, to see if any of the close or distance neighbours of this Node has a component assigned
    		//If so, that should be the component id
    		int assignedComp = searchForAssignedNeighbour(thisNode, nodeList, 0);
    		if(assignedComp > 0)
    		{
    			componentIndex = assignedComp;
    		}
    		if(!thisNode.isComponentAssigned())
    		{
    			assignComponentIDToNode(thisNode, nodeList, 0, componentIndex, indexFrequency, paths);
    			int checker = 0;
    		}

    		componentIndex++;

    	}



    	//now correct numbers have been assigned to each node
    	int maxFrequency = 0;
    	int biggestComponentIndex = 0;
    	for(int j = 0; j< indexFrequency.length;j++)
    	{
    		if(indexFrequency[j] > maxFrequency)
    		{
    			maxFrequency = indexFrequency[j];
    			biggestComponentSize = maxFrequency;
    			biggestComponentIndex = j;
    		}
    	}

//    	de-assign all components, now that the needed information has been calculated
    	for(int j = 1; j< nodeList.size();j++)
    	{
    		thisNode = nodeList.elementAt(j);
    		thisNode.resetComponentAssignedStatus();
    	}
    	return biggestComponentSize;


    }
	 */
	/**
	 * This method calculates and returns the size of the bigggest
	 * connected component in the network based on SHORTEST PATHS
	 * Network is assumed as UNDIRECTED 
	 * If there is a shortest path between two nodes, they are indirect neighbours
	 *
	 */
	public int getBiggestComponentSize_sp(boolean[][] pathMatrix)
	{
		int biggestComponentSize = 0;
		this.sortNodes();


		int[] indexFrequency = new int[pathMatrix.length+1];

		Vector<Node> nodeList =  new Vector<Node>();
		int check;


		//Copy nodes to a new list
		//Note that nodes are NOT CLONED, so each node is still part of the network
		for(int i = 0; i< nodes.size();i++)
		{
			nodeList.addElement(nodes.elementAt(i));
		}

		int componentIndex = 1;

		Node thisNode;

		//Now go through the list of nodes



		for(int j = 0; j< nodeList.size();j++)
		{
			thisNode = nodeList.elementAt(j);
			check= thisNode.getNumberOfLinks();
			if(!thisNode.isComponentAssigned())
			{

				thisNode.setComponentIndex(componentIndex);
				indexFrequency[componentIndex]++;
				//System.out.println("node degree: "+thisNode.getNumberOfLinks()+" ID: "+thisNode.getID()+ " assigned TO component = "+ componentIndex);
			}
			else
			{
				continue; //no need to check neighbours
			}
			// see if any of the close or distance neighbours of this Node has a component assigned
			//If so, that should be the component id
			//    		Check all neighbours
			//No need to check neighbours who are alreayd assigned:
			for(int k = j; k< nodeList.size();k++)
			{

				//THis is a neihgbour
				Node potential_neighbour = nodeList.elementAt(k);
				check= potential_neighbour.getNumberOfLinks();

				if(isShortestPath(pathMatrix, (int)thisNode.getID(),(int)potential_neighbour.getID()))
					//these two are indirect neighboursshould be in the same component
				{

					//        						
					if(!potential_neighbour.isComponentAssigned())
					{

						potential_neighbour.setComponentIndex(componentIndex);
						indexFrequency[componentIndex]++;
						//System.out.println("node degree: "+potential_neighbour.getNumberOfLinks()+" ID: "+potential_neighbour.getID()+ " assigned TO component = "+ componentIndex);
					}
					else
					{
						//        				    		/already assigned: do nothing
					}
				}

			}

			componentIndex++;


		}



		//now correct numbers have been assigned to each node
		int maxFrequency = 0;
		int biggestComponentIndex = 0;
		for(int j = 0; j< indexFrequency.length;j++)
		{
			if(indexFrequency[j] > maxFrequency)
			{
				maxFrequency = indexFrequency[j];
				biggestComponentSize = maxFrequency;
				biggestComponentIndex = j;
			}
		}

		//    	de-assign all components, now that the needed information has been calculated
		for(int j = 1; j< nodeList.size();j++)
		{
			thisNode = nodeList.elementAt(j);
			thisNode.resetComponentAssignedStatus();
		}
		return biggestComponentSize;


	}


	public void assignNeighboursToComponent(Node thisNode, Vector<Node>nodeList, int startIndex, int componentNumber, int[] indexFrequency)
	{
		for(int j = startIndex; j< nodeList.size();j++)
		{
			if(thisNode.isLinked(nodeList.elementAt(j)))
			{
				//assign if not assigned already
				if(!nodeList.elementAt(j).isComponentAssigned())
				{
					nodeList.elementAt(j).setComponentIndex(componentNumber);

					indexFrequency[componentNumber]++;
				}
			}
		}
	}

	/*This is a recursice definition
	 * WARNING: function calls itself on all neighbours.
	 * Each neighbour has the responsibility of assignning IDs to all its neighbours 'behind it' (in the 
	 * sorted list) before assignning itself
	 * escape clause: node alreasdy assigned
	 * NEED TO SUPPLY A SHORTEST PATH Matrix of the right size
	 * 
	 */
	public void assignComponentIDToNode(Node thisNode, Vector<Node> nodeList, int startIndex, int componentNumber, int[] indexFrequency, boolean[][] pathMatrix)
	{
		int check = thisNode.getNumberOfLinks();
		check = 0;
		//Check all neighbours
		//No need to check neighbours who are alreayd assigned: hence the use of startindex
		for(int j = 0; j< nodeList.size();j++)
		{
			//recursively assign all neighbours
			if(thisNode.isLinked(nodeList.elementAt(j)) && !((thisNode.getID()== (nodeList.elementAt(j)).getID() )) )
			{
				//THis is a neihgbour
				Node neighbour = nodeList.elementAt(j);
				check = neighbour.getNumberOfLinks();

				if(!nodeList.elementAt(j).isComponentAssigned())
				{

					if(isShortestPath(pathMatrix, (int)thisNode.getID(),(int)neighbour.getID()))
						//these two should be in the same component
					{
						if(nodeList.elementAt(j).getID() < thisNode.getID())
						{
							assignComponentIDToNode(nodeList.elementAt(j),nodeList, 0, componentNumber,indexFrequency, pathMatrix);
						}

					}


				}
				else
				{
					///already assigned: do nothing
				}
			}


		}

		//Now this node can be assigned the ID
		if(!thisNode.isComponentAssigned())
		{
			check = thisNode.getNumberOfLinks();
			thisNode.setComponentIndex(componentNumber);
			indexFrequency[componentNumber]++;
			System.out.println("node degree: "+ check+ " assigned TO component = "+ componentNumber);
		}
	}

	public int searchForAssignedNeighbour(Node thisNode, Vector<Node> nodeList, int startIndex)
	{
		int neighbourAssignedComponentNumber = 0; //A non zero number will mean neighbour has already assigned comp. no
		int check = thisNode.getNumberOfLinks();
		check = 0;
		//Check all neighbours
		//No need to check neighbours who are alreayd assigned: hence the use of startindex
		for(int j = 0; j< nodeList.size();j++)
		{
			//recursively search all neighbours
			if(thisNode.isLinked(nodeList.elementAt(j)) && !((thisNode.getID()== (nodeList.elementAt(j)).getID() )) )
			{
				//THis is a neihgbour
				Node neighbour = nodeList.elementAt(j);
				check = neighbour.getNumberOfLinks();

				if(!nodeList.elementAt(j).isComponentAssigned()) //Not already assigned
				{
					if(nodeList.elementAt(j).getNumberOfLinks() < thisNode.getNumberOfLinks())
					{
						neighbourAssignedComponentNumber = searchForAssignedNeighbour(nodeList.elementAt(j),nodeList, 0);
						if(neighbourAssignedComponentNumber > 0)
						{
							return neighbourAssignedComponentNumber;
						}
					}
					if(nodeList.elementAt(j).getNumberOfLinks() == thisNode.getNumberOfLinks())
					{
						if(nodeList.elementAt(j).getID()>  thisNode.getID())
						{
							neighbourAssignedComponentNumber = searchForAssignedNeighbour(nodeList.elementAt(j),nodeList, 0);
							if(neighbourAssignedComponentNumber > 0)
							{
								return neighbourAssignedComponentNumber;
							}
						}
					}
				}
				else
				{
					//already assigned
					neighbourAssignedComponentNumber = neighbour.getComponentIndex();
					if(neighbourAssignedComponentNumber > 0)
					{
						return neighbourAssignedComponentNumber;
					}
				}
			}


		}

		return neighbourAssignedComponentNumber;

	}

	/**
	 * This method returns TRUE if one or more paths exist between source Node and
	 * target node based on a matrix of number of shortest paths:
	 * the matrix should be supplied as an argument
	 * @param paths
	 * @param sourceID
	 * @param targetID
	 * @return
	 */

	public boolean isShortestPath(boolean[][] paths, int sourceID, int targetID)
	{



		if(paths[sourceID][targetID])
		{
			return true;
		}
		else
		{
			return false;
		}

	}

	public void populate_rho()
	{
		AssortativenessCalculator ac  =  new AssortativenessCalculator(this);
		double rho =0;
		for(int i=0; i< nodes.size();i++)
		{

			nodes.elementAt(i).clearLocalR();
		}

		for(int i=0; i< nodes.size();i++)
		{
			rho = ac.calculate_localR(nodes.elementAt(i));
			//rho = calculate_localR(nodes.elementAt(i));
			nodes.elementAt(i).setLocalR(rho);
		}
	}

	/**
	 * This calculates betweenness centrality in an innefficient way -  as opposed to the Brandes algorithm which is an efficient way.
	 *
	 */
	void populate_Namal_BC() 
	{

		/* size is incremented by 1 so that node's can start from 1 rather
          than 0, this eases the operation */
		int gSize = getSize();   	



		float[] BCDistribution = new float[gSize];

		for(int i = 0; i < gSize; i++)
		{

			try {
				calculateShortestPath(i);
				updateAllPairSPArrays(i, BCDistribution);

			} catch (Exception e) {
				e.printStackTrace();
			}

		}  


		float div = (float)(gSize - 1)*(gSize - 2);
		for(int v = 0; v < gSize; v++)
		{
			System.out.println(v + "\t" + (BCDistribution[v]/div));

			nodes.elementAt(v).setBCCentrality(BCDistribution[v]/div);

		}


	}

	/**
	 * This calculates Risk centrality in an innefficient way -  as opposed to the Brandes algorithm which is an efficient way.
	 *
	 */
	void populate_Namal_RC() 
	{

		/* size is incremented by 1 so that node's can start from 1 rather
          than 0, this eases the operation */


		//First calculates the 'total weights' that would be needed
		double total_riskweight = 0;
		for(int p=0; p< getSize(); p++)
		{
			try {
				total_riskweight = total_riskweight + getNode(p).getRiskLevel();
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}
		}




		int gSize = getSize();  

		float[] BCDistribution = new float[gSize];

		for(int i = 0; i < gSize; i++)
		{

			try {
				calculateShortestPath(i);
				updateAllPairSPArrays_rc(i, BCDistribution, total_riskweight);

			} catch (Exception e) {
				e.printStackTrace();
			}

		}  


		float div = (float)(gSize - 2);
		for(int v = 0; v < gSize; v++)
		{
			System.out.println(v + "\t" + (BCDistribution[v]/div));

			nodes.elementAt(v).setBCCentrality(BCDistribution[v]/div);

		}


	}

	/**
	 * This calculates Risk centrality for source -  target weighted version. This cannot be implemented as a modification of Brande's algorithm. The straght forward
	 * implementation of betweenness centrality must be modified.
	 * The TIME COMPLEXITY IS O(N^3). THIS IS STILL IMPRESSIVE.
	 */
	void populate_Namal_RC_st() 
	{



		//First calculates the 'total weights' that would be needed
		//first calculate total of Tramp(sw - tw) this will be needed for scaling.
		double total_st_diff = 0;
		double st = 0;
		for(int s=0; s< getSize(); s++)
		{
			for(int t=0; t< getSize(); t++)
			{
				try {
					st  = getNode(s).getRiskLevel() - getNode(t).getRiskLevel();
					if(st < 0)
					{
						st = 0;
					}

					total_st_diff  = total_st_diff  + st;

				} catch (Exception e) {
					System.out.println("node with given ID could not be found in sum of weight differences (s-t) calculation");
					e.printStackTrace();
				}

			}

		}



		int gSize = getSize();  

		float[] BCDistribution = new float[gSize];

		for(int i = 0; i < gSize; i++)
		{

			try {
				calculateShortestPath(i);
				updateAllPairSPArrays_rc_st(i, BCDistribution, total_st_diff);

			} catch (Exception e) {
				e.printStackTrace();
			}

		}  


		float div = (float)1.0; //No division needed for s-t risk centrality
		for(int v = 0; v < gSize; v++)
		{
			System.out.println(v + "\t" + (BCDistribution[v]/div));

			nodes.elementAt(v).setBCCentrality(BCDistribution[v]/div);

		}


	}




	public void populate_BC()
	{
		float[] bc;
		try {
			bc = this.calcBetweenessCentrality();


			for(int i=0; i< nodes.size();i++)
			{

				nodes.elementAt(i).clearBCCentrality();
			}

			for(int i=0; i< nodes.size();i++)
			{
				int thisID =  (int) nodes.elementAt(i).getID();
				double bcVal = bc[thisID];
				nodes.elementAt(i).setBCCentrality(bcVal);
			}

		} catch (Exception e) {
			System.out.println("problem calculating betweennes array");
		}

		//System.out.println("finished pop bc");
	}

	public void populate_CC()
	{
		float[] cc;
		try {
			cc = this.calcClosenessCentrality();


			for(int i=0; i< nodes.size();i++)
			{

				nodes.elementAt(i).clearCCCentrality();
			}

			for(int i=0; i< nodes.size();i++)
			{
				int thisID =  (int) nodes.elementAt(i).getID();
				double ccVal = cc[thisID];
				nodes.elementAt(i).setCCCentrality(ccVal);
			}

		} catch (Exception e) {
			System.out.println("problem calculating CLOENESS array");
		}

		//System.out.println("finished pop bc");
	}

	public void populate_RC()
	{
		float[] bc;
		try {
			bc = this.calcRiskCentrality();


			for(int i=0; i< nodes.size();i++)
			{

				nodes.elementAt(i).clearRCCentrality();
			}

			for(int i=0; i< nodes.size();i++)
			{
				int thisID =  (int) nodes.elementAt(i).getID();
				double bcVal = bc[thisID];
				nodes.elementAt(i).setRCCentrality(bcVal);
			}

		} catch (Exception e) {
			System.out.println("problem calculating betweennes array");
		}

		//System.out.println("finished pop rc");
	}
	/**
	 * A method to populate hop distance of all nodes in network
	 * A simple implementation with O(N^2) time
	 *
	 */

	public void populate_HopDistance()
	{



		for(int i=0; i< getSize();i++)
		{

			Node thisNode = getAllNodes().elementAt(i);
			int hd = 100000; //VERY BIG

			for(int j=0; j< getSize();j++)
			{

				Node otherNode = getAllNodes().elementAt(i);
			}


		}

		System.out.println("WARNING: CALCULATING HOP DISTANCES IS NOT FULLY IMPLEMENETED");
	}


	/**
	 * count and print the number of infected nodes according to a double threhold
	 */
	public int countInfectedNodes(double threshold)
	{
		int count = 0;
		for(int i=0; i<getSize(); i++)
		{
			if(getAllNodes().elementAt(i).getRiskLevel()>=threshold)
			{
				count++;
			}
		}


		return count;
	}

	/**
	 * This function reads sensor data from a file and returns it as a matrix
	 * @return
	 */
	public  double[][] readMatrix()
	{


		BufferedReader stdrd = null;
		BufferedReader counter = null;
		int lines = 0;
		int columns =0;

		String filename = ".\\SensorAssort\\caltech-facstatus.txt";

		//first count lines and columns
		try {

			counter = new BufferedReader(new FileReader(filename));


			String str  = counter.readLine();

			while(!(str==null))
			{

				if(str==null)
				{
					break;
				}

				if(str.equals(""))
				{
					break;
				}
				//There are lines to read
				lines++;
				if(lines <=1)
				{
					StringTokenizer st = new StringTokenizer(str);
					while (st.hasMoreTokens()) {
						columns++;
						st.nextToken();

					}
					// columns++;

					int  x = 1;
				}

				str  = counter.readLine();

			}


		}

		catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}

		lines = lines+10;
		int iLimit = 1;
		//Now the lines and colums have been counted
		double[][] matrix = new double[lines][iLimit];
		try {

			stdrd = new BufferedReader(new FileReader(filename));

			int linecount = 0;
			String str  = stdrd.readLine();

			while(!(str==null))
			{

				if(str==null)
				{
					break;
				}

				if(str.equals(""))
				{
					break;
				}
				//There are lines to read

				int i = 0;
				StringTokenizer st = new StringTokenizer(str);
				while (st.hasMoreTokens()) {
					if(i >= iLimit)
					{
						break;
					}

					//System.out.println(st.nextToken());
					double value = Double.parseDouble(st.nextToken());
					matrix[linecount][i] = Math.round(value);
					i++;

				}
				linecount++;
				str  = stdrd.readLine();

			}



		}

		catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}

		return matrix;
	}

	/**
	 * This function assigns sensor values at a given timestep to nodes
	 * @param matrix
	 */
	public void assignSensorValues(double[][] matrix, int timeStep)
	{
		for(int i=0; i< this.getSize(); i++)
		{
			Node thisNode =  this.getAllNodes().elementAt(i);
			thisNode.setReading(matrix[(int)thisNode.getID()][timeStep]);
			System.out.println("node "+ thisNode.getID() + " assighned "+matrix[(int)thisNode.getID()][timeStep]);

		}
	}

	/**
	 * This method updates sensor values to a pattern where
	 * nodes get positive/negative/random feedback from their neighbours
	 * 
	 */

	public void simNet(int timeSteps)
	{

		//prepare for printing

		BufferedWriter out = null;



		try {

			out = new BufferedWriter(new FileWriter(".\\..\\Output\\CONG\\cong e EXPRESSION DUMP.xls"));

		}

		catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}


		// print IDs

		this.sortNodes();


		try
		{
			for(int j=0; j < nodes.size(); j++)
			{
				Node thisNode = nodes.elementAt(j);

				out.write(thisNode.getID()+"\t");
				out.flush();

			}

			out.write( System.getProperty(("line.separator")));
			out.flush();
		}

		catch (IOException e)
		{
			System.out.println("cant write to file ");
			e.printStackTrace();
		}



		//similate network




		for (int i=0; i< timeSteps; i++)
		{
			//simNetStep_logic2();
			simNetStep_logicPageRank();
			//assignAllOnes();
			//one time step completed

			int y = i; //not needed except for debudding

			//print node states
			this.sortNodes(); //Make sure to sort based on ID
			for(int j=0; j < nodes.size(); j++)
			{
				Node thisNode = nodes.elementAt(j);


				try
				{
					out.write(thisNode.getReading()+"\t");
					out.flush();


				}
				catch (IOException e)
				{
					System.out.println("cant write to file ");
					e.printStackTrace();
				}

			}

			try
			{
				out.write( System.getProperty(("line.separator")));
				out.flush();
			}

			catch (IOException e)
			{
				System.out.println("cant write to file ");
				e.printStackTrace();
			}

			//printing finished

		}



	}

	//Assign note states to random binary values

	public void assignRandomBinaries()
	{


		//Make random asignment at the start
		for(int i=0; i< this.getSize(); i++)
		{
			Node thisNode = this.getAllNodes().elementAt(i);
			double state = 0;


			if(random.nextBoolean())
			{
				state = 1;
			}

			//	state =  random.nextDouble();

			thisNode.setReading(state);

		}
	}


	public void assignRandomBinaries_1(double oneProb)
	{


		//Make random asignment at the start
		for(int i=0; i< this.getSize(); i++)
		{
			Node thisNode = this.getAllNodes().elementAt(i);
			double state = 0;


			if(random.nextDouble()< oneProb)
			{
				state = 1;
			}



			thisNode.setReading(state);

		}
	}

	public void assignAllOnes()
	{


		for(int i=0; i< this.getSize(); i++)
		{
			Node thisNode = this.getAllNodes().elementAt(i);
			double state = 1;
			thisNode.setReading(state);

		}
	}

	public void assignAllZeros()
	{


		for(int i=0; i< this.getSize(); i++)
		{
			Node thisNode = this.getAllNodes().elementAt(i);
			double state = 0;
			thisNode.setReading(state);

		}
	}

	/**
	 * Simulate for a time step
	 *
	 */
	public void simNetStep_logic1()
	{
		this.shuffle(); //shuffle nodes

		for(int j=0; j < nodes.size(); j++)
		{
			//select nodes in random order

			Node thisNode = nodes.elementAt(j);

			double avgNeighbourReading =  thisNode.getAverageNeighbourReading();
			if(avgNeighbourReading >= 0.5) //neighbours are largely 'positive'
			{
				if(thisNode.getReading() < 0.5 )
				{
					//This node has negative reading


					//Swap with the 'weakest' (in terms of degree) neighbour

					Node sn =  thisNode.getSmallestNeighbour();

					double temp = thisNode.getReading();

					thisNode.setReading(sn.getReading());
					sn.setReading(temp);
				}
			}

		}

	}

	/**
	 * Another logic
	 */

	/**
	 * Simulate for a time step
	 * Tries to FOLLOW neighbour states
	 *
	 */
	public void simNetStep_logic2()
	{
		this.shuffle(); //shuffle nodes

		for(int j=0; j < nodes.size(); j++)
		{
			//select nodes in random order

			Node thisNode = nodes.elementAt(j);

			double avgNeighbourReading =  thisNode.getAverageNeighbourReading();
			double changeProb = 0;

			if(avgNeighbourReading >= 0.5) //neighbours are largely 'positive'
			{
				changeProb =  0.25*(avgNeighbourReading - 0.5);
				if(thisNode.getReading() < 0.5 )
				{
					//With a probability, set it to '1'
					if(random.nextDouble() < changeProb)
					{
						thisNode.setReading(1.0);
					}
				}
			}

			else //neighbours are largely 'negative'
			{
				changeProb = 0.25*(0.5 -  avgNeighbourReading);
				if(thisNode.getReading() > 0.5 )
				{
					//With a probability, set it to '0'
					if(random.nextDouble() < changeProb)
					{
						thisNode.setReading(0.0);
					}
				}
			}

			//Logic completed



		}

	}

	/**
	 * Another logic
	 */

	/**
	 * Simulate for a time step
	 *
	 */
	public void simNetStep_logic3()
	{
		this.shuffle(); //shuffle nodes

		for(int j=0; j < nodes.size(); j++)
		{
			//select nodes in random order

			Node thisNode = nodes.elementAt(j);

			double nodeStrength = (double)thisNode.getNumberOfLinks() / (double)this.getHighestDegree(); //Node's probability
			//nodeStrength = nodeStrength + 0.5;
			//to influence change is prop. to its degree

			double avgNeighbourReading =  thisNode.getAverageNeighbourReading();
			double changeProb = 0;

			if(avgNeighbourReading >= 0.5) //neighbours are largely 'positive'
			{
				changeProb =  -nodeStrength+(avgNeighbourReading - 0.5);
				if(thisNode.getReading() < 0.5 )
				{
					//With a probability, set it to '1'
					if(random.nextDouble() < changeProb)
					{
						thisNode.setReading(0.0);
					}
				}
			}

			else //neighbours are largely 'negative'
			{
				changeProb = -nodeStrength+(0.5 -  avgNeighbourReading);
				if(thisNode.getReading() > 0.5 )
				{
					//With a probability, set it to '0'
					if(random.nextDouble() < changeProb)
					{
						thisNode.setReading(1.0);
					}
				}
			}

			//Logic completed



		}

	}


	public void simNetStep_logic4()
	{
		this.shuffle(); //shuffle nodes

		for(int j=0; j < nodes.size(); j++)
		{
			//select nodes in random order

			Node thisNode = nodes.elementAt(j);

			Node bN = thisNode.getBiggestNeighbour();

			if(thisNode.getID()== bN.getID())
			{
				//no need to do anything
				break;
			}

			else
			{
				//probabilistically change to the opposite of biggest neighbour
				if(random.nextDouble() < 0.45)
				{
					thisNode.setReading((int)((bN.getReading()+1.0))%2);
				}
			}

			//Logic completed



		}

	}


	public void simNetStep_logicPageRank()
	{
		double alpha =  0.8;
		this.shuffle(); //shuffle nodes
		double ave = this.updatePageRank(alpha);

		for(int j=0; j < nodes.size(); j++)
		{
			Node thisNode = nodes.elementAt(j);
			//thisNode.setReading(thisNode.getPageRank());
			if( (ave -  thisNode.getPageRank()) > 0)
			{
				thisNode.setReading(1);

			}
			else
			{
				thisNode.setReading(0);
			}
		}
	}
	/**
	 * Simply counts nodes whose status is '1' (relevant only in boolean context)
	 * @return
	 */

	public int  countOnes()
	{
		int maxReading = 1000;
		double[] dist = new double[maxReading];
		int count = 0;

		for(int i=0;i<nodes.size();i++) 
		{			


			if((int)nodes.elementAt(i).getReading()==1)
			{
				count++;
			}
		}

		return count;
	}


	public double calculate_MI_straight_sensor() {

		double[][] ejk = this.getPjkDistribution();
		double[] qk = this.getqkDistribution_sensor();
		double i = 0;

		for (int j = 0; j < qk.length; j++)
		{
			for (int k = 0; k < qk.length; k++)
			{  double addition =  ejk[j][k] / (qk[j]*qk[k]);
			if(addition != 0)
			{
				if( (qk[j] != 0) && (qk[k] != 0) )
				{
					i = i + (ejk[j][k])*Math.log(  ejk[j][k] / (qk[j]*qk[k])  );
				}
			}
			}

		}
		i  =  i /Math.log(2);

		return i;
	}





	public double calculate_H1_straight_sensor() {


		double[] qk = this.getpkDistribution_sensor();
		double i = 0;

		for (int j = 0; j < qk.length; j++)
		{

			if( (qk[j] != 0) )
			{
				i = i + Math.log(  qk[j]  );
			}



		}
		i  =  i /Math.log(2);

		return -i;
	}

	public double calculate_H2_straight_sensor() {


		double[] qk = this.getqkDistribution_sensor();
		double i = 0;

		for (int j = 0; j < qk.length; j++)
		{

			if( (qk[j] != 0) )
			{
				i = i + Math.log(  qk[j]  );
			}



		}
		i  =  i /Math.log(2);

		return -i;
	}

	/**
	 * Delete the linkless nodes from network
	 */

	public void pruneNetwork()
	{
		for(int k=0; k< nodes.size(); k++)
		{
			if(((Node)(nodes.elementAt(k))).getNumberOfLinks()==0) //link less node
			{
				try {
					removeNode((int)  ( (nodes.elementAt(k)).getID()) );
				} catch (Exception e) {
					System.out.println("Trying to delete node with non existing ID");
					e.printStackTrace();
				}
			}
		}
	}


	//==================================================================================

	//Code from this point onward is for directed Dynamic networks

	//TODO: All the code below should be checked.

	/**The following method returns e(J_tcurrent, K_tprevious), where J,K are remaining degree
	 * The neighbours remaining degree is measured from the 'previous' degree of the neighbour
	 * All nodes in the network must belog to the subclass 'TimeNodes', which remembers the previous degree in  a dynamic network
	 * This is a 'directed' distribution, ( though e(J,K) is normally undirected ) with K being the neighbour.
	 */

	public double[][]  get_ejk_conditional_on_time()
	{
		double[][] dist = new double[maxDegree][maxDegree];

		for(int i=0;i<links.size();i++) 
		{	
			int a1 = links.elementAt(i).iNode().getNumberOfLinks()-1;
			int b1 = (int)((TimeNode)(links.elementAt(i).jNode())).getPrevDegree()-1;

			int a2 = links.elementAt(i).jNode().getNumberOfLinks()-1;

			int b2 = (int)((TimeNode)(links.elementAt(i).iNode())).getPrevDegree()-1;


			if(a1<0)
				a1 = 0;

			if(a2<0)
				a2 = 0;

			if(b1<0)
				b1 = 0;

			if(b2<0)
				b2 = 0;

			dist[a1][b1]= (0.5/links.size())+  dist[a1][b1];		
			dist[a2][b2]= (0.5/links.size())+  dist[a2][b2];		


		}

		return dist;
	}

	/**
	 * Remaining degre distribution for a node at previous timestep in a dynamic network
	 */
	public double[]  get_qk_conditional_on_time()
	{
		double[] dist = new double[maxDegree];
		for(int i=0;i<links.size();i++) 
		{

			//int a1 = links.elementAt(i).iNode().getNumberOfLinks()-1;
			int b1 = (int)((TimeNode)(links.elementAt(i).jNode())).getPrevDegree()-1;

			//int a2 = links.elementAt(i).jNode().getNumberOfLinks()-1;

			int b2 = (int)((TimeNode)(links.elementAt(i).iNode())).getPrevDegree()-1;


			//if(a1<0)
			//a1 = 0;

			//	if(a2<0)
			//a2 = 0;

			if(b1<0)
				b1 = 0;

			if(b2<0)
				b2 = 0;

			dist[b2]= (0.5/(double)(1.0*links.size()))+dist[b2];
			dist[b1]= (0.5/(double)(1.0*links.size()))+dist[b1];
		}

		return dist;
	}

	/**
	 * Calculate the conditional Entropy of of dynamic network
	 * Answer in Bits
	 */
	public double calculate_Hqc_conditional_on_time()
	{
		double in =  0;
		double[][] ejk_t = get_ejk_conditional_on_time();
		double[] q_t = get_qk_conditional_on_time();
		for (int i = 0; i < ejk_t.length; i++) {
			for (int j = 0; j < ejk_t[i].length; j++) {
				if ( (q_t[j] != 0) &&  (ejk_t[i][j] != 0))
				{
					in = in - (  ejk_t[i][j] * Math.log(ejk_t[i][j] / q_t[j])  ); 
				}
			}
		}

		return ( in / Math.log(2)); /* Convert to bits before returning */
	}


	/**
	 * Calculate the conditional Entropy of of STATIC network based on Formula 9, Sole
	 * Answer in Bits
	 */
	public double calculate_Hqc_2()
	{
		double in =  0;
		double[][] ejk_t = getLinkDistribution();
		double[] q_t = getRemainingDegreeDistribution();
		for (int i = 0; i < ejk_t.length; i++) {
			for (int j = 0; j < ejk_t[i].length; j++) {
				if ( (q_t[j] != 0) &&  (ejk_t[i][j] != 0))
				{
					in = in - (  ejk_t[i][j] * Math.log(ejk_t[i][j] / q_t[j])  ); 
				}
			}
		}

		return ( in / Math.log(2)); /* Convert to bits before returning */
	}

	public double calculate_MI_Dynamic()
	{
		return ( calculateHq() - calculate_Hqc_conditional_on_time() );
	}


	/**
	 * Calculate the Entropy of this Q distribution, H(Q)for prev. time step
	 * Answer in Bits
	 */
	public double calculateHq_prev() {
		double h = 0;
		double[] qd = get_qk_conditional_on_time();
		for (int i = 0; i < qd.length; i++) {
			if (qd[i] != 0)
			{
				h = h - (qd[i] * Math.log(qd[i]));
			}
		}
		return (h/Math.log(2)); /* Convert to bits befreo returning */
	}


	/**
	 * This method implements the PARG model 
	 * @param net
	 * @param finalSize
	 * @return
	 */
	public void growPARGNet(int finalSize)
	{

		int initialSize = this.getSize();
		System.out.println("Number of Nodes now  = " + this.getSize());
		System.out.println("Number of links now  = " + this.getAllLinks().size());
		Random rand  =  new Random(System.currentTimeMillis());
		//	
		double EjoinSTART = 1.0000;
		double EdeleteSTART =1.002;
		double Ejoin = EjoinSTART; //Expectation of links created per node- a parameter
		double Edelete = EdeleteSTART; // Expected links deleted per node addition. A parameter

		int kcounter = 0;

		for(int i = 1; i <= finalSize; i++)
		{

			if(i%1000==0)
			{
				kcounter++;
				System.out.println(kcounter +" thousand Nodes grown");
			}

			Ejoin = EjoinSTART; // + (0.02*(double)(net.getSize())); //increase addded links with growing network size
			Edelete = EdeleteSTART; // + (0.02*(double)(net.getSize())); //increase addded links with growing network size
			//Node node = new Node(initialSize + i);
			TimeNode node = new TimeNode(initialSize + i);
			double probs[] =  new double[this.getSize()];
			int noOfLinks = this.getNoOfLinks();

			for(int j = 0; j <this.getSize();j++)
			{

				//Node destNode = this.getAllNodes().elementAt(j);
				TimeNode destNode = (TimeNode)this.getAllNodes().elementAt(j);
				double prob = Ejoin;
				if(this.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
				{
					prob = 0.5*Ejoin*((double)destNode.getNumberOfLinks()) / (double)this.getNoOfLinks();
				}
				else //No links made yet
				{
					prob = 1.0;//If no links are made, this node node has to join with the existing isolated nodes
					//isolated nodes are not allowed to exist.

				}
				probs[j] = prob;
			}


			for(int j = 0; j <this.getSize();j++)
			{
				//Node destNode = this.getAllNodes().elementAt(j);
				TimeNode destNode = (TimeNode)this.getAllNodes().elementAt(j);
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

			noOfLinks = this.getNoOfLinks();
			//System.out.println("Number of links now  = " + net.getAllLinks().size());
			int linksNow = this.getAllLinks().size();

			//Now, stochastically delete some links that connect hubs
			//			Now, stochastically delete some links that connect hubs

			double raw =0.99999;


			//If this is set to 1.0, stochastic deleting is ALWAYS done.



			//net.sortLinksBasedOnCombinedDegree();
			this.sortNodes();

			if( rand.nextDouble() <=  raw ) 
			{
				double probsum = 0;
				//destroy 'Edeleted' number of the existing links
				double probs2[] =  new double[this.getNoOfLinks()];
				for(int k = 0; k <this.getNoOfLinks();k++)
				{

					Link destLink = this.getAllLinks().elementAt(k);
					double prob = ((double)Edelete)*0.25*0.01*(double)this.getNoOfLinks()*((double)((destLink.iNodeLinks() )+( destLink.jNodeLinks())) / (double)this.getNoOfLinks());
					// prob = ((double)Edelete)*0.25*((double)((destLink.iNodeLinks() )+( destLink.jNodeLinks())) / (double)net.getNoOfLinks());
					//More aggresive method has to be implemented here to selectively delete hub-to-hub links
					/*
					 prob = 0;
					 int Np = net.getHighestDegree();
					 if((destLink.iNodeLinks()+3 > Np )&& (destLink.jNodeLinks()+8 > Np) )
					 {
						 prob = 1;
					 }
					 if((destLink.iNodeLinks()+7 > Np)&& (destLink.jNodeLinks()+8 > Np) )
					 {
						 pro
						 }
					 */
					if(k >= (this.getNoOfLinks() -  (20.0*Edelete))) //This is a BIG (hub-hub) link, should be deleted.
						//ASSORTATIVE conenctions should be discouraged here; BOTH hub-to-hub as well as periphery-to-periphery
						//1+2+3+4+....+N = N*0.5*(N-1)
					{
						//double rank = net.getNoOfLinks() - k;
						//double weight = (20.0 - rank+ 1.0)/ (10.0*21.0);
						//prob = weight;
						//prob = Edelete*2.0*(double)k / ((double)net.getNoOfLinks()*(double)net.getNoOfLinks()-1.0);
						//prob = 1.0;
						//probsum = probsum + prob;
					}
					else
					{
						//prob = 0;
						//prob = Edelete*2.0*(double)k / ((double)net.getNoOfLinks()*(double)net.getNoOfLinks()-1.0);
					}
					//prob = 0;

					TimeNode nodei = (TimeNode)destLink.iNode();
					TimeNode nodej = (TimeNode)destLink.jNode();
					TimeNode hubNode;
					TimeNode nonHubNode;
					int bestPosition;

					int iPosition = this.getPosition((int)nodei.getID());
					int jPosition = this.getPosition((int)nodej.getID());

					if(iPosition < jPosition)
					{
						bestPosition  = iPosition;
						hubNode = nodei;
						nonHubNode = nodej;
					}
					else
					{

						bestPosition  = jPosition;
						hubNode = nodej;
						nonHubNode = nodei;
					}

					double neighbourSum = 0;
					//Calculate sum of degreed for all neighbours of hub
					for(int y=0; y < hubNode.getNumberOfLinks(); y++)
					{
						Link thisLink = hubNode.getLinks().elementAt(y);

						if(thisLink.iNodeID() == hubNode.getID())
						{
							neighbourSum = neighbourSum + thisLink.jNode().getNumberOfLinks();
						}
						else //jNode is the HubNode
						{
							neighbourSum = neighbourSum + thisLink.iNode().getNumberOfLinks();
						}
					}


					double CUTOFF = 1.0 + (0.0065*this.getSize());

					if(bestPosition < CUTOFF) //This Link has a node that is one of the hubs
					{
						double pp = (1.0 / (double)hubNode.getNumberOfLinks()); // one over number of Links for that Node
						//double w = 200.0*(double) nonHubNode.getNumberOfLinks() / neighbourSum;
						//System.out.println("w "+w);
						prob = (1.0 / (double)CUTOFF)* pp * Edelete; // Stochastically delete 'Edelete' number of links from the hub

					}
					else
					{
						prob = 0;
					}


					//If this Link has a Node which is one of the hubs, should be deleted stochastically



					probs2[k] = prob;

					//randomly select some links from the 'assortative links' and replace them
				}



				for(int k = 0; k <this.getNoOfLinks();k++)
				{
					Link destLink = this.getAllLinks().elementAt(k);
					//Throw dice
					if(rand.nextDouble() <=  probs2[k])
					{
						try 
						{
							//Delete Link

							//Deleting immediately from network will cause problems with size, therefore mark and delete later.
							this.getAllLinks().elementAt(k).markForDeletion();
							Node iNode = this.getAllLinks().elementAt(k).iNode();
							Node jNode = this.getAllLinks().elementAt(k).jNode();
							//System.out.println("Link REPLACED "+ iNode.getNumberOfLinks()+ " "+jNode.getNumberOfLinks());
							//System.out.println("BETWEEN NODES "+ iNode.getID()+ " "+jNode.getID());
							iNode.deleteLink(jNode);
							jNode.deleteLink(iNode);


						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}


				//Now delete all the marked links
				//net.removeMarkedLinks();
				this.replaceMarkedLinks(); //replace deleted links by connecting the corresponding nodes to other nodes selected preferentially,
				//so that a deleted Link does not result in breaking the network into components
				noOfLinks = this.getNoOfLinks();
				int la = 0;

				//Now the links have been deleted and second mechanism completed
				//System.out.println("Number of Nodes now  = " + net.getSize());
				//System.out.println("Number of links now  = " + net.getAllLinks().size());
				//System.out.println("DELETED Links in this node addition  = " + (linksNow - net.getAllLinks().size()));

			}


		}

		System.out.println("Number of Nodes now  = " + this.getSize());
		System.out.println("Number of links now  = " + this.getAllLinks().size());

		//return net;
	}

	/**
	 * This method randomly grows a network without changing the total number of nodes.
	 * The links are re assigned
	 *
	 */
	public void  growNetwork()
	{
		int ranChanges = 500; //The number of random changes to be made
		int linksno = this.getNoOfLinks();


		int linkSize1 = getNoOfLinks();

		for(int i=0; i < ranChanges; i++)
		{
			int l1 = Math.abs(random.nextInt())%linksno;
			int l2 = Math.abs(random.nextInt())%linksno;
			Link L1 = this.getAllLinks().elementAt(l1);
			//L1.markForDeletion();
		}

		this.removeMarkedLinks();

		int linkSize2 = getNoOfLinks();

		//for(int i=0; i < (linkSize1 -  linkSize2); i++)
		for(int i=0; i < ranChanges; i++)
		{
			int l1 = Math.abs(random.nextInt())%getSize();
			int l2 = Math.abs(random.nextInt())%getSize();


			try {
				this.addLink(l1, l2);
			} catch (Exception e) {
				System.out.println("Node with given ID  does not exist in growNetwork method");
				e.printStackTrace();
			}
		}

		System.out.println("Number of Nodes now  = " + this.getSize());
		System.out.println("Number of links now  = " + this.getAllLinks().size());


	}




	/**
	 * Grow a network of a given size using preferential attachment
	 * @param finalSize denotes the size of the network after growth 
	 * @return
	 */
	public void growNetworkPA(int finalSize, double joiningp)
	{
		BufferedWriter out = null;


		Random rand = new Random();


		TimeNode firstNode  = new TimeNode(getSize());
		firstNode.assignWeight(rand.nextDouble());
		this.addNode(firstNode);

		//Add the rest of the nodes acccording to PA
		double p = joiningp; //Joining probability - a parameter
		for(int i = 1; i < finalSize -  getSize(); i++)
		{
			TimeNode node = new TimeNode(i);
			node.assignWeight(rand.nextDouble());
			double probs[] =  new double[this.getSize()];

			for(int j = 0; j <this.getSize();j++)
			{
				TimeNode destNode = (TimeNode) (this.getAllNodes().elementAt(j));
				double prob = p;
				if(this.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
				{
					prob = p*((double)destNode.getNumberOfLinks()) / (double)this.getNoOfLinks();
					// prob = (prob*destNode.getWeight()* (double)net.getNoOfLinks()) / (double)net.getWeightedDegreeSum(); //For weighted PA
				}
				else //No links made yet
				{
					prob = p;
					//prob = prob*destNode.getWeight(); //For weighted PA
				}
				probs[j] = prob;
			}

			for(int j = 0; j <this.getSize();j++)
			{
				TimeNode destNode = (TimeNode) this.getAllNodes().elementAt(j);
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
		System.out.println("net size nodes "+this.getSize()+ " links "+ this.getNoOfLinks());

	}

	/**
	 * This method implements the interactive growth model proposed by Zhou
	 * @param network
	 * @param Mrand
	 * @return
	 */
	public void growIGNet(int size)
	{
		//start with a random net of m0 nbodes and links
		//we take that m0 = 10

		int m0 = getSize();

		Random rand  =  new Random(System.currentTimeMillis());	

		for(int i = (m0+1) ;i <= size ; i++)
		{


			if ( rand.nextDouble() <= 0.4) //with 40% probability
			{
				TimeNode node = new TimeNode(i); 
				this.addNode(node);

				//connnect it to one host node
				int linkIDMax = this.getNoOfLinks();
				int NodeIDMax = this.getSize()-2;
				long hostID = (this.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				//linear preference is implemented by randomizing link based
				//int hostID = 1 + (Math.abs((rand.nextInt()))%NodeIDMax); //implement linear preference
				try {
					this.addLink(i, (int)hostID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}

				//coonect the host to twp peers
				long peer1ID = (this.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				long peer2ID = (this.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 


				try {
					this.addLink((int)peer1ID, (int)hostID);
					this.addLink((int)peer2ID, (int)hostID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
			}

			else //with 60 % probability
			{
				TimeNode node = new TimeNode(i); 
				this.addNode(node);

				//connnect it to two host node
				int NodeIDMax = this.getSize()-2;
				int linkIDMax = this.getNoOfLinks();
				long host1ID = (this.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				long host2ID = (this.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 



				try {
					this.addLink(i, (int)host1ID);
					this.addLink(i, (int)host2ID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}



				//Connect one of the host nodes to peers
				//Here we choose the first hhost node, since both host nodes were randomly selected
				long peerID = (this.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 


				try {
					this.addLink((int)peerID,(int) host1ID);

				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}

			}

			//keep doing this until we get the required size
		}

		System.out.println("Number of Nodes now  = " + this.getSize());
		System.out.println("Number of links now  = " + this.getAllLinks().size());

	}

	public void initialize_pagerank(double r)
	{
		for(int i=0; i< this.getSize(); i++)
		{
			Node thisNode =  nodes.elementAt(i);
			thisNode.p_t = r;
			thisNode.p_t_1 = r;

			System.out.println("Node "+ thisNode.getID()+" Rank "+ thisNode.p_t);
		}
	}

	public void initialize_pagerank_random()
	{

		for(int i=0; i< this.getSize(); i++)
		{
			double r =  random.nextDouble();
			Node thisNode =  nodes.elementAt(i);
			thisNode.p_t = r;
			thisNode.p_t_1 = r;

			System.out.println("Node "+ thisNode.getID()+" Rank "+ thisNode.p_t);
		}
	}

	public double updatePageRank(double alpha)
	{
		double average = 0;
		double median = 0;
		double[] scores =  new double[this.getSize()];
		for(int i=0; i< this.getSize(); i++)
		{
			Node thisNode =  nodes.elementAt(i);
			thisNode.updatePageRank(alpha,this.getSize());
		}

		//Only after all p_t are updated, we must update p_t-1, to avoid synchronization problems

		for(int i=0; i< this.getSize(); i++)
		{
			Node thisNode =  nodes.elementAt(i);
			thisNode.p_t_1 =  thisNode.p_t;
			//System.out.println("Node "+ thisNode.getID()+" Rank "+ thisNode.p_t);
			average = average + thisNode.p_t_1;
			scores[i] = thisNode.p_t_1;
		}

		Arrays.sort(scores);
		int middleIndex =  this.getSize()/2;

		average =  average/(double)this.getSize();
		median  =  scores[middleIndex];

		System.out.println("average"+ average+" median "+ median);
		//return average;
		return median;

	}

	/**
	 * Method regarding new node assortativity definitions
	 */
	public double sumOfDDAverage()
	{
		double S = 0 ;
		for(int i=0; i< this.getSize(); i++)
		{
			Node thisNode =  nodes.elementAt(i);
			S =  S +thisNode.getAverageNeighbourDifference();
		}


		return S;
	}

	public Network growNetwork(int size, double joiningp){

		Network net = new Network();
		Network authorNet = new Network();
		BufferedWriter out = null;

		Random rand = new Random();
		//First Node
		Node firstNode  = new Node(0);
		firstNode.assignWeight(generateRandom(0.0,20.0)); //weight is equivalent to the impact factor
		net.addNode(firstNode);
		
		Node firstAuthor = new Node(0);//add a new author to authornet
		firstAuthor.assignWeight(100);
		authorNet.addNode(firstAuthor);

		//Add the rest of the nodes acccording to PA
		double p = joiningp; //Joining probability - a parameter
		for(int i = 1; i < size; i++)
		{
			Node node = new Node(i);
			node.assignWeight(generateRandom(0.0,20.0)); //weight is equivalent to the impact factor
			int[] authorsArray = generateAuthors(i,node.authorArray.length);//get the new authors and old authors count
			int[] randAraay = new int [(node.authorArray.length-authorsArray[node.authorArray.length])];
			for (int y=0;y<(node.authorArray.length-authorsArray[node.authorArray.length]);y++){//pick old authors randomly
				Integer auth= authorNet.getAllNodes().elementAt(generateRandom(1,authorNet.getAllNodes().size())).ID;
				while(Arrays.asList(node.authorArray).contains(auth)){
					auth= authorNet.getAllNodes().elementAt(generateRandom(1,authorNet.getAllNodes().size())).ID;
				}
				node.authorArray[y]=auth;
			}
			
			for(int z=0;z<authorsArray[node.authorArray.length];z++){//create nodes for new authors
				Node author = new Node(authorNet.networkSize+1);
				authorNet.addNode(author);
			}
			
			for (int j=0;j<node.authorArray.length;j++){
				for (int k=j+1;k<node.authorArray.length;k++){
					if(!(authorNet.getAllNodes().elementAt(node.authorArray[j]).isLinked(authorNet.getAllNodes().elementAt(node.authorArray[k])))){
						try 
						{
							net.addLink(authorNet.getAllNodes().elementAt(node.authorArray[j]), authorNet.getAllNodes().elementAt(node.authorArray[k]));
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
			}
			
			double probs[] =  new double[net.getSize()];
			int citationsToNetwork = (int)(citationToNetwork(i)*node.numberOfCitations);//citations to network as a function of time



			for(int j = 0; j <citationsToNetwork;j++)
			{
				Node destNode = net.getAllNodes().elementAt(j);
				double prob = p;
				if(net.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
				{
					//prob = p*((double)destNode.getNumberOfLinks()) / (double)net.getNoOfLinks();
					prob = (prob*destNode.getWeight()* (double)destNode.getNumberOfLinks()) / (double)net.getWeightedDegreeSum(); //For weighted PA
				}
				else //No links made yet
				{
					//prob = p;
					prob = prob*destNode.getWeight(); //For weighted PA
				}
				probs[j] = prob;
			}

			for(int j = 0; j <citationsToNetwork;j++)
			{
				Node destNode = net.getAllNodes().elementAt(j);
				//Throw dice
				if(rand.nextDouble() <=  probs[j])
				{
					try 
					{
						net.addLink(node, destNode);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}

			net.addNode(node);

			net.createTextFile("f://CNRes/paperNet_"+i+".txt");
			authorNet.createTextFile("f://CNRes/authorNet_"+i+".txt");

		}
		System.out.println("net size nodes "+net.getSize()+ " links "+ net.getNoOfLinks());

		return net;

	}

	public int generateRandom(int min, int max){
		return min + (int)(Math.random() * ((max - min) + 1)); //return a random number of citations to be initiated.
	}

	public double generateRandom(double min, double max){
		return min + (Math.random() * ((max - min) + 1)); //return a random number of citations to be initiated.
	}

	public double citationToNetwork(int i){

		return 0.7/500 * i;
	}

	public int[] generateAuthors(int i,int n){//n=number of authors for the paper, i=iteration variable to determine the percentage of old authors
		int [] authors = new int[n+1];
		int newA=0;
		int oldA=0;
		for(int j=0;j<n;j++){
			double r = generateRandom(0,1);
			if(r<(0.2*i/500)){
				oldA++;
				authors[j]=0;//pick an old author randomly
			}
			else
				newA++;
				authors[j]=1;//generate a new author
		}
		
		authors[n]=newA;
		return authors;
		
	}
	
	public int[] defineAuthors(int i, int n){
		int [] authors = generateAuthors(i,n);
		for(int j=0; j<n;j++){
			if (authors[j]==0){
				//pick an old author randomly
			}
			else if (authors[j]==1){
				//generate a new authors
			}
					
		}
		
		return authors;
	}




}
