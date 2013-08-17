package Basic;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.Vector;

import Distributions.EDistribution;
import Distributions.EDistributionGenerator;
import Distributions.QDistribution;
import Distributions.QDistribution_PowerLaw;
import oldClassesArchive.AssortativenessCalculator;

/**
 * This class is used to generate a network with given topology
 * Networks are hard coded
 * @author Piraveenan
 *
 */
public class NetworkGeneratorRefined {
	
	public Network generateStar(int size)
	{
		Network net = new Network();
		
		Node centrenode = new Node(size,(size-1) ); //Peripheral Nodes
		net.addNode(centrenode);
		
		for(int  i=1; i<= (size-1); i++)
		{
			Node node = new Node(i,1); //Peripheral Nodes
			
			net.addNode(node);
			try {
				net.addLink(node, centrenode);
			} catch (Exception e) {
				
				e.printStackTrace();
			}
		}
		
		
		
		return net;
	}
	
	public Network generateHubLinkStar()
	{
		//Size HardCoded
		Network net = new Network();
		
		Node centrenode = new Node(1,5 ); //Central Hub
		net.addNode(centrenode);
		
		for(int  i=1; i<= (5); i++)
		{
			Node pronode = new Node(i+1,6); //Provincial Hubs
			Node hubLinknode = new Node(i+6,2); //Hub Linkers
			Node proOutnode = new Node(i+11,5); //Outer Provincial Hubs;
		
			
			net.addNode(pronode);
			net.addNode(hubLinknode);
			net.addNode(proOutnode);
			try {
				net.addLink(pronode, centrenode);
				net.addLink(pronode, hubLinknode);
				net.addLink(proOutnode, hubLinknode);
			} catch (Exception e) {
				
				e.printStackTrace();
			}
			
			//Now add four peripheral Nodes for each Provincial or Outer Provincial Hub
			for(int  j=1; j<= (4); j++)
			{
				Node perinode = new Node(j+16+8*(i-1),1); //inner peripheral Nodes
				Node outerperinode = new Node(j+20+8*(i-1),1); //Outer Peripheral Nodes
				net.addNode(perinode);
				net.addNode(outerperinode);
				try {
					net.addLink(pronode, perinode);
					net.addLink(proOutnode, outerperinode);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
			}
			
		}
		
		net.createTextFile("testText");
		
		return net;
	}
	
	/** Generates a fully connected network of  N nodes
	 * 
	 * @param size
	 * @return
	 */
	public Network generateFullyConnected(int size)
	{
		//Size HardCoded
		Network net = new Network();
		
		for(int  i=1; i<= (size); i++)
		{
			Node node = new Node(i,size ); //Central Hub
			net.addNode(node);
		}
		
	
		Random rand = new Random();
		for(int  i=1; i<= size; i++)
		{
			int idx = i;
			
			for(int  j=(i+1); j<= size; j++)
			{
				int idy = j;
			
			
			
			try {
				net.addLink(idx, idy);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
			}
		}
		
		System.out.println("no of links "+net.getNoOfLinks());
		double N =  net.getSize();
		System.out.println("no of links should be "+(N-1)*N*0.5); //checking if the network is fully connected.
		
		return net;
	}

	/** Generates a fully connected  DIRECTED network of  N nodes
	 * All pairs of nodes are connected BOTH ways.
	 * @param size
	 * @return
	 */
	public Network generateFullyConnected_directed(int size)
	{
		//Size HardCoded
		Network net = new Network();
		
		for(int  i=1; i<= (size); i++)
		{
			Node node = new Node(i,size ); //Central Hub
			net.addNode(node);
		}
		
	
		Random rand = new Random();
		for(int  i=1; i<= size; i++)
		{
			int idx = i;
			
			for(int  j=(i+1); j<= size; j++)
			{
				int idy = j;
			
			
			
			try {
				net.addLink(idx, idy);
				net.addLink(idy, idx);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
			}
		}
		
		System.out.println("no of links "+net.getNoOfLinks());
		double N =  net.getSize();
		System.out.println("no of links should be "+(N-1)*N*0.5); //checking if the network is fully connected.
		
		return net;
	}

	public Network generateERRandom(int size)
	{
		//Size HardCoded
		Network net = new Network();
		
		for(int  i=1; i<= (size); i++)
		{
			Node node = new Node(i,size ); //Central Hub
			net.addNode(node);
		}
		
		double p = 0.5;
		int links = (int)Math.round(p*size*(size-1)*0.5);
		links = 3*size;
		
		Random rand = new Random();
		for(int  i=1; i<= (links); i++)
		{
			int idx = 1 + (Math.abs(rand.nextInt()) % size);
			int idy = 1 + (Math.abs(rand.nextInt()) % size);
			
			try {
				net.addLink(idx, idy);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
		}
		
		return net;
	}
	
	/**
	 * This method implements the interactive growth model proposed by Zhou
	 * @param network
	 * @param Mrand
	 * @return
	 */
	public Network growIGNet(int size)
	{
		//start with a random net of m0 nbodes and links
		//we take that m0 = 10
		
		int m0 = 20;
		
		Network net = generateERRandom(m0, m0);
		net  = generateStar(m0);
		Random rand  =  new Random(System.currentTimeMillis());	
		
		for(int i = (m0+1) ;i <= size ; i++)
		{
			
			
			if ( rand.nextDouble() <= 0.4) //with 40% probability
			{
				Node node = new Node(i); 
				net.addNode(node);
				
				//connnect it to one host node
				int linkIDMax = net.getNoOfLinks();
				int NodeIDMax = net.getSize()-2;
				long hostID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				//linear preference is implemented by randomizing link based
				//int hostID = 1 + (Math.abs((rand.nextInt()))%NodeIDMax); //implement linear preference
				try {
					net.addLink(i, (int)hostID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
				
				//coonect the host to twp peers
				long peer1ID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				long peer2ID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				
				
				try {
					net.addLink((int)peer1ID, (int)hostID);
					net.addLink((int)peer2ID, (int)hostID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
			}
				
			else //with 60 % probability
			{
				Node node = new Node(i); 
				net.addNode(node);
				
				//connnect it to two host node
				int NodeIDMax = net.getSize()-2;
				int linkIDMax = net.getNoOfLinks();
				long host1ID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				long host2ID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				
				
			
				try {
					net.addLink(i, (int)host1ID);
					net.addLink(i, (int)host2ID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
				
				
				
				//Connect one of the host nodes to peers
				//Here we choose the first hhost node, since both host nodes were randomly selected
				long peerID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				
				
				try {
					net.addLink((int)peerID,(int) host1ID);
			
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
				
			}
			
			//keep doing this until we get the required size
		}
		
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		return net;
	}
	
	/**
	 * This method implements the PARG model 
	 * @param net
	 * @param Npa
	 * @return
	 */
	public Network growPARGNet(int Npa)
	{
		Network net = generateERRandom(50);
		net.growPARGNet(Npa);
		
		return net;
	}

	
	/**
	 * This method implements the PFP growth model proposed by Zhou
	 * PRE 70-066108, 2004 and PRE 74-016124, 2006
	 * @param network
	 * @param Mrand
	 * @return
	 */
	public Network growPFPNet(int size, double delta)
	{
		//start with a random net of m0 nbodes and links
		//we take that m0 = 10
		
		int m0 = 20;
	
		
		
		
		Network net = generateERRandom(m0, m0);
		net  = generateStar(m0);
		Random rand  =  new Random(System.currentTimeMillis());	
		
		System.out.println("starts growth");
		int kcounter = 0;
		for(int i = (m0+1) ;i <= size ; i++)
		{
			
			if(i%1000==0)
			{
				kcounter++;
				System.out.println(kcounter +" thousand Nodes grown");
			}
			
			if ( rand.nextDouble() <= 0.4) //with 40% probability
			{
				Node node = new Node(i); 
				net.addNode(node);
				
				//connnect it to one host node
				int linkIDMax = net.getNoOfLinks();
				int NodeIDMax = net.getSize()-2;
				
				double pfpSUM = 0;
				
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
				}
				
				double randPFP = rand.nextDouble()*pfpSUM;
				//Now convert this to a node ID
				long hostID = 0;
				
				pfpSUM = 0;
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
					if(pfpSUM >= randPFP )
					{
						//This is the id we should choose
						hostID = net.getAllNodes().elementAt(ind).getID();
						break;
					}
				}
				
				
				
				
				
				//coonect the host to two peers
				long peer1ID =0;
				long peer2ID=0;
				
				
				pfpSUM = 0;
				
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
				}
				
				 randPFP = rand.nextDouble()*pfpSUM;
				//Now convert this to a node ID
				peer1ID = 0;
				
				pfpSUM = 0;
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
					if(pfpSUM >= randPFP )
					{
						//This is the id we should choose
						peer1ID = net.getAllNodes().elementAt(ind).getID();
						break;
					}
				}
				

				pfpSUM = 0;
				
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
				}
				
				 randPFP = rand.nextDouble()*pfpSUM;
				//Now convert this to a node ID
				peer2ID = 0;
				
				pfpSUM = 0;
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
					if(pfpSUM >=randPFP )
					{
						//This is the id we should choose
						peer2ID = net.getAllNodes().elementAt(ind).getID();
						break;
					}
				}
				
			
				try {
					net.addLink(i, (int)hostID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
				
				
				try {
					net.addLink((int)peer1ID, (int)hostID);
					net.addLink((int)peer2ID, (int)hostID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
			}
				
			else //with 60 % probability
			{
				Node node = new Node(i); 
				net.addNode(node);
				
				//connnect it to two host node
				int NodeIDMax = net.getSize()-2;
				int linkIDMax = net.getNoOfLinks();
				long host1ID = 0;
				long host2ID = 0;
				

				double pfpSUM = 0;
				
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
				}
				
				double  randPFP = rand.nextDouble()*pfpSUM;
				//Now convert this to a node ID
				host1ID = 0;
				
				pfpSUM = 0;
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
					if(pfpSUM >= randPFP )
					{
						//This is the id we should choose
						host1ID = net.getAllNodes().elementAt(ind).getID();
						break;
					}
				}
				

				pfpSUM = 0;
				
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
				}
				
				 randPFP = rand.nextDouble()*pfpSUM;
				//Now convert this to a node ID
			
				
				pfpSUM = 0;
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
					if(pfpSUM >= randPFP )
					{
						//This is the id we should choose
						host2ID = net.getAllNodes().elementAt(ind).getID();
						break;
					}
				}
				
				
			
				
				
				
				
				//Connect one of the host nodes to peers
				//Here we choose the first hhost node, since both host nodes were randomly selected
				long peerID = 0; 
				

				pfpSUM = 0;
				
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
				}
				
				 randPFP = rand.nextDouble()*pfpSUM;
				//Now convert this to a node ID
				peerID = 0;
				
				pfpSUM = 0;
				for(int ind = 0; ind < net.getSize(); ind++)
				{
					double pfpDegree =   net.getAllNodes().elementAt(ind).getNumberOfLinks();
					if(net.getAllNodes().elementAt(ind).getID() == node.getID())
					{
						continue; //This is the node we are adding NOW. We dont want to consider this!!
					}
					if(pfpDegree!= 0)
					{
						pfpSUM =  pfpSUM + Math.pow(pfpDegree, 1.0 + delta*Math.log(pfpDegree));
					}
					if(pfpSUM >= randPFP )
					{
						//This is the id we should choose
						peerID = net.getAllNodes().elementAt(ind).getID();
						break;
					}
				}
				
				try {
					net.addLink(i, (int)host1ID);
					net.addLink(i, (int)host2ID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
				
				
				try {
					net.addLink((int)peerID,(int) host1ID);
			
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
				
			}
			
			//keep doing this until we get the required size
		}
		
		System.out.println("Number of Nodes now  = " + net.getSize());
		//System.out.println("Number of links now  = " + net.getAllLinks().size());
	
		return net;
	}
	
	
	
	public Network generateERRandom(int size, int links)
	{
		//Size HardCoded
		Network net = new Network();
		
		for(int  i=0; i<= (size); i++)
		{
			Node node = new Node(i,size ); //Central Hub
			net.addNode(node);
		}
		
	
		Random rand = new Random();
		for(int  i=0; i<= (links); i++)
		{
			int idx = 1 + (Math.abs(rand.nextInt()) % size);
			int idy = 1 + (Math.abs(rand.nextInt()) % size);
			
			try {
				net.addLink(idx, idy);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
		}
		
		return net;
	}
	
	/**
	 * This method generates a small-world network as described by Watts & Strogartz: 'Collective Dynamics of Small world networks' -  Nature 1998.
	 * @param N Number of nodes
	 * @param K Number of nearest neighbours connected in the orginical ring -  Needs to be EVEN NUMBER
	 * @param p rewiring probability
	 * @return small world network
	 * Note well: Number of Links M = 0.5*N*K
	 */
	public Network generateSmallWorld(int N, int K, double p)
	{
		Network net = new Network();
		Random rand = new Random();
		int halfK = (int)(0.5*K);
		int M = halfK*N;
		int[][] links =  new int[N][N];
		
		// add nodes
		for(int  i=0; i< N; i++)
		{
			Node node = new Node(i,N); 
			net.addNode(node);
		}
		
		//creat ring
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
						net.addLink(x,y);
					} catch (Exception e) {
						System.out.println("SMALL WORLD NET: COULD NOT CREATE LINK BETWEEN "+x+ " and "+y);
						e.printStackTrace();
					}
				}
			}
		}
		
		System.out.println("SMALL WORLD NET CREATED WITH "+net.getSize()+" NODES AND "+net.getNoOfLinks()+ " LINKS");
		return net;
	}


	/**
	 * Generate a network from given datafile
	 * @param size
	 * @return
	 */

	public Network generateFromFile()
	{
		
		//String filename = ".\\..\\Data\\node rho\\PA 5000 0-75 for PLOSONE REVIEW.txt";
		//String filename = ".\\..\\Data\\kenneth\\net2001.txt"; as20040105
		
		String filename = ".\\..\\Data\\hi\\smallworld_1000_4_2.txt";
		//String filename = ".\\..\\Data\\hi\\sample SW 0-5.txt";
		Network net = new Network();
		FileProcessor  x= new FileProcessor();
	    x.SetFile(filename);
	    x.Process(); 
	    x.showall();
	    x.buildlinks();
	    net = x.getNetwork();
	    
		net.createTextFile("Latest Net Created");

		return net;
	}
	
	public Network generateFromMultipleFiles()
	{
		
		String filename = ".\\..\\Data\\node rho\\test";
		String fname="";
	
		Network net = new Network();
		FileProcessor  x= new FileProcessor();
		
		for(int i=0; i< 3; i++)
		{
	
			fname =  filename+""+i+".txt";
		    x.SetFile(fname);
		    x.Process(); 
		    x.showall();
		    x.buildlinks();
		    net = x.getNetwork();
		}
	    
		net.createTextFile("Latest Net Created");
		net.printToConsole();

		return net;
	}
	
	public Network generateDynamicNetFromFile()
	{
		
		String filename = ".\\..\\Data\\schoolnet\\Schoolnet I above 50.txt";
		//String filename = ".\\Data\\Guimera_internet\\ecoli-properIDs.txt";
		//String filename = ".\\Data\\PA\\PADENSE3.txt";
	
		Network net = new Network();
		FileProcessor  x= new FileProcessor();
	    x.SetFile(filename);
	    x.Process_ForTimeNodes(); 
	    x.showall();
	    x.buildlinks();
	    net = x.getNetwork();
	    
		net.createTextFile("Lstest Net Created");

		return net;
	}
	
	
	public Network generateFromFile(String NAME)
	{

		String filename =NAME+".txt";
		
		Network net = new Network();
		FileProcessor  x= new FileProcessor();
	    x.SetFile(filename);
	    x.Process(); 
	    x.showall();
	    x.buildlinks();
	    net = x.getNetwork();
	    
		net.createTextFile("LAtest Net Created");

		return net;

		
	}
	
	/**
	 * Grow a network of a given size using preferential attachment
	 * @param size
	 * @return
	 */
	public Network growNetworkPA(int size, double joiningp)
	{
		BufferedWriter out = null;

		Random rand = new Random();
		
		Network net = new Network();
		//First Node
		Node firstNode  = new Node(0);
		firstNode.assignWeight(rand.nextDouble());
		net.addNode(firstNode);
		
		//Add the rest of the nodes acccording to PA
		double p = joiningp; //Joining probability - a parameter
		for(int i = 1; i < size; i++)
		{
			Node node = new Node(i);
			node.assignWeight(rand.nextDouble());
			double probs[] =  new double[net.getSize()];
			
			for(int j = 0; j <net.getSize();j++)
			{
				Node destNode = net.getAllNodes().elementAt(j);
				double prob = p;
				if(net.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
				{
					 prob = p*((double)destNode.getNumberOfLinks()) / (double)net.getNoOfLinks();
					// prob = (prob*destNode.getWeight()* (double)net.getNoOfLinks()) / (double)net.getWeightedDegreeSum(); //For weighted PA
				}
				else //No links made yet
				{
					prob = p;
					//prob = prob*destNode.getWeight(); //For weighted PA
				}
				probs[j] = prob;
			}
			
			for(int j = 0; j <net.getSize();j++)
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
			net.createTextFile("f://CNRes/PA_"+i+".txt");
			
		}
		System.out.println("net size nodes "+net.getSize()+ " links "+ net.getNoOfLinks());
		return net;
	}
	
	public Network growNetworkAPA(int size, int Np, double gamma, double r)
	{
		
		EDistributionGenerator e = new EDistributionGenerator();
		double calculated_r = 0;
		double calculated_i = 0;
		
		int rTotalSteps = 20;
		
		QDistribution_PowerLaw qk = new QDistribution_PowerLaw(Np , gamma, true);
		//QDistribution_Poisson qk = new QDistribution_Poisson(Np , (Np/2.0));
		if(!checkQDistValidity(qk)){
			System.out.println("INVALID Qk");
			qk.print();
			System.out.println("PROGRAM EXITS");
			System.exit(0);
		}
		
		/* Generate E dist for the given  q, r */
		EDistribution ejk = e.generateEDistribution_class_A(qk, r, rTotalSteps);
		//EDistribution ejk = e.generateEDistribution_class_B(qk, r, rTotalSteps);
		//System.out.println("NEW TEMPLATE Ejk");
		//ejk.print();
		//qk.print();
		//System.out.println("avg k"+ qk.getAverageDegree());
		double[] q = qk.getDist();
		for (int i = 0; i < qk.getLength(); i++) 
		{
			for (int j = 0; j < qk.getLength(); j++) 
			{
				if(i==j)
				{
					//System.out.print(" "+q[i]*q[j]);
				}
				else if(i==(q.length - j-1))
				{
					//System.out.print(" "+q[i]*q[q.length - j-1]);
				}
				else
				{
					
				}
				
			}
			//System.out.println("");
		}
		boolean b = checkEDistValidity(ejk, qk);
		if(!b)
		{
			System.out.println("INVALID Ejk");
			ejk.print();
			qk.print();
			System.out.println("PROGRAM EXITS");
			System.exit(0);
		}
		
		//Now we have ejk, and qk
		
		return buildAPANetwork(qk, ejk, size);
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
		
		/**
		 * To force single giant component, connect together existing nodes in order
		 * 
		 */
		
//		for(int y=0; y < network.getSize(); y++)
//		{
//			Node node1 =  network.getAllNodes().elementAt(y);
//			Node node2 =  network.getAllNodes().elementAt((y+1)%network.getSize() );
//			try {
//				network.addLink(node1, node2);
//			} catch (Exception e) {
//				System.out.println("trouble adding giant component guarantee links");
//				e.printStackTrace();
//			}
//			
//		}

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
			// Node node = new Node(m0+i, qDist.getRandomDegreeDUMMY(i));
			// /*TODO: Hck random change back */
			Node node = joiningNodes.elementAt(i); // new Node(m0+i,
													// qDist.getRandomDegree());
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

			/* Measure r */
			/*
			 * if( (eDist.Expectation()<4.48) && (eDist.Expectation()>4.38) ){
			 * AssortativenessCalculator ac = new
			 * AssortativenessCalculator(network); double ra = ac.calculate_r();
			 * System.out.println("r now: "+ra); }
			 */

			/* update display */
			// update_NETEDistDataset(eDist, network.getLinkDistribution());
			// update_QDistDataset(qDist, network.getNodeDistribution());

		}
		// network.printArray(network.getLinkDistribution(),"link dist");
		return network;

	}
	
	/**
	 * 'Scale' a Edistribution
	 */
	private double[] getScaledEDist(int joiningNodeQ, EDistribution eDist,
			QDistribution_PowerLaw qDist) {
		if (joiningNodeQ < 0) {
			System.out
					.println("A Node Tries to join with Unsuitable remaining Degree: Error");
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


	
	public static boolean checkQDistValidity(QDistribution q) {
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
					if(Math.abs(eArr[j][k]) < Math.pow(10, -10))
					{
						eArr[j][k] = 0; //Too small to bother
						System.out.println(" < 0,  but too small to bother");
					}
					else
					{

						System.out.println(" < 0");
						return false;
					}

				}

				if (eArr[j][k] > 1) {
					System.out.println(" > 1");
					return false;
				}

				if (eArr[k][j] != eArr[j][k]) {
					if (Math.abs(eArr[k][j] - eArr[j][k]) < Math.pow(10, -10) )
					{
						System.out.println("non symmetry - too small to bother");
						eArr[k][j] = eArr[j][k];
					}
					else
					{
						System.out.println("non symmetry");
						return false;
					}
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
	
	public Network growNetwork(int size, double joiningp) throws Exception{

		Network net = new Network();
		Network authorNet = new Network();
		BufferedWriter out = null;

		Random rand = new Random();
		//First Node
		Node firstNode  = new Node(0);
		firstNode.assignWeight(generateRandom(0.0,20.0)); //weight is equivalent to the impact factor
		net.addNode(firstNode);
		firstNode.authorArray = new Integer[3];
		
		Node firstAuthor = new Node(0);//add a new author to authornet
		authorNet.addNode(firstAuthor);
		Node secondAuthor = new Node(1);
		authorNet.addNode(secondAuthor);
		authorNet.addLink(firstAuthor,secondAuthor);
		Node thirdAuthor = new Node(2);
		authorNet.addNode(thirdAuthor);
		authorNet.addLink(secondAuthor,thirdAuthor);
		authorNet.addLink(firstAuthor, thirdAuthor);
		Integer [] au = {0,1,2};
		firstNode.authorArray=au;
		firstAuthor.papers.add(firstNode.getID());
		secondAuthor.papers.add(firstNode.getID());
		thirdAuthor.papers.add(firstNode.getID());
		int newAuthors=0;

		//Add the rest of the nodes acccording to PA
		double p = joiningp; //Joining probability - a parameter
		for(int i = 1; i < size; i++)
		{
			Node node = new Node(i);
			node.assignWeight(generateRandom(0.0,20.0)); //weight is equivalent to the impact factor
			int[] authorsArray = generateAuthors(i,node.authorArray.length);//get the new authors and old authors count
			newAuthors=authorsArray[node.authorArray.length];
			//Integer[] randArray = new Integer [node.authorArray.length];
			for (int y=0;y<(node.authorArray.length-authorsArray[node.authorArray.length]);y++){//pick old authors randomly
				int gen = generateRandom(0,(authorNet.getSize()-1));
				//System.out.println(gen);
				//System.out.println(authorNet.getNode(gen).getID());
				Integer auth= authorNet.getAllNodes().elementAt(gen).ID;
				//System.out.println("auth is: "+auth + " gen is: "+gen);
				while(Arrays.asList(node.authorArray).contains(auth)){
					auth= authorNet.getNode(generateRandom(1,authorNet.getSize()-1)).getID();
				}
				node.authorArray[y]=auth;
				authorNet.getNode(auth).papers.add(node.getID());
			}
			
			for(int z=(node.authorArray.length-authorsArray[node.authorArray.length]);z<node.authorArray.length;z++){//create nodes for new authors
				node.authorArray[z] = authorNet.getSize();
				Node author = new Node(authorNet.getSize());
				authorNet.addNode(author);
				author.papers.add(node.getID());
			}
			
			
			//node.authorArray=randArray;
			
			System.out.println("Iteration: "+i+" AuthorNet Size: "+authorNet.getSize()+" New Authors: " + newAuthors);
			
			for (int j=0;j<node.authorArray.length;j++){
				for (int k=j+1;k<node.authorArray.length;k++){
					//System.out.println(j+" and "+k);
					//System.out.println(authorNet.getNode(j).getID());
					//System.out.println(node.authorArray[1]);
					if(!(authorNet.getNode(node.authorArray[j]).isLinked(authorNet.getNode(node.authorArray[k])))){
						try 
						{
							//System.out.println("WOOT");
							authorNet.addLink(authorNet.getAllNodes().elementAt(node.authorArray[j]), authorNet.getAllNodes().elementAt(node.authorArray[k]));
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
			}
			
			double probs[] =  new double[net.getSize()];
			int citationsToNetwork = (int)(citationToNetwork(i,node.numberOfCitations)*node.numberOfCitations);//citations to network as a function of time

			for(int j = 0; j <net.getSize();j++)
			{
				Node destNode = net.getAllNodes().elementAt(j);
				double prob = p;
				if(net.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
				{
					//prob = p*((double)destNode.getNumberOfLinks()) / (double)net.getNoOfLinks();
					prob = (prob*destNode.getWeight()* (double)destNode.getNumberOfInLinks()) / (double)net.getWeightedDegreeSum(); //For weighted PA
				}
				else //No links made yet
				{
					//prob = p;
					prob = prob*destNode.getWeight(); //For weighted PA
				}
				probs[j] = prob;
			}
			int links=0;
			for(int j = 0; j <probs.length && links<citationsToNetwork;j++)
			{
				Node destNode = net.getAllNodes().elementAt(j);
				//Throw dice
				if(rand.nextDouble() <=  probs[j])
				{
					try 
					{
						net.addLink(node, destNode);
						System.out.println(node.getID()+" cited "+destNode.getID());
						links++;
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				//System.out.println("CitationsToNetwork = "+j);
			}

			net.addNode(node);
			System.out.println("Number of Authors:" + node.authorArray.length);
			System.out.println("Number of Citations:" + node.numberOfCitations + " Number of Citations to the Network:" + citationsToNetwork);

			net.createTextFile("f://CNRes/paperNet_"+i+".txt");
			authorNet.createTextFile("f://CNRes/authorNet_"+i+".txt");
			
			//calculate page rank
			ArrayList <Integer> hIndex=calculatehIndex(authorNet,net);//calculate h-index
			System.out.println(hIndex.toString());
			for (int f=0;f<authorNet.getSize();f++){
				ArrayList <Integer> papers = authorNet.getNode(f).papers;
				//System.out.println(papers.toString());
			}
			
		}
		System.out.println("net size nodes "+net.getSize()+ " links "+ net.getNoOfLinks());
		System.out.println("net size nodes "+authorNet.getSize()+ " links "+ authorNet.getNoOfLinks());
		
		
		

		return net;

	}
	
	public ArrayList<Integer> calculatehIndex(Network authors, Network papernet) throws Exception{
		
		ArrayList <Integer> hIndex = new ArrayList<Integer>();
		
		for (int i=0;i<authors.getSize();i++){
			ArrayList <Integer> papers = authors.getNode(i).papers;
			ArrayList <Integer> nCite = new ArrayList <Integer>();
			for (int j=0;j<papers.size();j++){
				Node paperNode=papernet.getNode(papers.get(j));
				nCite.add(paperNode.getNumberOfInLinks());
			}
			Collections.sort(nCite);
			Collections.reverse(nCite);
			System.out.println("Author: "+i+" "+nCite.toString());
			//int hInd=0;
			int k=0;
			for (k=0;k<nCite.size();k++){
				if (k==0 && nCite.get(k)==0) break;
				if (k >=nCite.get(k)) break;
			}
			authors.getNode(i).hIndex=k;
			hIndex.add(k);
		}
		return hIndex;
	}

	public int generateRandom(int min, int max){
		return min + (int)(Math.random() * ((max - min) + 1)); //return a random number of citations to be initiated.
	}
	
	public int generateRandome(int min, int max){
		return min + (int)(Math.random() * ((max - min))); //return a random number of citations to be initiated exclusive of max.
	}

	public double generateRandom(double min, double max){
		return min + (Math.random() * ((max - min) + 1)); //return a random number of citations to be initiated.
	}

	public double citationToNetwork(int i, int j){
		if (i<500){
			return 0.7/500 * Math.max(i, j);
		}
		else
			return 0.7;
	}

	public int[] generateAuthors(int i,int n){//n=number of authors for the paper, i=iteration variable to determine the percentage of old authors
		int [] authors = new int[n+1];
		int newA=0;
		int oldA=0;
		for(int j=0;j<n;j++){
			double r = generateRandom(0,1);
			if(i<500){
				if(r<(0.2*i/500)){
					oldA++;
					authors[j]=0;//pick an old author randomly
				}
				else
					newA++;
					authors[j]=1;//generate a new author
			}
			else if (i>=500){
				if(r<0.2){
					oldA++;
					authors[j]=0;//pick an old author randomly
				}
				else
					newA++;
					authors[j]=1;//generate a new author
			}
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
