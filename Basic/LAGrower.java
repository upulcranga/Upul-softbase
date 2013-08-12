package Basic;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;


import Distributions.QDistribution_PowerLaw;

public class LAGrower {
	
	public static void main(String[] args) {
		

		LAGrower sim = new LAGrower();
		sim.simulate();

	}
	
	public void simulate()
	{
		//Set filenames
		String infilename = ".\\Data\\simplenet\\lamotif10.txt";
		//String infilename = ".\\Data\\Guimera_internet\\cocomac whole.txt";
		//String infilename = ".\\Data\\internetCaida\\as20040105.txt";
		String outfilenameONE = ".\\LA-Extended\\GrowthSimulation\\GROWNMODEL-1.xls";
		String outfilenameTWO = ".\\LA-Extended\\GrowthSimulation\\GROWNMODEL-2";
		  BufferedWriter out = null;
		try {
			out = new BufferedWriter(new FileWriter(outfilenameONE));	
		}	
		 catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace(); 
		}
		 
		 
		//Read a network from file
		Network net = new Network();
		FileProcessor  x= new FileProcessor();
	    x.SetFile(infilename);
	    x.Process(); 
	    x.showall();
	    x.buildlinks();
	    net = x.getNetwork();
	    
	      NetworkGenerator gen = new NetworkGenerator();
	     net  = gen.generateERRandom(50);
		
	  System.out.println("INPUT FILE READ");
	  System.out.println("Number of Nodes  = " + net.getSize());
	System.out.println("Number of links now  = " + net.getAllLinks().size());
		

		 
		//calculate and print local assortativeness distribution
	    double[] rl_dist = get_rl_dist(net);
	    boolean ZERO_OFF = false;
				try {	
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
					  System.out.println("LA dist printed");
					
				} catch (IOException e) {
					System.out.println("cant write");
					e.printStackTrace();
				}
				
				
		
		//Grow in Phase two
				int extraNodes = 3500;
				
				
				int extraLinks = 5000;
				
				
				for(int p=0; p< 1; p++)
				{
					//extraLinks = (40000);
					//net = growPhaseTwo_PA(net, extraNodes);
					//net = growIGNet(net, extraNodes);
					net = growASPA(net, extraNodes);
					net.createTextFile(".\\LA-Extended\\GrowthSimulation\\GROWNMODELFILE");
				
					
					System.out.println("Phase Two grown");
				}
				
		
		//claculate and print LA distribution again
				try {
					out = new BufferedWriter(new FileWriter(outfilenameTWO+".xls"));	
				}	
				 catch (IOException e) {
					System.out.println("cant write to file");
					e.printStackTrace();
				}
			rl_dist = get_rl_dist(net);
							try {	
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
								 System.out.println("LA dist printed after phase 2 growth");
								
							} catch (IOException e) {
								System.out.println("cant write");
								e.printStackTrace();
							}
				//}
	}
	
	public double[] get_rl_dist(Network network)
	{
		int Np = network.getDegreeDistribution().length;
		double[] rl = new double[Np];
		double times[] = new double[Np];
		for(int y=0; y<network.getAllNodes().size();y++)
		{
			Node thisNode = network.getAllNodes().elementAt(y);
			times[thisNode.getNumberOfLinks()] = times[thisNode.getNumberOfLinks()]+1;
			rl[thisNode.getNumberOfLinks()] = rl[thisNode.getNumberOfLinks()]+calculate_localR(network, thisNode); //Accumulated r for that degree
			
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
	
	public double calculate_localR(Network network, Node node)
	{
		if(node.getNumberOfLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
		double N = network.getSize();
		double Nlinks = network.getNoOfLinks();
		QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(network);
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
	
	/**
	 * THis method adds N number of links randomly
	 * @param network
	 * @return
	 */
	
	public Network growPhaseTwo(Network network, int Mrand)
	{
		System.out.println("Number of links now  = " + network.getAllLinks().size());
		Random rand  =  new Random(System.currentTimeMillis());
		int NodeIDMax = network.getSize();
		
		for(int i=0; i< Mrand; i++)
		{
			int sourceID = Math.abs((rand.nextInt()))%NodeIDMax;
			int destID = Math.abs((rand.nextInt()))%NodeIDMax;
			try {
				network.addLink(sourceID, destID);
			} catch (Exception e) {
				System.out.println("problem adding random link");
				e.printStackTrace();
			}
		}
		System.out.println("Number of links now  = " + network.getAllLinks().size());
		return network;
	}
	
	/**
	 * This method implements the interactive growth model proposed by Zhou
	 * @param network
	 * @param Mrand
	 * @return
	 */
	public Network growIGNet( Network net,int size)
	{
		int initialSize = net.getSize();
	
		Random rand  =  new Random(System.currentTimeMillis());	
		
		for(int i = 1 ;i <= size ; i++)
		{

			if ( rand.nextDouble() <= 0.4) //with 40% probability
			{
				Node node = new Node(initialSize + i); 
				net.addNode(node);
				
				//connnect it to one host node
				int linkIDMax = net.getNoOfLinks();
				int NodeIDMax = net.getSize()-2;
				long hostID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				//linear preference is implemented by randomizing link based
				//int hostID = 1 + (Math.abs((rand.nextInt()))%NodeIDMax); //implement linear preference
				try {
					net.addLink(initialSize + i, (int)hostID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
				
				//coonect the host to two peers
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
				Node node = new Node(initialSize + i); 
				net.addNode(node);
				
				//connnect it to two host node
				int NodeIDMax = net.getSize()-2;
				int linkIDMax = net.getNoOfLinks();
				long host1ID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				long host2ID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				
				
			
				try {
					net.addLink(initialSize + i, (int)host1ID);
					net.addLink(initialSize + i, (int)host2ID);
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
	 * This method adds a given number of nodes in the Preferential Attachment way.
	 * i.e the chance of a node getting a link is proportional to the links it already has
	 * NOTE: It seems that PA is not a good enough model for internet growth!!
	 * @param network
	 * @param Npa the number of extra ndoes to be added
	 * @return
	 */
	
	public Network growPhaseTwo_PA(Network net, int Npa)
	{
		int initialSize = net.getSize();
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		Random rand  =  new Random(System.currentTimeMillis());
//		
		double p = 0.3; //Joining probability - a parameter
		for(int i = 1; i <= Npa; i++)
		{
			Node node = new Node(initialSize + i);
			//node.assignWeight(rand.nextDouble());
			double probs[] =  new double[net.getSize()];
			
			for(int j = 0; j <net.getSize();j++)
			{
				Node destNode = net.getAllNodes().elementAt(j);
				double prob = p;
				if(net.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
				{
					 prob = p*((double)destNode.getNumberOfLinks()) / (double)net.getNoOfLinks();
					 //prob = (prob*destNode.getWeight()* (double)net.getNoOfLinks()) / (double)net.getWeightedDegreeSum(); //For weighted PA
				}
				else //No links made yet
				{
					prob = p*(0.00/(double)net.getNoOfLinks());
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
		}
		
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		return net;
	}
	
	
	
	/**
	 * This method adds a given number of *LINKS* in the Preferential Attachment way.
	 * i.e the chance of a node getting a link is proportional to the links it already has
	 * @param network
	 * @param Mpa the number of extra ndoes to be added
	 * @return
	 */
	
	public Network growPhaseTwo_PA2(Network net, int Mpa)
	{
	
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		Random rand  =  new Random(System.currentTimeMillis());	
	
		for(int i = 1; i <= Mpa; i++)
		{
			int linkIDMax = net.getNoOfLinks();
			//Choose a link ID randomly, and choose the corrersponding node
			//Therefore the probability of a node being chosen is proportional to the links it already has.
			long sourceID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
			long destID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).jNodeID(); 
			try {
				net.addLink((int)sourceID, (int)destID);
			} catch (Exception e) {
				System.out.println("problem adding PA link");
				e.printStackTrace();
			}
			
		}
		
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		return net;
	}
	
	/**
	 * This method is the implementation of the new algorithm for AS network growth
	 * It has two mechanisms of growth. One is a Preferential attachemnt based mechanism,
	 * the other is an 'adjustment' mechanism where hub-to-hub connections are lost.
	 * @param network
	 * @param Npa the number of extra ndoes to be added
	 * @return
	 */
	
	public Network growASIG(Network net, int Npa)
	{
		int initialSize = net.getSize();
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		Random rand  =  new Random(System.currentTimeMillis());
//		
		double EjoinSTART = 0.8;
		double EdeleteSTART =3.02;
		double Ejoin = EjoinSTART; //Expectation of links created per node- a parameter
		double Edelete = EdeleteSTART; // Expected links deleted per node addition. A parameter
								
		for(int i = 1; i <= Npa; i++)
		{
			Ejoin = EjoinSTART; // + (0.02*(double)(net.getSize())); //increase addded links with growing network size
			Edelete = EdeleteSTART; // + (0.002*(double)(net.getSize())); //increase addded links with growing network size
			Node node = new Node(initialSize + i);
			double probs[] =  new double[net.getSize()];
			int noOfLinks = net.getNoOfLinks();
			
			if ( rand.nextDouble() <= 0.4) //with 40% probability
			{
			
				net.addNode(node);
				
				//connnect it to one host node
				int linkIDMax = net.getNoOfLinks();
				int NodeIDMax = net.getSize()-2;
				long hostID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				//linear preference is implemented by randomizing link based
				//int hostID = 1 + (Math.abs((rand.nextInt()))%NodeIDMax); //implement linear preference
				try {
					net.addLink(initialSize + i, (int)hostID);
				} catch (Exception e) {
					System.out.println("problem adding random link");
					e.printStackTrace();
				}
				
				//connect the host to two peers: again, peers with bigger degrees are more likely to be chosen.
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
			
				net.addNode(node);
				
				//connnect it to two host nodes
				int NodeIDMax = net.getSize()-2;
				int linkIDMax = net.getNoOfLinks();
				long host1ID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				long host2ID = (net.getAllLinks().elementAt(Math.abs((rand.nextInt()))%linkIDMax)).iNodeID(); 
				
				
			
				try {
					net.addLink(initialSize + i, (int)host1ID);
					net.addLink(initialSize + i, (int)host2ID);
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
			
	
	
		
			
			 noOfLinks = net.getNoOfLinks();
				System.out.println("Number of links now  = " + net.getAllLinks().size());
				int linksNow = net.getAllLinks().size();
			
			//Now, stochastically delete some links that connect hubs
			
			double raw = 0.99; //A parameteter, determining whether stochastic deleted needs to be done
			//If this is set to 1.0, stochastic deleting is ALWAYS done.
			
			
			
			net.sortLinksBasedOnCombinedDegree();
			
			if( rand.nextDouble() <=  raw ) 
			{
				//destroy 'Edeleted' number of the existing links
				double probs2[] =  new double[net.getNoOfLinks()];
				for(int k = 0; k <net.getNoOfLinks();k++)
				{
					Link destLink = net.getAllLinks().elementAt(k);
					double prob = ((double)Edelete)*0.25*0.01*(double)net.getNoOfLinks()*((double)((destLink.iNodeLinks() )+( destLink.jNodeLinks())) / (double)net.getNoOfLinks());
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
					if(k >= (net.getNoOfLinks() -  Edelete)) //This is a BIG (hub-hub) link, should be deleted.
						//ASSORTATIVE conenctions should be discouraged here; BOTH hub-to-hub as well as periphery-to-periphery
					{
						prob = 1.0;
					}
					else
					{
						prob = 0;
					}
					//prob = 0;
					 probs2[k] = prob;
					 
					 //randomly select some links from the 'assortative links' and replace them
				}
				
				
				
				for(int k = 0; k <net.getNoOfLinks();k++)
				{
					Link destLink = net.getAllLinks().elementAt(k);
					//Throw dice
					if(rand.nextDouble() <=  probs2[k])
					{
						try 
						{
							//Delete Link
							
							//Deleting immediately from network will cause problems with size, therefore mark and delete later.
							net.getAllLinks().elementAt(k).markForDeletion();
							Node iNode = net.getAllLinks().elementAt(k).iNode();
							Node jNode = net.getAllLinks().elementAt(k).jNode();
							System.out.println("Link REPLACED "+ iNode.getNumberOfLinks()+ " "+jNode.getNumberOfLinks());
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
				net.replaceMarkedLinks(); //replace deleted links by connecting the corresponding nodes to other nodes selected preferentially,
				//so that a deleted Link does not result in breaking the network into components
				 noOfLinks = net.getNoOfLinks();
				 int la = 0;
				
				//Now the links have been deleted and second mechanism completed
					//System.out.println("Number of Nodes now  = " + net.getSize());
					//System.out.println("Number of links now  = " + net.getAllLinks().size());
				//System.out.println("DELETED Links in this node addition  = " + (linksNow - net.getAllLinks().size()));
				
			}
			
			
		}
		
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		return net;
	}
	
	/**
	 * This method is the implementation of the new algorithm for AS network growth with the IG growth incorporated.
	 * It has two mechanisms of growth. The first one is the
	 * IG growth mechanism suggested by Zhou (rather than the Preferential attachemnt based mechanism in the method above),
	 * the other is an 'adjustment' mechanism where hub-to-hub connections are lost.
	 * @param network
	 * @param Npa the number of extra ndoes to be added
	 * @return
	 */
	
	public Network growASNet(Network net, int Npa)
	{
		int initialSize = net.getSize();
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		Random rand  =  new Random(System.currentTimeMillis());
//		
		double EjoinSTART = 0.9;
		double EdeleteSTART =1.02;
		double Ejoin = EjoinSTART; //Expectation of links created per node- a parameter
		double Edelete = EdeleteSTART; // Expected links deleted per node addition. A parameter
								
		for(int i = 1; i <= Npa; i++)
		{
			Ejoin = EjoinSTART; // + (0.02*(double)(net.getSize())); //increase addded links with growing network size
			Edelete = EdeleteSTART; // + (0.02*(double)(net.getSize())); //increase addded links with growing network size
			Node node = new Node(initialSize + i);
			double probs[] =  new double[net.getSize()];
			int noOfLinks = net.getNoOfLinks();
			
			for(int j = 0; j <net.getSize();j++)
			{
				Node destNode = net.getAllNodes().elementAt(j);
				double prob = Ejoin;
				if(net.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
				{
					 prob = 0.5*Ejoin*((double)destNode.getNumberOfLinks()) / (double)net.getNoOfLinks();
				}
				else //No links made yet
				{
					prob = 1.0;//If no links are made, this node node has to join with the existing isolated nodes
					//isolated nodes are not allowed to exist.
					
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
			
			 noOfLinks = net.getNoOfLinks();
				System.out.println("Number of links now  = " + net.getAllLinks().size());
				int linksNow = net.getAllLinks().size();
			
			//Now, stochastically delete some links that connect hubs
			
			double raw = 0.05; //A parameteter, determining whether stochastic deleted needs to be done
			//If this is set to 1.0, stochastic deleting is ALWAYS done.
			
			net.sortLinksBasedOnCombinedDegree();
			
			if( rand.nextDouble() <=  raw ) 
			{
				//destroy 'Edeleted' number of the existing links
				double probs2[] =  new double[net.getNoOfLinks()];
				for(int k = 0; k <net.getNoOfLinks();k++)
				{
					Link destLink = net.getAllLinks().elementAt(k);
					double prob = ((double)Edelete)*0.25*0.01*(double)net.getNoOfLinks()*((double)((destLink.iNodeLinks() )+( destLink.jNodeLinks())) / (double)net.getNoOfLinks());
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
					if(k >= (net.getNoOfLinks() -  Edelete)) //This is a BIG (hub-hub) link, should be deleted.
					{
						prob = 1.0;
					}
					else
					{
						prob = 0;
					}
					//prob = 0;
					 probs2[k] = prob;
				}
				
				
				
				for(int k = 0; k <net.getNoOfLinks();k++)
				{
					Link destLink = net.getAllLinks().elementAt(k);
					//Throw dice
					if(rand.nextDouble() <=  probs2[k])
					{
						try 
						{
							//Delete Link
							
							//Deleting immediately from network will cause problems with size, therefore mark and delete later.
							net.getAllLinks().elementAt(k).markForDeletion();
							Node iNode = net.getAllLinks().elementAt(k).iNode();
							Node jNode = net.getAllLinks().elementAt(k).jNode();
							iNode.deleteLink(jNode);
							jNode.deleteLink(iNode);
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
				
			
				//Now delete all the marked links
				net.removeMarkedLinks();
				 noOfLinks = net.getNoOfLinks();
				 int la = 0;
				
				//Now the links have been deleted and second mechanism completed
					System.out.println("Number of Nodes now  = " + net.getSize());
					System.out.println("Number of links now  = " + net.getAllLinks().size());
					System.out.println("DELETED Links in this node addition  = " + (linksNow - net.getAllLinks().size()));
				
			}
			
			
		}
		
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		return net;
	}
	
	public Network growASPA(Network net, int Npa)
	{
		int initialSize = net.getSize();
		System.out.println("Number of Nodes now  = " + net.getSize());
		System.out.println("Number of links now  = " + net.getAllLinks().size());
		Random rand  =  new Random(System.currentTimeMillis());
//		
		double EjoinSTART = 2.0;
		double EdeleteSTART =3.0002;
		double Ejoin = EjoinSTART; //Expectation of links created per node- a parameter
		double Edelete = EdeleteSTART; // Expected links deleted per node addition. A parameter
								
		for(int i = 1; i <= Npa; i++)
		{
			Ejoin = EjoinSTART; // + (0.02*(double)(net.getSize())); //increase addded links with growing network size
			Edelete = EdeleteSTART; // + (0.02*(double)(net.getSize())); //increase addded links with growing network size
			Node node = new Node(initialSize + i);
			double probs[] =  new double[net.getSize()];
			int noOfLinks = net.getNoOfLinks();
			
			for(int j = 0; j <net.getSize();j++)
			{
				Node destNode = net.getAllNodes().elementAt(j);
				double prob = Ejoin;
				if(net.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
				{
					 prob = 0.5*Ejoin*((double)destNode.getNumberOfLinks()) / (double)net.getNoOfLinks();
				}
				else //No links made yet
				{
					prob = 1.0;//If no links are made, this node node has to join with the existing isolated nodes
					//isolated nodes are not allowed to exist.
					
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
			
			 noOfLinks = net.getNoOfLinks();
				System.out.println("Number of links now  = " + net.getAllLinks().size());
				int linksNow = net.getAllLinks().size();
			
			//Now, stochastically delete some links that connect hubs
//				Now, stochastically delete some links that connect hubs
				
				double raw = 0.999; //A parameteter, determining whether stochastic deleted needs to be done
				//If this is set to 1.0, stochastic deleting is ALWAYS done.
				
				
				
				net.sortLinksBasedOnCombinedDegree();
				
				if( rand.nextDouble() <=  raw ) 
				{
					double probsum = 0;
					//destroy 'Edeleted' number of the existing links
					double probs2[] =  new double[net.getNoOfLinks()];
					for(int k = 0; k <net.getNoOfLinks();k++)
					{
						
						Link destLink = net.getAllLinks().elementAt(k);
						double prob = ((double)Edelete)*0.25*0.01*(double)net.getNoOfLinks()*((double)((destLink.iNodeLinks() )+( destLink.jNodeLinks())) / (double)net.getNoOfLinks());
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
						if(k >= (net.getNoOfLinks() -  (20.0*Edelete))) //This is a BIG (hub-hub) link, should be deleted.
							//ASSORTATIVE conenctions should be discouraged here; BOTH hub-to-hub as well as periphery-to-periphery
							//1+2+3+4+....+N = N*0.5*(N-1)
						{
							double rank = net.getNoOfLinks() - k;
							double weight = (20.0 - rank+ 1.0)/ (10.0*21.0);
							prob = weight;
							//prob = Edelete*2.0*(double)k / ((double)net.getNoOfLinks()*(double)net.getNoOfLinks()-1.0);
							//prob = 1.0;
							probsum = probsum + prob;
						}
						else
						{
							prob = 0;
							//prob = Edelete*2.0*(double)k / ((double)net.getNoOfLinks()*(double)net.getNoOfLinks()-1.0);
						}
						//prob = 0;
						 probs2[k] = prob;
						 
						 //randomly select some links from the 'assortative links' and replace them
					}
					
					
					
					for(int k = 0; k <net.getNoOfLinks();k++)
					{
						Link destLink = net.getAllLinks().elementAt(k);
						//Throw dice
						if(rand.nextDouble() <=  probs2[k])
						{
							try 
							{
								//Delete Link
								
								//Deleting immediately from network will cause problems with size, therefore mark and delete later.
								net.getAllLinks().elementAt(k).markForDeletion();
								Node iNode = net.getAllLinks().elementAt(k).iNode();
								Node jNode = net.getAllLinks().elementAt(k).jNode();
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
					net.replaceMarkedLinks(); //replace deleted links by connecting the corresponding nodes to other nodes selected preferentially,
					//so that a deleted Link does not result in breaking the network into components
					 noOfLinks = net.getNoOfLinks();
					 int la = 0;
					
					//Now the links have been deleted and second mechanism completed
						//System.out.println("Number of Nodes now  = " + net.getSize());
						//System.out.println("Number of links now  = " + net.getAllLinks().size());
					//System.out.println("DELETED Links in this node addition  = " + (linksNow - net.getAllLinks().size()));
					
				}
				
				
			}
			
			System.out.println("Number of Nodes now  = " + net.getSize());
			System.out.println("Number of links now  = " + net.getAllLinks().size());
			return net;
	}



}
