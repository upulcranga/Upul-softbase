/**
 * @author Piraveenan
 */
package Basic;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Vector;

public class Node implements Comparable {

	protected Integer ID;
	protected int indended_degree;
	protected Vector<Link> links;
	protected Random random;
	protected double weight; //For weighted models
	int componentIndex; //to indicate membership of component in fractured networks
	protected boolean componentAssigned;
	protected double local_r; //local assortativity can be calculated and assigned to this
	protected double BCCentrality;  //Betweenness centrality can be calculated and assigned to this variable
	protected double RCCentrality;  //Risk centrality can be calculated and assigned to this variable
	protected double CCCentrality;  //Closeness centrality can be calculated and assigned to this variable
	protected double hopDistance; // The number of hops from a noe with a certain property, such as an infected state
	//NB: The above three quantities are NOT calculated during initialization of Network. They need to be calculated and assigned 
	//explicitly before being used.
	
	/* Namal senaratne's code */
	protected int distance; //This node's distance from a pre difined node
	protected List<Integer> predecessor;
	protected NodeState nodeState; // * Current state of the node
	protected List<Integer> adjacency; //The list which contains all the adjacent nodes of this node
	
	/*variables for network dynamics*/
    double nodeReading; //a generic variable to indicate sensor reading, expression level etc
    double riskFactor; //A variable that indicates a node's 'Risk Level' in decease networks
    int birth_timestep; //The timestep in which this node has been added -  needed to analyze the 'history' of node in PA type scenarios
    protected double p_t; //variables to calculate pagerank
    protected double p_t_1;
    protected int noOfStalkers; //variables to calculate stalking coefficients
    protected double noOfVictims; 
    public Vector<Node> victimSet =  new Vector<Node>(); //A set of victims that this node intentionally stalks  //TODO: Make this private and have a proper mechanism for adding nodes
    
    //Darshana's code for Game simulation
    //Variable names modified an incorporated by Piraveen on 07.05.2013
    public boolean isCoordinator = false; //previous variable name 'Collaborator': denotes if a node is coordinator or non-coordinator in a coordinatio game
  //  public Vector<Double> coordinationReturns = new Vector<Double>(); //previous variable name 'oordinationReturns'
    public double gameTotalPayoff = 0; //previous variable name 'coordinationValue': the total pay-off a player has accumulated during the evolutionary game process
    public int payOffLag = 20; //the lag time after which the payoffs are taken into account
    public double[] payOffBuffer =  new double[payOffLag+1]; //The extra cell is given so that the curren value can be stored
    public double researchAbility = generateRandom(0,100); //for authorNet
    public double percentile = 0; //percentile rank
    //public double impactFactor = 0; //for papernet
	public int numberOfCitations = generateRandom(10,50); //for papernet
	public Integer [] authorArray = new Integer[generateRandom(1,5)]; //for papernet
	public ArrayList<Integer> papers = new ArrayList<Integer>();//papers authored/co-authored by each author
	public int hIndex = 0;
    
    
    public Node()
	{
    	
	}
    
    public Node(int id)
	{
		ID = id;
		indended_degree = 1000;
		links = new Vector<Link>();
		weight = 1;
		componentIndex = 0;
		componentAssigned = false;
		local_r = 0;
		BCCentrality =0;
		RCCentrality =0;
		CCCentrality=0;
		random = new Random(System.currentTimeMillis());
		
		//Namal Senaratne's code
		nodeState = NodeState.NOT_ENQUEUED;
        adjacency = new ArrayList<Integer>();
        predecessor = new ArrayList<Integer>();
        
        nodeReading = 0;
        riskFactor=0.0;
        birth_timestep = 0;
        p_t_1 = 0;
        p_t = 0;
        noOfStalkers=0; 
        noOfVictims=0;
	
	}
	

    
	public Node(int id, int d)
	{
		ID = id;
		indended_degree = d;
		links = new Vector<Link>();
		weight = 1;
		componentIndex = 0;
		componentAssigned = false;
		local_r = 0;
		BCCentrality =0;
		RCCentrality =0;
		CCCentrality=0;
		random = new Random(System.currentTimeMillis());
		
		//Namal Senaratne's code
		nodeState = NodeState.NOT_ENQUEUED;
        adjacency = new ArrayList<Integer>();
        predecessor = new ArrayList<Integer>(); 
        
        nodeReading = 0;
        riskFactor=0.5;
        birth_timestep = 0;
        p_t_1 = 0;
        p_t =0;
        noOfStalkers=0; 
        noOfVictims=0;
	}
	
	public Integer getID()
	{
		return ID;
	}
	
	public void reassignID(int id)
	{
		ID =  id;
	}
	
	public void setBirthTimeStep(int ts)
	{
		birth_timestep =  ts;
	}
	
	public int getBirthTimeStep()
	{
		return birth_timestep;
	}
	
	public int getAssignedDegree()
	{
		return indended_degree;
	}
	
	public int getAssignedRemainingDegree()
	{
		return (indended_degree-1);
	}
	
	public void addLink(Node node, boolean source)
	{
		if(source)
		{
			//this node is the source node
			links.addElement(new Link(this,node));
		}
		else
		{
			//this node is the target node
			links.addElement(new Link(node,this));
		}
		
		//degree++;
	}
	
	public void deleteLink(Node node)
	{
		//Delete a link with this node
		for(int i=0; i< links.size(); i++)
		{
			if( (links.elementAt(i).iNodeID()== node.getID()) || (links.elementAt(i).jNodeID()== node.getID() ) )
			{
				links.removeElementAt(i);  //delete link
				break;
			}
		}
	}
	
	
	
	public Vector<Link> getLinks()
	{
		return links;
	}
	
	public Node getRandomNeighbour()
	{
		if(links.size()==0)
		{
			return this;
		}
		Link randomLink =  links.elementAt(Math.abs(random.nextInt()% links.size()));
		Node nei =  randomLink.iNode();
		if(nei.getID()== this.getID())
		{
			nei =  randomLink.jNode();
		}
		
		return nei; //TODO: be aware that because of self links, 'this' node might be returned as random neighbour still.
	}
	
	
	public int getNumberOfInLinks()
	{
		int no = 0;
		

		for(int i=0; i< links.size(); i++)
		{
			Node thisNode = links.elementAt(i).jNode();
			
			if(thisNode.getID() == ID)
			{
				//This node is destination node for this Link
				no++;
			}
			
		}
		
		return no;
	}
	
	public int getNumberOfOutLinks()
	{
		int no = 0;
		
		for(int i=0; i< links.size(); i++)
		{
			Node thisNode = links.elementAt(i).iNode();
			
			if(thisNode.getID() == ID)
			{
				//This node is source node for this Link
				no++;
			}
			
		}
		
		return no;
	}
	
	public int getNumberOfLinks()
	{
		return links.size();
	}
	
	/**
	 * In this function, self-links and double links are not counted
	 * @return
	 */
	public Vector<Node> getAllNeighbours()
	{
		Vector<Node> neighbours =  new Vector<Node>();
		for(int i=0; i< links.size(); i++)
		{
			Link thisLink = links.elementAt(i);
			if(thisLink.iNodeID()== thisLink.jNodeID()) //self link
			{
				// Dont do anything
				continue;
				
			}
			
			//Identify the neighbour in this Link
			Node nei1 =  thisLink.iNode();
			if(nei1.getID() == this.getID())
			{
				nei1 = thisLink.jNode(); //Now we have identified the neighbour in this link
			}
			
			boolean alreadyThere  =  false;
		
			
			for(int j=0; j< i; j++) //check already traversed elements for double entry
			{
				Link otherLink = links.elementAt(j);
				Node nei2 =  otherLink.iNode();
				if(nei2.getID() == this.getID())
				{
					nei2 = otherLink.jNode(); //Now we have identified the neighbour in this link
				}
				
				//Is this neighbour already identified?
				if(nei1.getID() ==nei2.getID()) //double link
				{
					alreadyThere =  true;
					break;
				}
			}
			
			if(!alreadyThere)
			{
				//nei1 is not already there
				neighbours.addElement(nei1);
			}
			
		}
		return neighbours;
	}
	
	public int generateRandom(int min, int max){
		return min + (int)(Math.random() * ((max - min) + 1)); //return a random number of citations to be initiated.
	}

	/**
	 * In this function, self-links and double links are not counted
	 * @return
	 */
	public int getNumberOfNeighbours()
	{
		int noOfNeighbours = links.size();
		for(int i=0; i< links.size(); i++)
		{
			Link thisLink = links.elementAt(i);
			if(thisLink.iNodeID()== thisLink.jNodeID()) //self link
			{
				noOfNeighbours--;
				continue;
			}
			
			//Identify the neighbour in this Link
			Node nei1 =  thisLink.iNode();
			if(nei1.getID() == this.getID())
			{
				nei1 = thisLink.jNode(); //Now we have identified the neighbour in this link
			}
		
			
			for(int j=0; j< i; j++) //check already traversed elements for double entry
			{
				Link otherLink = links.elementAt(j);
				Node nei2 =  otherLink.iNode();
				if(nei2.getID() == this.getID())
				{
					nei2 = otherLink.jNode(); //Now we have identified the neighbour in this link
				}
				
				//Is this neighbour already identified?
				if(nei1.getID() ==nei2.getID()) //double link
				{
					noOfNeighbours--;
					break;
				}
			}
			
		}
		return noOfNeighbours;
	}
	

	/**
	 * This function gives the average DIFFERENCE IN DEGREES between a node and its neighbours.
	 * Self links are NOT counted. Double links ARE counted.
	 * The difference is taken as a modulus.
	 * @return
	 */
	public double getAverageNeighbourDifference()
	{
		double dbar = 0;
		double ind_links = 0;
	
		for(int i=0; i< links.size(); i++)
		{
			Link thisLink = links.elementAt(i);
			if(thisLink.iNodeID()== thisLink.jNodeID()) //self link
			{
				System.out.println("Self Links in new NA methods");
				continue;
			}
			
			//Identify the neighbour in this Link
			int d1 =  thisLink.iNode().getNumberOfLinks();
			int d2 =  thisLink.jNode().getNumberOfLinks();
			
			dbar = dbar + Math.abs(d1-d2);
			ind_links++;
		
			
			
			
		}
		
		if(ind_links==0)
		{
			System.out.println("ZERO RETURNED in new NA methods");
			return 0;
			
		}
		
		 dbar =  dbar / ind_links;
		
		return dbar;
	}
	
	/**
	 * This method returns a link list without self links and double links
	 */
	
	public Vector<Link> returnUniqueLinkList()
	{
		 Vector<Link> l = new Vector<Link>();
		 boolean noSelfLink =  true;
		 for(int i=0; i<links.size();i++)
		 {
			 noSelfLink =  true;
			 Link thisLink = links.elementAt(i);
			 if(thisLink.iNodeID()== thisLink.jNodeID()) //self link
				{
					// do not add to the new list
				 	continue;
				}
			 
				//Identify the neighbour in this Link
				Node nei1 =  thisLink.iNode();
				if(nei1.getID() == this.getID())
				{
					nei1 = thisLink.jNode(); //Now we have identified the neighbour in this link
					if(nei1.getID() == this.getID()) //Just a double check
					{
						System.out.println("SELF LINK: Shouldn't have occured here");
					}
				}
				
				for(int j=0; j< i; j++) //check already traversed elements for double entry
				{
					Link otherLink = links.elementAt(j);
					Node nei2 =  otherLink.iNode();
					if(nei2.getID() == this.getID())
					{
						nei2 = otherLink.jNode(); //Now we have identified the neighbour in this link
						if(nei2.getID() == this.getID()) //Just a double check
						{
							System.out.println("SELF LINK: Shouldn't have occured here");
						}
					}
					
					//Is this neighbour already identified?
					if(nei1.getID() == nei2.getID()) //double link
					{
						noSelfLink =  false;
						break;
					}
				}
				
				//Self-links have been checked for ith element
				if(noSelfLink)
				{
					l.addElement(links.elementAt(i));
				}
				
			 
		 }
		 
		 //Now all links which were not self links and double links have been copied into l
		 
		 return l;
	}

	public void assignAuthors(Integer [] authors){
		authorArray = authors;
	}
	 
	public void assignWeight(double w){
		weight = w;
	}
	public double getWeight()
	{
		return weight;
	}
	public Node clone()
	{
		 Node copy = new Node((int)this.getID(), this.getAssignedDegree()); //It is O.K to shallow copy logic function, as it cannot be changed once the ndoe is created
		 copy.indended_degree = this.indended_degree;
		 copy.assignWeight(this.weight);
		 //Vector<Link> l = new Vector<Link>(this.links.size());
		 // copy the links from the original ndoe
		// for(int y=0; y<links.size();y++)
		// {
		//	 l.addElement(((DirectedLink)links.elementAt(y)).clone());
		// }
		// copy.links = l;
		 
		 //Namals' code
		 
		 copy.setDistance(this.distance);
		 copy.nodeState = this.nodeState;
		 copy.setReading(this.getReading());
	        
	        // creating the new adjacency copy
	        Integer val;
	        Integer newVal;
	        for(int i = 0; i < adjacency.size(); i++)
	        {
	            val = adjacency.get(i);
	            newVal = new Integer(val.intValue());
	            copy.adjacency.add(newVal);
	        }
	        
	        // creating the new predecessor copy
	        for(int j = 0; j < predecessor.size(); j++)
	        {
	            val = predecessor.get(j);
	            newVal = new Integer(val.intValue());
	            copy.predecessor.add(newVal);
	        }
	        
	    copy.BCCentrality = this.BCCentrality ;
		copy.RCCentrality = this.RCCentrality ;
		copy.CCCentrality = this.CCCentrality;
			
		 return copy;
	}
	
	public SIRNode cloneToSIRNode()
	{
		 SIRNode copy = new SIRNode((int)this.getID(), this.getAssignedDegree()); //It is O.K to shallow copy logic function, as it cannot be changed once the ndoe is created
		 copy.indended_degree = this.indended_degree;
		 copy.assignWeight(this.weight);
		 //Vector<Link> l = new Vector<Link>(this.links.size());
		 // copy the links from the original ndoe
		// for(int y=0; y<links.size();y++)
		// {
		//	 l.addElement(((DirectedLink)links.elementAt(y)).clone());
		// }
		// copy.links = l;
		 
		 //Namals' code
		 
		 copy.setDistance(this.distance);
		 copy.nodeState = this.nodeState;
		 copy.setReading(this.getReading());
	        
	        // creating the new adjacency copy
	        Integer val;
	        Integer newVal;
	        for(int i = 0; i < adjacency.size(); i++)
	        {
	            val = adjacency.get(i);
	            newVal = new Integer(val.intValue());
	            copy.adjacency.add(newVal);
	        }
	        
	        // creating the new predecessor copy
	        for(int j = 0; j < predecessor.size(); j++)
	        {
	            val = predecessor.get(j);
	            newVal = new Integer(val.intValue());
	            copy.predecessor.add(newVal);
	        }
		 return copy;
	}
	
	public boolean isLinked(Node node)
	{
		boolean b = false;
		for(int y=0; y<links.size();y++)
		{
			if(links.elementAt(y).iNodeID() == node.getID())
			{
				b = true;
			}
			
			if(links.elementAt(y).jNodeID() == node.getID())
			{
				b = true;
			}
				
		}
		
		return b;
	}
	
	public int compareTo222(Object on)
	{
		Node otherNode = (Node)on;
		
		if(getNumberOfLinks() < otherNode.getNumberOfLinks())
		{
			return 1;
		}
		else 
		{
			return -1;
		}
		
	}
	
	/**
	 * Node ID based comparison
	 * @param on
	 * @return
	 */
	
	public int compareTo6(Object on)
	{
		Node otherNode = (Node)on;
		
		if(getID() < otherNode.getID())
		{
			return 1;
		}
		else 
		{
			return -1;
		}
		
	}
	
	/**
	 * Local assortativity based comparison
	 * WARNING: local assortativity must be assigned before comparing!!
	 * @param on
	 * @return
	 */
	// > should not be the case unless reverse order desired!!!!
	public int compareTo5(Object on)
	{
		Node otherNode = (Node)on;
		
		if( Math.abs(getLocalR()) <  Math.abs(otherNode.getLocalR()))
		//if( getLocalR() < otherNode.getLocalR())
		{
			return 1;
		}
		else 
		{
			return -1;
		}
		
	}
	
	/**
	 * Betweenness centrality based comparison
	 * WARNING: BC must be assigned before comparing!!
	 * @param on
	 * @return
	 */
	
	public int compareTo(Object on)
	{
		Node otherNode = (Node)on;
		
		if(Math.abs(getBCCentrality()) < Math.abs(otherNode.getBCCentrality()))
		{
			return 1;
		}
		else 
		{
			return -1;
		}
		
	}
	
	/**
	 * Risk centrality based comparison
	 * WARNING: Risk Centrality must be assigned before comparing!!
	 * @param on
	 * @return
	 */
	
	public int compareTo4444(Object on)
	{
		Node otherNode = (Node)on;
		
		if(Math.abs(getRCCentrality()) < Math.abs(otherNode.getRCCentrality()))
		{
			return 1;
		}
		else 
		{
			return -1;
		}
		
	}
	
	public int compareTo22(Object on)
	{
		Node otherNode = (Node)on;
		
		if(Math.abs(getHopDist()) < Math.abs(otherNode.getHopDist()))
		{
			return 1;
		}
		else 
		{
			return -1;
		}
		
	}
	
//	FROM HERE ONWARDS, THE FUNCTIONS IMPLEMENTED ARE MODIFIED VERSIONS OF NAMAL SENARATNE'S CODE
	//These functions are needed to implement shortest path calculations and betweenness centrality calculations

	/**
     * Returns the distance of the node from a pre-specified source node in the 
     * graph
     *
     * @throw   the distance
     */
    public int getDistance() 
    {
        return distance;
    }

    /**
     * Sets the distance of the node from a pre-specified source node in the 
     * graph
     *
     * @param   distance    the distance calculated
     */
    public void setDistance(int distance) 
    {
        this.distance = distance;
    }
    

    /**
     * Returns the list of predecessor nodes of this node
     *
     * @throw   a list which holds Integers which specify the IDs of the predecessor
     *          nodes
     */
    public List<Integer> getPredecessor()
    {
        return predecessor;
    }
    
    /**
     * Add a newly identified predecessor to this node's predecessor list
     *
     * @param   NodeID  the ID of the new predecessor node
     */
    public void addPredecessorNode(int NodeID)
    {
        Integer nodeId = new Integer(NodeID);
        if(!predecessor.contains(nodeId))
        {
            predecessor.add(nodeId);
        }
    }
    
    /**
     * Clears the predecessor list
     */
    public void clearPredecessors()
    {
        predecessor.clear();
    }
    
    /**
     * Sets the node state
     *
     * @param   nodeState   the state of the node defined by NodeState
     */
    public void setNodeState(NodeState nodeState)
    {
        this.nodeState = nodeState;
    }
    
    /**
     * Returns this nodes state
     *
     * @return  the NodeState enum object specifying the node states
     */
    public NodeState getNodeState()
    {
        return nodeState;
    }
    

    /**
     * Returns the list of adjacent nodes of this node
     *
     * @throw   a list which holds Integers which specify the IDs of the adjacent
     *          nodes
     */
    public ArrayList<Integer> getAdjacency() 
    {
        return (ArrayList)adjacency;
    }

    /**
     * Add a newly identified adjacenct node to this node's adjacency list
     *
     * @param   nodeID  the ID of the new adjacent node
     */
    public void addAdjacendyNode(int nodeID)
    {
        Integer nodeId = new Integer(nodeID);
        if(!adjacency.contains(nodeId))
        {
            adjacency.add(nodeId);
        }
    }
    
    
    public int getComponentIndex()
    {
        return componentIndex;
    }
    
    public void setComponentIndex(int index)
    {
        componentIndex = index;
        componentAssigned = true;
    }
    
    public boolean isComponentAssigned()
    {
    	return componentAssigned;
    }
    
    public void resetComponentAssignedStatus()
    {
    	componentIndex  =0;
        componentAssigned = false;
    }
    
    
    public double getLocalR()
    {
        return local_r;
    }
    
    public void setLocalR(double rho)
    {
    	local_r = rho;
    }
    
    public void clearLocalR()
    {
    	local_r = 0;
    }
    
    public double getRCCentrality()
    {
        return RCCentrality;
    }
    
    public void setRCCentrality(double bc)
    {
    	RCCentrality = bc;
    }
    
    
    public double getHopDist()
    {
        return hopDistance;
    }
    
    public void setHopDist(double d)
    {
    	hopDistance = d;
    }
    
    
    public void clearRCCentrality()
    {
    	RCCentrality = 0;
    }
    
    public double getBCCentrality()
    {
        return BCCentrality;
    }
    
    public void setBCCentrality(double bc)
    {
    	BCCentrality = bc;
    }
    
    
    public double getCCCentrality()
    {
        return CCCentrality;
    }
    
    public void setCCCentrality(double cc)
    {
    	CCCentrality = cc;
    }
    
    
    public void clearBCCentrality()
    {
    	BCCentrality = 0;
    }
    
    public void clearCCCentrality()
    {
    	CCCentrality = 0;
    }
    
    
    /**
     * For consistency: a function is provided which treats degree as 'degree centrality'
     * @return
     */
    public double getDCCentrality()
    {
        return this.getNumberOfLinks();
    }
    
 
    public void setReading(double x)
    {
    	nodeReading = x;
    }
    
    public double getReading()
    {
    	return nodeReading;
    }
    
    public void setRiskLevel(double x)
    {
    	riskFactor = x;
    }
    
    public double getRiskLevel()
    {
    	return riskFactor;
    }
    
    /**
     * This function gets the 'average' of neighbour readings of a given node
     * @return
     */
    
    public double getAverageNeighbourReading()
    {
    	double nbar = 0;
    	double count = 0;
    	

		for(int i=0; i< links.size(); i++)
		{
			if( (links.elementAt(i).iNodeID()!= this.getID())  )
			{
				//this is a neghbour
				Node neighbourNode = links.elementAt(i).iNode();
				nbar =  nbar + neighbourNode.getReading();
				count++;
			}
			
			if(  (links.elementAt(i).jNodeID()!= this.getID() ) )
			{
				//this is a neghbour
				Node neighbourNode = links.elementAt(i).jNode();
				nbar =  nbar + neighbourNode.getReading();
				count++;
			}
		}
		
		nbar =  nbar / count;
    	return nbar;
    }
    
    /**
     * Fimd the smallest neighbour in terms of degrees
     * In addition: The 'smallest' neighbour should have a reading abigger than a THRESHOLD
     * @return
     */
    public Node getSmallestNeighbour()
    {
    	double SmallestNeighbourDegree = 100000; //set it to very high value at the beginning
    	Node sn =  this;
    	
    	double THRESHOLD = 0.5;

		for(int i=0; i< links.size(); i++)
		{
			if( (links.elementAt(i).iNodeID()!= this.getID())  )
			{
				//this is a neghbour
				Node neighbourNode = links.elementAt(i).iNode();
				if(neighbourNode.getNumberOfLinks() < SmallestNeighbourDegree)
				{
					if(neighbourNode.getReading() > THRESHOLD)
					{
						SmallestNeighbourDegree = neighbourNode.getNumberOfInLinks();
						sn = neighbourNode;
					}
				}
			}
			
			if(  (links.elementAt(i).jNodeID()!= this.getID() ) )
			{
				//this is a neghbour
				Node neighbourNode = links.elementAt(i).jNode();
				if(neighbourNode.getNumberOfLinks() < SmallestNeighbourDegree)
				{
					if(neighbourNode.getReading() > THRESHOLD)
					{
						SmallestNeighbourDegree = neighbourNode.getNumberOfInLinks();
						sn = neighbourNode;
					}
				}
			}
		}
		
	
    	return sn;
    }
    
    /**
     * Fimd the biggest neighbour in terms of degrees
     * In addition: The 'biggest' neighbour should have a reading bigger than a THRESHOLD
     * @return
     */
    public Node getBiggestNeighbour()
    {
    	double BiggestNeighbourDegree = 0; //set it to very high value at the beginning
    	Node sn =  this;
    	
    	double THRESHOLD = -0.5;

		for(int i=0; i< links.size(); i++)
		{
			if( (links.elementAt(i).iNodeID()!= this.getID())  )
			{
				//this is a neghbour
				Node neighbourNode = links.elementAt(i).iNode();
				if(neighbourNode.getNumberOfLinks() > BiggestNeighbourDegree)
				{
					if(neighbourNode.getReading() > THRESHOLD)
					{
						BiggestNeighbourDegree = neighbourNode.getNumberOfInLinks();
						sn = neighbourNode;
					}
				}
			}
			
			if(  (links.elementAt(i).jNodeID()!= this.getID() ) )
			{
				//this is a neghbour
				Node neighbourNode = links.elementAt(i).jNode();
				if(neighbourNode.getNumberOfLinks() > BiggestNeighbourDegree)
				{
					if(neighbourNode.getReading() > THRESHOLD)
					{
						BiggestNeighbourDegree = neighbourNode.getNumberOfInLinks();
						sn = neighbourNode;
					}
				}
			}
		}
		
	
    	return sn;
    }


    /**
     * This function checks if this node has a 'bigger' neighbour in terms of degree.
     * If there is a bigger neighbour, it will return true; otherwise it will return false.
     */
    public boolean isBiggerNeighbour()
    {
    	boolean bool = true;
    	

		for(int i=0; i< links.size(); i++)
		{
			if( (links.elementAt(i).iNodeID()!= this.getID())  )
			{
				//this is a neghbour
				Node neighbourNode = links.elementAt(i).iNode();
				if(neighbourNode.getNumberOfLinks() > this.getNumberOfLinks())
				{
					bool = true;
				}
			}
			
			if(  (links.elementAt(i).jNodeID()!= this.getID() ) )
			{
				//this is a neghbour
				Node neighbourNode = links.elementAt(i).jNode();
				if(neighbourNode.getNumberOfLinks()> this.getNumberOfLinks())
				{
					bool = true;
				}
			}
		}
		
	
    	return bool;
    }
    
    public void setNoOfStalkers(int n)
    {
    	noOfStalkers = n;
    }
    
    public double getNoOfStalkers()
    {
    	return noOfStalkers;
    }
    
    
    public void setNoOfVictims(int n)
    {
    	noOfVictims = n;
    }
    
    public double getNoOfVictims()
    {
    	return noOfVictims;
    }

    public double getPageRank()
    {
    	return p_t;
    }
    
    public double setPageRank(double r)
    {
    	return p_t = r;
    }
    
    public void updatePageRank(double alpha, int N)
    {
    	
    	double sum = 0;
    	for(int i=0; i< this.getLinks().size(); i++)
    	{
    		Link thisLink = this.getLinks().elementAt(i);
    		Node neighNode =  thisLink.iNode();
    		
    		if(neighNode.getID()== this.getID())
    		{
    			neighNode =  thisLink.jNode();
    			
    		}
    		if(neighNode.getNumberOfOutLinks() > 0)
    		{
    			sum =  sum + (neighNode.p_t_1 / (double)neighNode.getNumberOfOutLinks());
    		}
    	}
    	
    	sum = alpha*sum;
    	
    	sum =  sum + ((1-alpha)/(double)N);
    	
    
    	this.p_t =  sum;
    }
    
    //GAME THEORY methods
    
    public double[] getPayoffBuffer()
    {
    	return payOffBuffer;
    }
    
    public void updatePayoffBuffer(double newVal)
    {
 
    	
    	for(int i=(payOffBuffer.length-1); i>=1; i--)
    	{
    		payOffBuffer[i] = payOffBuffer[i-1];
    	}
    	payOffBuffer[0] = newVal;
    }
    
    

}
