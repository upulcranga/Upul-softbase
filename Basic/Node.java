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
    
    int birth_timestep; //The timestep in which this node has been added -  needed to analyze the 'history' of node in PA type scenarios
    protected double p_t; //variables to calculate pagerank
    protected double p_t_1;
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
   
        birth_timestep = 0;
        p_t_1 = 0;
        p_t = 0;
       
	
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
        
 
        birth_timestep = 0;
        p_t_1 = 0;
        p_t =0;
   
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

	@Override
	public int compareTo(Object arg0) {
		// TODO Auto-generated method stub
		return 0;
	}    

}
