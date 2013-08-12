package Basic;

import java.util.Vector;

public class SIRNode extends Node {
	
	private String state; //Needs to be between S,I,R as strings. This is the 'current' state, synchronized or otherwise.
	public String synchState; // Also Needs to be between S,I,R as strings. 
									//Used to denote the 'synchronized' states while states are updated for the entire network
									//TODO: Need to make private, but the owning network must have access to change this state.
	
	
	public SIRNode()
	{
		super();
		state = "S";
		synchState =  "S";
	}
	
	public SIRNode(int id)
	{
		super(id);
		state = "S";
		synchState =  "S";
		
	
	}
	
	public SIRNode(int id, int d)
	{
		super(id,d);
		state = "S";
		synchState =  "S";
		
	
	}

	public boolean setSIRState(String st)
	{
		if(st.equals("S") || st.equals("I") || st.equals("R") || st.equals("V"))
		{
			state =  st;
			return true;
		}
		
		return false;
	}
	
	public String getSIRState()
	{
		return state;
	}
	
	
	/**
	 * Note that there is only a 'get' method to get the synchronized state.
	 * Synchronized states cannot be 'set' directly, but will be updated by a network level method when synchronization is causally appropriate.
	 * @return
	 */
	public String getSynchronizedSIRState()
	{
		return synchState;
	}
	
	/**
	 * Count the number of neighbour ndoes in 'I' state based on the synchronized states
	 * @return
	 */
	public int noOfInfectedNeoghbours()
	{
		int count = 0;
		 Vector<Node>v = this.getAllNeighbours();
		 for(int i=0; i<v.size(); i++)
		 {
			 SIRNode neighbour = (SIRNode)v.elementAt(i);
			 if(neighbour.getSynchronizedSIRState().equals("I"))
			 {
				 count++;
			 }
		 }
		return count;
	}
	
	public SIRNode clone()
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
	        
	    copy.BCCentrality = this.BCCentrality ;
		copy.RCCentrality = this.RCCentrality ;
		copy.CCCentrality = this.CCCentrality;
		copy.state = this.state;
		copy.synchState = this.synchState;
		return copy;
	}
}
