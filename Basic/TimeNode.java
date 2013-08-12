package Basic;

import java.util.ArrayList;
import java.util.Vector;

public class TimeNode extends Node {
	
	private long prevDegree; // The degree of this node in the previous timestep of the dynamic network
	private boolean  prevDegreeFlag; //A falg to state if prev. degree has been set
	
	
	public TimeNode(int id)
	{
		ID = id;
		indended_degree = 1200;
		links = new Vector<Link>();
		weight = 1;
		componentIndex = 0;
		componentAssigned = false;
		local_r = 0;
		BCCentrality =0;
		


		//Namal Senaratne's code
		nodeState = NodeState.NOT_ENQUEUED;
        adjacency = new ArrayList<Integer>();
        predecessor = new ArrayList<Integer>();
        
        
        nodeReading = 0;
        riskFactor=0.0;
        
        prevDegree = 0;
        prevDegreeFlag = false;
	
	}
	
	  
	    public void setPrevDegree(long d)
	    {
	    	 prevDegree =  d;
	    	 prevDegreeFlag = true;
	
	    }
	    
	    public long getPrevDegree()
	    {
	        return prevDegree;
	    }
	    
	    public boolean isPrevDegreeSet()
	    {
	        return prevDegreeFlag;
	    }


}
