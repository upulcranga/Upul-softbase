package Basic;

import java.util.ArrayList;
import java.util.Vector;

public class CityNode extends Node {
	
	private long houseID;
	private int age;
	private long schoolID;
	private long jobID;
	private long hospitalID;
	private long wardID;
	
	public CityNode(int id)
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
        
	    houseID=0;
	    age=0;
	    schoolID=0;
	    jobID=0;
	    hospitalID=0;
	    wardID=0;
	
	}
	
	   public long getSchoolID()
	    {
	        return schoolID;
	    }
	    
	    public void setSchoolID(long id)
	    {
	    	schoolID = id;
	
	    }
	    
	    public long getHouseID()
	    {
	        return houseID;
	    }
	    
	    public void setHouseID(long id)
	    {
	    	houseID = id;
	
	    }

	    
	    public long getHospitalID()
	    {
	        return hospitalID;
	    }
	    
	    public void setHospitalID(long id)
	    {
	    	wardID = id;
	
	    }

	    
	    public long getJobID()
	    {
	        return jobID;
	    }
	    
	    public void setJobID(long id)
	    {
	    	jobID = id;
	
	    }

	    
	    public long getAge()
	    {
	        return age;
	    }
	    
	    public void setAge(int a)
	    {
	    	age = a;
	
	    }


}
