package Basic;

public class SecurityAgent {
	
	final int LIMIT1  =  50;
	final int LIMIT2  =  40;
	final int LIMIT3  =  30;
	
	int numberOgAgents = 5;
	int agentID;
	
	double temp = 0;
	double raTemp = 0;
	
	double[] rat =  new double[12];
	int ratIndex = 0;
	
	public SecurityAgent(int ID)
	{
		agentID = ID;
		
	}
	
	public void updateTemp()
	{
		temp = getSensorReading();
		rat[ratIndex] = temp;
		ratIndex = (ratIndex+1)%12;
		
	}
	
	public void updateAverageTemp()
	{
		double sum = 0 ;
		for(int i=0; i< rat.length; i++)
		{
			sum = sum + rat[i];
		}
		raTemp =  sum / 12.0;
	}
	
	public void checkAlarm()
	{
		if(temp > LIMIT1)
		{
			soundAlarm();
		}
		
		if(temp > LIMIT2)
		{
			//	CHECK OTHER AGENTS
			
			int count = 0;
			
			for(int i=0; i< rat.length; i++)
			{
				int tempOther = getAverageTemperature(i);
				
				if(tempOther > LIMIT3)
				{
					count++;
				}
				
			}
			
			if(count > 3) //Note that including this agent the number is three
			{
				soundAlarm();
			}
		}
		
		
	}
	

}
