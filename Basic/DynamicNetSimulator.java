package Basic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Random;

public class DynamicNetSimulator {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	
		DynamicNetSimulator sim = new DynamicNetSimulator();
		sim.simulate();

	}
	
	public void simulate()
	{
		System.out.println("Executed class is DynamicNetSimulator -  SIM method");
		
		Random random = new Random(System.currentTimeMillis());
		NetworkGenerator gen = new NetworkGenerator();
		
		Network network;

		
		BufferedWriter out2 = null;
		BufferedReader in2 = null;
	
	
	
		//network = gen.generateDynamicNetFromFile();
		
		network  =  new Network();
		network.growNetworkPA(10, 1.85);
		
		network.createTextFile("Dynamic Net");
		

		System.out.println("Network size = "+ network.getSize()+ " links "+ network.getNoOfLinks());
		System.out.println("Network Gamma = "+ network.getGamma());
		
		int TIMESTEPS = 2;
		
		for(int loops =0; loops < TIMESTEPS; loops++)
		
		{
			
			//Store the prev. degrees
			for(int n=0; n< network.getSize(); n++)
			{
				TimeNode s =  (TimeNode)network.getAllNodes().elementAt(n);
				s.setPrevDegree(s.getNumberOfLinks());
			}
			
			
			
			
			
			//Simulate network so that changes are made in links
			if(loops <10)
			{
				network.growNetworkPA(50*loops + 50, 1.85);
			}
			else
			{
				network.growIGNet(100*loops + 100);
			}
			
		
			
			
			//Calculate dynamic information content
			
			double dynamicInfo = network.calculate_MI_Dynamic();
			double staticInfo = network.calculateMI();
			
			System.out.println("static MI "+ staticInfo+ " Dynamic MI "+ dynamicInfo);
			//System.out.println("static cond-ent "+ network.calculate_Hqc_2()+ " Dynamic cond-ent "+ network.calculate_Hqc_conditional_on_time());
			//System.out.println("ent "+ network.calculateHq()+ " ent t " + network.calculateHq_prev());
			System.out.println("Future link info: "+ (dynamicInfo - staticInfo) );
			System.out.println(" ");
			
	
		} 
		
		
		try {
//			Store the prev. degrees
			for(int n=0; n< network.getSize(); n++)
			{
				TimeNode s =  (TimeNode)network.getAllNodes().elementAt(n);
				s.setPrevDegree(s.getNumberOfLinks());
			}
			
			
		
				
		} catch (Exception e1) {
			System.out.println("Nodes in the network are not of type TimeNode");
			e1.printStackTrace();
		}
		
		

	} //END OF SIM METHOD
	
}
