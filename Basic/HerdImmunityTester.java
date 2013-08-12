package Basic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.Vector;

public class HerdImmunityTester {
	
	public static void main(String[] args) {

		
		HerdImmunityTester sim = new HerdImmunityTester();
	
		 sim.mul_simulatespread();

		
	}
	
	public void mul_simulatespread()
	{
		System.out.println("Executed class is HerdImmunityTester");
		
		BufferedWriter out = null;
		BufferedWriter out_g20 = null;
		BufferedWriter out_g50 = null;
		BufferedWriter out_g100 = null;
		try {
			out = new BufferedWriter(new FileWriter(".\\..\\Output\\HI\\smallworld_1000_4_2- THIRD.xls"));
			//out_g20 = new BufferedWriter(new FileWriter(".\\..\\Output\\HI\\SW nat imm g20 t 31.xls"));
			//out_g50 = new BufferedWriter(new FileWriter(".\\..\\Output\\HI\\SW nat imm g50 t 31.xls"));
			//out_g100 = new BufferedWriter(new FileWriter(".\\..\\Output\\HI\\SW nat imm g100 t 31.xls"));
		}
			
		 catch (IOException e) {
			System.out.println("cant write to file");
			e.printStackTrace();
		}

		NetworkGenerator gen = new NetworkGenerator();
		SIRNetwork network;
		//Network network;
		//Network net =  gen.growNetworkPA(1000, 10.05);
		//Network net =  gen.generateERRandom(1000,2000);
		//Network net =  gen.generateSmallWorld(200, 56, 0.5);
		//Network net = gen.growNetworkAPA(5000,200,2,-0.5 );
		//net.createTextFile("sw 200 4000");
		//System.exit(0);
		
		//network = (gen.growNetworkPA(5000, 1)).cloneToSIRNet();
		//network = (gen.generateERRandom(1000,2400)).cloneToSIRNet();
		//network = (gen.generateERRandom(1000,2400));
		network = gen.generateFromFile().cloneToSIRNet();
		
		System.out.println("Network size = "+ network.getSize()+ " links "+ network.getNoOfLinks());
		System.out.println("Network Gamma = "+ network.getGamma());
	
		for(int k=0; k<100; k++)
		{
			//network = simulatespread_firstspread(network);
			network = immunizeNetwork_BC(network, (k+1));
			//network =  immunizeNetwork_Random(network,(double)(k+1)*0.1);
			network.forceRecoveryofAllNodes(); //To make sure that no nodes are now in 'I" state
			network.confirmVaccination();
			double vacPer=0;
			double avgSurPer=0;
				for(int y=0; y<20; y++)
					{
						Vector v = simulatespread_secondspread(network);
						//Print the contents of vector v into a file
						vacPer = ((Percentages)v.elementAt(0)).immunizedPercentage;
						double surPer = ((Percentages)v.elementAt(0)).survivingPercentage;
						avgSurPer = avgSurPer + surPer;
						if( (surPer - vacPer) >=80)
							
						{
							System.out.println("EFFECTIVE NATURAL IMMUNIZATION");
						}

						network.SynchronizeAllSIRStates();
						network.reInititializeNonVaccinatedNodes();
					}
				
				avgSurPer =  avgSurPer / 20;
				try {
					out.write(""+vacPer+"\t"+avgSurPer);
					out.flush();
					out.write( System.getProperty(("line.separator")));
					out.flush();
				} catch (IOException e) {
					System.out.println("Could not write percentage outputs to the file");
					e.printStackTrace();
				}
				
//				//Now grow network by a certain percentage and reproduce the results
//				//*******************************************************************
//				double growthPercent = 20;
//				int n0 = (int)(0.01*network.getSize()*growthPercent);
//				int m0 = 2*n0;
//				
//				SIRNetwork gnetwork = (SIRNetwork)network.clone();
//				gnetwork = gnetwork.growNetwork_SW(n0, 4, 0.5);
//				
//				for(int y=0; y<20; y++)
//				{
//					Vector v = simulatespread_secondspread(gnetwork);
//					//Print the contents of vector v into a file
//					vacPer = ((Percentages)v.elementAt(0)).immunizedPercentage;
//					double surPer = ((Percentages)v.elementAt(0)).survivingPercentage;
//					avgSurPer = avgSurPer + surPer;
//					if( (surPer - vacPer) >=80)
//						
//					{
//						System.out.println("EFFECTIVE NATURAL IMMUNIZATION");
//					}
//
//					gnetwork.SynchronizeAllSIRStates();
//					gnetwork.reInititializeNonVaccinatedNodes();
//					}
//					
//					avgSurPer =  avgSurPer / 20;
//					try {
//						out_g20.write(""+vacPer+"\t"+avgSurPer);
//						out_g20.flush();
//						out_g20.write( System.getProperty(("line.separator")));
//						out_g20.flush();
//					} catch (IOException e) {
//						System.out.println("Could not write percentage outputs to the file");
//						e.printStackTrace();  
//					}
//			
//				
//					
//				//Growing part of code completed
//				//*********
//				//Repeat
//					 growthPercent = 50;
//					 n0 = (int)(0.01*network.getSize()*growthPercent);
//					 m0 = 2*n0;
//					
//					 gnetwork = (SIRNetwork)network.clone();
//					 gnetwork = gnetwork.growNetwork_SW(n0, 4, 0.5);
//					
//					for(int y=0; y<20; y++)
//					{
//						Vector v = simulatespread_secondspread(gnetwork);
//						//Print the contents of vector v into a file
//						vacPer = ((Percentages)v.elementAt(0)).immunizedPercentage;
//						double surPer = ((Percentages)v.elementAt(0)).survivingPercentage;
//						avgSurPer = avgSurPer + surPer;
//						if( (surPer - vacPer) >=80)
//							
//						{
//							System.out.println("EFFECTIVE NATURAL IMMUNIZATION");
//						}
//
//						gnetwork.SynchronizeAllSIRStates();
//						gnetwork.reInititializeNonVaccinatedNodes();
//						}
//						
//						avgSurPer =  avgSurPer / 20;
//						try {
//							out_g50.write(""+vacPer+"\t"+avgSurPer);
//							out_g50.flush();
//							out_g50.write( System.getProperty(("line.separator")));
//							out_g50.flush();
//						} catch (IOException e) {
//							System.out.println("Could not write percentage outputs to the file");
//							e.printStackTrace();
//						}
//				
//					
//						
//					//Growing part of code completed
//					//*********
//						
//						//Repeat
//						 growthPercent = 100;
//						 n0 = (int)(0.01*network.getSize()*growthPercent);
//						 m0 = 2*n0;
//						
//						 gnetwork = (SIRNetwork)network.clone();
//						 gnetwork = gnetwork.growNetwork_SW(n0, 4, 0.5);
//						
//						for(int y=0; y<20; y++)
//						{
//							Vector v = simulatespread_secondspread(gnetwork);
//							//Print the contents of vector v into a file
//							vacPer = ((Percentages)v.elementAt(0)).immunizedPercentage;
//							double surPer = ((Percentages)v.elementAt(0)).survivingPercentage;
//							avgSurPer = avgSurPer + surPer;
//							if( (surPer - vacPer) >=80)
//								
//							{
//								System.out.println("EFFECTIVE NATURAL IMMUNIZATION");
//							}
//
//							gnetwork.SynchronizeAllSIRStates();
//							gnetwork.reInititializeNonVaccinatedNodes();
//							}
//							
//							avgSurPer =  avgSurPer / 20;
//							try {
//								out_g100.write(""+vacPer+"\t"+avgSurPer);
//								out_g100.flush();
//								out_g100.write( System.getProperty(("line.separator")));
//								out_g100.flush();
//							} catch (IOException e) {
//								System.out.println("Could not write percentage outputs to the file");
//								e.printStackTrace();
//							}
//					
//						
							
						//Growing part of code completed
						//*********
						
						
			
			network.reInititializeAllNodes();
			//System.out.println("net size now "+gnetwork.getSize()+ " "+ gnetwork.getNoOfLinks());
			
		}
	}
	

	
	/**
	 * This method simulates the spread of infection for a network that is read in, for the first time, and updates the 'status' of the immunized nodes.
	 *
	 */
	public SIRNetwork simulatespread_firstspread(SIRNetwork network)
	{
		Random random = new Random(System.currentTimeMillis());
		BufferedWriter out = null;
		BufferedWriter out2 = null;
		BufferedReader in = null;
	
	
		
		//Now simulate spread of infection.
		int startID = 830;  //Pick a node that is likely to be peripheral -  The higher the ID, the more peripheral node is.
		startID =  (int) (Math.abs(random.nextInt()) %network.getSize()); // find a start node at random
		 try {
		//	startID =  (int) (network.getHighestDegreeNodeID()); // Use the hub as starting node
		} catch (Exception e2) {

			e2.printStackTrace();
		}
		
		try {
			System.out.println("First infected node "+ startID+ " with degree "+ network.getNode(startID).getNumberOfLinks() );
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		
		final int TIMESTEPS = 31;
		final double beta = 0.4;  //makes prob = 1- exp(-1*0.1)= 0.09 for infection to take off from a peripheral node at each timestep. Comparable with Gamma.
		final double gamma = 0.1;

		
		try {
			((SIRNode)(network.getNode(startID))).setSIRState("I");
			
			 out = new BufferedWriter(new FileWriter(".\\..\\Output\\HI\\trial.xls"));
			
			 for(int i=0; i<TIMESTEPS;i++)
			 {

				 //simulate spread 
				 for(int j=0; j<network.getSize();j++)
				 {
					 SIRNode thisNode =  ( SIRNode)network.getAllNodes().elementAt(j);
					 
					 ////////>>>>>>>>>>>>>>>>
					 
					// Check the current status of the node 
					 // make sure that the SYNCHORIZED state is used for comparison, not the (partially updated) asynchronized state
					 
					 //If the node is susceptible, it can be infected with the probability 1 - e^(-bI), where I is the number of Infected nodes it is 
					 //connected to.
					 if( ((SIRNode)(thisNode)).getSynchronizedSIRState().equals("S"))
					 {
						 double prob = 1 - Math.exp(-beta*thisNode.noOfInfectedNeoghbours());
					
						 if(random.nextDouble() < prob)
						 {
							 thisNode.setSIRState("I");
						 }
					 }
					 
					 //If the node is Infected, it can become Recovered (removed) at a given timestep with probability gamma
					 if( ((SIRNode)(thisNode)).getSynchronizedSIRState().equals("I"))
					 {
						 if( random.nextDouble() <gamma)
						 {
							 thisNode.setSIRState("R");
						 }
					 }
					 
					 //If the node is recovered, its state never changes thereafter, as nodes are considered to have lifelong immunity
					
					 if( ((SIRNode)(thisNode)).getSynchronizedSIRState().equals("R"))
					 {
						 //DO NOTHING
					 }
					 

					 if( ((SIRNode)(thisNode)).getSynchronizedSIRState().equals("V"))
					 {
						 //DO NOTHING
					 }
					 
					 
				 }
				 
				 //Now synchronize all nde states.
				 
				 network.SynchronizeAllSIRStates();
				 
				 //THE SPREAD HAS BEEN COMPLETED FOR ONE TIMESTEP. PRINT RESULTS
				//System.out.println(i+" One round of spread completed: Infetced nodes "+ network.countInfectedNodes());
				//System.out.println(i+" One round of spread completed: Immunized nodes "+ network.countImmunizedNodes());
				
			
				 
			 }
			
			
			
			
			
		} catch (Exception e) {
			
			e.printStackTrace();
		}
		
		network.forceRecoveryofAllNodes(); //Make sure all "I" nodes become "R" nodes
		network.confirmVaccination(); //Now recovered nodes have state "V"
		System.out.println("Immunized nodes "+ network.countImmunizedNodes());
		System.out.println("recovered non immunized nodes "+ network.countRecoveredNonImmunizedNodes());
		System.out.println("Infected nodes now nodes "+ network.countInfectedNodes());
		System.out.println("Susceptible nodes now "+ network.countSusceptibleNodes());
		
		
		System.out.println("END method: simulatespread_firstspread()");
		return network;
	}
	
	
	
	/**
	 * This method simulates the spread of infection for a network that is parameter, for the SECOND time, and updates the 'status' of the immunized nodes.
	 * The network which is a parameter is therefore residual network
	 * But it still contains immunized nodes, with status R
	 * We calculate the fraction of the nodes which escape infection, at the end of this secondary epidemic.
	 * @return a vector containing the  vaccinated percentages and surviving
	 */
	public Vector simulatespread_secondspread(SIRNetwork network)
	{
		System.out.println("SECONDARY epidemic begins" );
		Random random = new Random(System.currentTimeMillis());
		BufferedWriter out = null;
		BufferedWriter out2 = null;
		BufferedReader in = null;
	
	
		double residualNetworkSize = network.getSize() - network.countImmunizedNodes();
		
		//Now simulate spread of infection.
		int startID = 4300;  //Pick a node that is likely to be peripheral -  The higher the ID, the more peripheral node is.
		int index = 0;
		Vector<SIRNode> susNodes = network.getAllSusceptibleNodes();
		index =  (int) (Math.abs(random.nextInt()) %susNodes.size()); // find a start node at random WHICH IS AT SUSCEPTIBLE STATE.
		startID = (int)susNodes.elementAt(index).getID(); //Using index of the vector, we find the ID of the node at that position, which might be higher than the size of vector
		 try {
		//	startID =  (int) (network.getHighestDegreeNodeID()); // Use the hub as starting node

		} catch (Exception e2) {

			e2.printStackTrace();
		}
		
		try {
			System.out.println("First infected node "+ startID+ " with degree "+ network.getNode(startID).getNumberOfLinks() );
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		
		final int TIMESTEPS = 200;
		final double beta = 0.4;  //makes prob = 1- exp(-1*0.1)= 0.09 for infection to take off from a peripheral node at each timestep. Comparable with Gamma.
		final double gamma = 0.1;

		
		try {
			((SIRNode)(network.getNode(startID))).setSIRState("I");
			
			 out = new BufferedWriter(new FileWriter(".\\..\\Output\\HI\\trial.xls"));
			
			 for(int i=0; i<TIMESTEPS;i++)
			 {
					
				 //simulate spread 
				 for(int j=0; j<network.getSize();j++)
				 {
					 SIRNode thisNode =  ( SIRNode)network.getAllNodes().elementAt(j);
					 
					 ////////>>>>>>>>>>>>>>>>
					 
					// Check the current status of the node 
					 // make sure that the SYNCHORIZED state is used for comparison, not the (partially updated) asynchronized state
					 
					 //If the node is susceptible, it can be infected with the probability 1 - e^(-bI), where I is the number of Infected nodes it is 
					 //connected to.
					 if( ((SIRNode)(thisNode)).getSynchronizedSIRState().equals("S"))
					 {
						 double prob = 1 - Math.exp(-beta*thisNode.noOfInfectedNeoghbours());
					
						 if(random.nextDouble() < prob)
						 {
							 thisNode.setSIRState("I");
						 }
					 }
					 
					 //If the node is Infected, it can become Recovered (removed) at a given timestep with probability gamma
					 if( ((SIRNode)(thisNode)).getSynchronizedSIRState().equals("I"))
					 {
						 if( random.nextDouble() <gamma)
						 {
							 thisNode.setSIRState("R");
						 }
					 }
					 
					 //If the node is recovered, its state never changes thereafter, as nodes are considered to have lifelong immunity
					
					 if( ((SIRNode)(thisNode)).getSynchronizedSIRState().equals("R"))
					 {
						 //DO NOTHING
					 }
					 
					 if( ((SIRNode)(thisNode)).getSynchronizedSIRState().equals("V"))
					 {
						String ss =  ((SIRNode)(thisNode)).getSIRState();
						//System.out.println("ss "+ss);
					 }
					
				 }
				 
				 //Now synchronize all nde states.
				 
				 network.SynchronizeAllSIRStates();
				 
				 //THE SPREAD HAS BEEN COMPLETED FOR ONE TIMESTEP. PRINT RESULTS
				//System.out.println(i+" One round of spread completed: Infetced nodes "+ network.countInfectedNodes());
				//System.out.println(i+" One round of spread completed: Immunized nodes "+ network.countImmunizedNodes());
				
			
				 
			 }
			
			
			
			
			
		} catch (Exception e) {
			
			e.printStackTrace();
		}
		
		//calculate the number of susteptible nodes now
		double immunizedPercentage = 100*(network.getSize() -  residualNetworkSize)/ network.getSize();
		double sus =  network.countSusceptibleNodes();
		double vac = network.getSize() -  residualNetworkSize; //Number of nodes already vaccinated
		//double survivingPercentage = 100*sus / residualNetworkSize;
		double survivingPercentage = 100*(sus+vac) / network.getSize(); //The percentage of nodes which are never infetced, including the vaccinated ones.
		
		System.out.println("Residual network size "+(residualNetworkSize));
		System.out.println("NEVER infected nodes now "+(sus));
		System.out.println("immunized percentage "+(immunizedPercentage));
		System.out.println("susviving percentage "+(survivingPercentage));
		System.out.println("END method: simulatespread_secondspread()");
		
		//At this point, immununized nodes are already lost -  how?

		Vector v =  new Vector();

		Percentages perc =  new Percentages();
		perc.immunizedPercentage = immunizedPercentage;
		perc.survivingPercentage = survivingPercentage;
		v.addElement(perc);
		return v;
	}
	
	
	/**
	 * This method immunizes a certain propotion of the network based on BC  ordering
	 * The immunized nodes are given R state
	 * All other nodes will still be in S state
	 * @param network
	 */
	public SIRNetwork immunizeNetwork_BC(SIRNetwork network, double vp)
	{
		//first, calculate BC for all nodes
		 if(vp==1)
		 {
			 network.populate_BC();
			 System.out.println("BC populated");
		 }
		 else
		 {
			 //should be already populated
		 }
		 
		final double VAC_PERCENT = vp; //In percentage
		int count = 0;
		network.sortNodes(); //Nodes have to have BC as comparing quantity
		
		 int threshold_node_id = (int)( VAC_PERCENT*0.01*network.getSize()); //tHE '-1' is because node IDs start from 0
		 double nodesToBeVaccinated = VAC_PERCENT*0.01*network.getSize();
		
		 //calculate the BC threshold above which nodes need to be vaccinated
		 double bc_threshold = network.getAllNodes().elementAt(threshold_node_id).getBCCentrality();
		System.out.println("Threshodl BC = "+ bc_threshold );
		 
		 for(int yy=0; yy<network.getSize();yy++)
		 {
			 SIRNode thisnode =  (SIRNode)network.getAllNodes().elementAt(yy);
			 if(thisnode.getBCCentrality() > bc_threshold)
			 {
				 thisnode.setSIRState("R");
				 count++;
			 }
		 }
		 
		 if( (nodesToBeVaccinated  - count) >= 2)  //More nodes nees to be vaccinated until the count is satisfied
		 {
			
			 for(int yy=0; yy<network.getSize();yy++)
			 {
				 if(count>= nodesToBeVaccinated)
				 {
					 break;
				 }
				 SIRNode thisnode =  (SIRNode)network.getAllNodes().elementAt(yy);
				 if(thisnode.getBCCentrality() == bc_threshold)
				 {
					 thisnode.setSIRState("R");
					 count++;
				 }
				 
		
			 }
		 }
		 
		 System.out.println("immununized nodes="+ count);
		 ((SIRNetwork)network).SynchronizeAllSIRStates();
		 
		 return  ((SIRNetwork)network);
	}
	
	
	/**
	 * This method immunizes a certain propotion of the network by choosing nodes randomly
	 * The immunized nodes are given R state
	 * All other nodes will still be in S state
	 * @param network
	 */
	public SIRNetwork immunizeNetwork_Random(SIRNetwork network, double vp)
	{
		Random random = new Random(System.currentTimeMillis());
		final double VAC_PERCENT = vp; //In percentage
		int count = 0;
	    final double rand_threshold = 0.01*VAC_PERCENT;  //Making the percentage a fraction
		
	

		 
		 for(int yy=0; yy<network.getSize();yy++)
		 {
			 if(random.nextDouble()<= rand_threshold)
			 {
				 SIRNode thisnode =  (SIRNode)network.getAllNodes().elementAt(yy);
				 thisnode.setSIRState("R");
				 count++;
			 
			 }
		 }
		 
		 System.out.println("immununized nodes="+ count);
		 ((SIRNetwork)network).SynchronizeAllSIRStates();
		 
		 return  ((SIRNetwork)network);
	}
	
	/**
	 * A private class to structure the percentages returned from simulatespread
	 * 
	 */
	
	private class Percentages
	{
		public double immunizedPercentage;
		public double survivingPercentage;
	}
	
}//end class
