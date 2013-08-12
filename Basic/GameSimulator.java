package Basic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.Vector;

import Basic.Link;
import Basic.Network;
import Basic.NetworkGenerator;
import Basic.Node;

public class GameSimulator

{

	/**
	 * @param args
	 * @throws IOException
	 * Written originally by Darshana around 14.04.2013
	 * Modified by Piraveen on 18.04.2013
	 * Integrated into Basic package of RC project by piraveen on 07.05.2013
	 */
	public static void main(String[] args) throws IOException
	{

		NetworkGenerator ng = new NetworkGenerator();
		final String root = "C:\\Documents and Settings\\piraveenan\\Desktop\\SOFTBASE1\\APASimulation\\Output\\GAME\\beta\\"; //Piraveens code to remember the output folder
		

		Network network = ng.generateFromFile("C:\\Documents and Settings\\piraveenan\\Desktop\\SOFTBASE1\\APASimulation\\Data\\hi\\ER 1000 2000"); //P
		Node node = null;
		Vector<Link> links = null;
		Link link = null;

		int coordinatorCount = 0;
		int defectorCount = 0;


		long fromTime;
		long toTime;

		double totalCoordinationValue = 0;
		double totalDefectorValue = 0;


		BufferedWriter betaout = null;
		BufferedWriter betaout2 = null;



		try
		{
		
			betaout = new BufferedWriter(new FileWriter(root+"beta SW 200 NOPOINT.xls"));
			
			betaout2 = new BufferedWriter(new FileWriter(root+"evolution time ER 1000 2000 beta 2-1.xls"));

		} 
		catch (IOException e1)
		{
			System.out.println("The file is not found");
			e1.printStackTrace();
		}


		double beta = 2.1;
		//while(beta <3.0)
		//for (double beta = 1.0; beta < 4.00; beta = beta + 0.2)
		//for(double alpha = 0.0; alpha<=1.0; alpha = alpha +0.1)
		//for(int lambda = 0; lambda<2000; lambda=lambda+50)
		//{
			if( (beta >=1.6) && (beta <2.2))
			{
				beta =  beta + 0.05;
			}
			else
			{
				beta =  beta + 0.2;
			}
			//double beta = 2.1;
			double alpha =  1.0;
			int lambda = 0;
			
			for (int i = 0; i < 1; i++)
			
			{

				//network = ng.generateERRandom(1000, 8000);
				//network.createTextFile("ER 1000 8000");
				
				fromTime = System.currentTimeMillis();


				System.out.println(network.getSize() + " "
						+ network.getAllLinks().size()+" "+network.realMaxDegree());

				
				//System.exit(0);
				Random random = new Random(System.currentTimeMillis());
				
				double avg_CordCount = 0;
				double avg_NonCount = 0;
				
				double avg_CordVal = 0; //avg over multiple iterations
				double avg_NonVal = 0;

				for (int q = 0; q <2; q++)
				{


					coordinatorCount = 0;
					defectorCount = 0;

					//Randomly initializing the initial state
					for (int u = 0; u < network.getSize(); u++)
					{
						network.nodes.get(u).isCoordinator = false;
						network.nodes.get(u).gameTotalPayoff = 0;
						
					}
					
					for (int j = 0; j < network.getSize(); j++) 
					{
						if (random.nextDouble() < 0.5)
						{
							network.nodes.get(j).isCoordinator = true;
							coordinatorCount++;
						}
					}

					Node iNode, jNode = null;

					defectorCount = network.getSize() - coordinatorCount;  

					Node node3 = null;


					//To count node degrees
					//for(int m=0; m<network.getSize(); m++){
					//	node3 = network.nodes.get(m);
					//	betaout.write(node3.getNumberOfLinks() + "\t" + node3.isCoordinator);
					//	betaout.flush();
						
					//	betaout.write(System.getProperty(("line.separator")));
					//	betaout.flush();
						
					//}
					
					for (int k = 0; k < 1500; k++)
					{

						/****/

						if (k != 0) 
						{
							GameSimulator cg = new GameSimulator();

							double prob = 0;

							// for each link, run the update strategy
							for (int x = 0; x < links.size(); x++) 
							{
								link = links.get(x);

								iNode = link.iNode();
								jNode = link.jNode(); //We do not choose a random neighbour!
								
								
								//PIRAVEENS CODE
								double kmax =  iNode.getNumberOfNeighbours();
								
								
								double probab = (1.0 / kmax);
								if(random.nextDouble() > probab)
								{
									continue;  //Stochastic way of ensuring only one neighbour is selected.
								}
								
								//END PIRAVEENS CODE

								prob = random.nextDouble();

								if (iNode.isCoordinator != jNode.isCoordinator)
								{
									if (prob < cg.getUpdateStrategyProbability(iNode, jNode, beta, alpha, lambda))
									{
										iNode.isCoordinator = jNode.isCoordinator;

										if (jNode.isCoordinator)
										{
											coordinatorCount++;
											defectorCount--;
										} 
										else 
										{
											defectorCount++;
											coordinatorCount--;
										}
										
									}
									else
									{

										prob = random.nextDouble();
										if (prob < cg
												.getUpdateStrategyProbability(
														jNode, iNode, beta, alpha, lambda))
										{
											jNode.isCoordinator = iNode.isCoordinator;

											if (iNode.isCoordinator)
											{
												coordinatorCount++;
												defectorCount--;
											}
											else
											{
												defectorCount++;
												coordinatorCount--;
											}
										}
									}
								}

							}
						} //End if k!=0

						if (defectorCount < 0 || coordinatorCount < 0)
						{
							System.out.println();
						}
						/****/

				

						links = network.getAllLinks();
			
						for (int t = 0; t < network.getSize(); t++)
						{
							network.nodes.get(t).gameTotalPayoff = 0; //CumulativePayoff gets reassigned for every evolutionary step.
						}

						for (int l = 0; l < links.size(); l++)
						{
							link = links.get(l);

							iNode = link.iNode();
							jNode = link.jNode();

							if (iNode.isCoordinator)
							{
								if (jNode.isCoordinator)
								{
									iNode.gameTotalPayoff += beta;
									jNode.gameTotalPayoff += beta;
								}
								else
								{
									jNode.gameTotalPayoff += 1;
								}
							}
							else
							{
								if (jNode.isCoordinator)
								{
									iNode.gameTotalPayoff += 1;
								}
								else
								{
									iNode.gameTotalPayoff += 1;
									jNode.gameTotalPayoff += 1;
								}
							}
						}
						
						//Payoffs have been upgraded.  Now upgrade the buffer
						for (int n = 0; n < network.getSize(); n++)
						{

							node = network.nodes.get(n);
							node.updatePayoffBuffer(node.gameTotalPayoff);
						}

						
					
						//Print the trend of evolution in some cases
						if(q==0) //randomly choose a case
						{
							
						
							System.out.println("EVOL: "+defectorCount + "\t" + coordinatorCount+"\n");
							betaout2.write(k+"\t"+defectorCount + "\t" + coordinatorCount);
							betaout2.flush();
							betaout2.write(System.getProperty(("line.separator")));
							betaout2.flush();
						}

					} //End k- (multiple timesteps) for-loop
					

					//Calculate the total cumulative pay-offs of coordinators and non-coordinators at this stage.
					totalCoordinationValue = 0;
					totalDefectorValue = 0;

					for (int n = 0; n < network.getSize(); n++)
					{

						node = network.nodes.get(n);
					

						if (node.isCoordinator)
						{
							totalCoordinationValue += node.gameTotalPayoff;
						}
						else
						{
							totalDefectorValue += node.gameTotalPayoff;
						}

					}

					
					//System.out.println("end of evolution " +"beta "+ beta+" "+coordinatorCount + " " + defectorCount+ "\n");

					//network.createTextFile(root+"" + nlr
							//+ "\\" + String.valueOf(beta).substring(0, 3)
							//+ "\\" + i + "\\network");
					
					avg_CordCount =  avg_CordCount +  coordinatorCount;
					avg_NonCount = avg_NonCount + defectorCount;
					
					avg_CordVal =  avg_CordVal +  totalCoordinationValue;
					avg_NonVal = avg_NonVal + totalDefectorValue;
		


					//betaout.write(coordinatorCount + " " + defectorCount+ "\n");
				
					//betaout.flush();
				

				} //End q (multiple initializations in a network) for-loop
				
				avg_CordCount =  avg_CordCount / 200; //average over multiple ijitializations
				avg_NonCount =  avg_NonCount /200;
				
				avg_CordVal =  avg_CordVal / (200*avg_CordCount); //average over multiple initializations and all nodes, but cumulative payoff
				avg_NonVal =  avg_NonVal /(200*avg_NonCount);
				
				System.out.println("end of evolution " +"beta "+ beta+" "+avg_CordCount + " " + avg_NonCount+" "+avg_CordVal + " " + avg_NonVal+ "\n");
				betaout.write(beta+"\t"+avg_CordCount + "\t" + avg_NonCount+"\t"+avg_CordVal + "\t" + avg_NonVal);
				betaout.flush();
				betaout.write(System.getProperty(("line.separator")));
				betaout.flush();
				
				//System.out.println("end of evolution " +"alpha "+alpha+" "+avg_CordCount + " " + avg_NonCount+ "\n");
				//betaout.write(alpha+"\t"+avg_CordCount + "\t" + avg_NonCount);
				//betaout.flush();
				//betaout.write(System.getProperty(("line.separator")));
				//betaout.flush();
				
			
				//System.out.println("end of evolution " +"lambda "+ lambda+" "+avg_CordCount + " " + avg_NonCount+" "+avg_CordVal + " " + avg_NonVal+ "\n");
				//betaout.write(lambda+"\t"+avg_CordCount + "\t" + avg_NonCount+"\t"+avg_CordVal + "\t" + avg_NonVal);
				//betaout.flush();
				//betaout.write(System.getProperty(("line.separator")));
				//betaout.flush();
		  
				//For degree count
				//Node node3 =  null;
				//	for(int m=0; m<network.getSize(); m++){
				//	node3 = network.nodes.get(m);
				//	betaout2.write(node3.getNumberOfLinks() + "\t" + node3.isCoordinator);
				//	betaout2.flush();
					
				//	betaout2.write(System.getProperty(("line.separator")));
				//	betaout2.flush();
					
					
				//}
				
				toTime = System.currentTimeMillis();	
				System.out.println("time take is : " + (toTime - fromTime));

			} //End i-for loop (multiple networks for-loop)
			

		//} //End alpha-beta For loop

		System.out.println("MAIN METHOD COMPLETES");


	} //End Main

	public double getUpdateStrategyProbability(Node iNode, Node jNode,
			double beta)
	{
		double updateProb = 0;

		double Py = jNode.gameTotalPayoff;
		double Px = iNode.gameTotalPayoff;
		double T, S;

		T = 1;
		S = 0;

		double k;

		if (iNode.getAllNeighbours().size() > jNode.getAllNeighbours().size()) {
			k = iNode.getAllNeighbours().size();
		} else {
			k = jNode.getAllNeighbours().size();
		}

		updateProb = (Py - Px) / (k * beta);

		if (updateProb < 0) {
			updateProb = 0;
		}

		return updateProb;
	}
	
	public double getUpdateStrategyProbability(Node iNode, Node jNode, double beta, double alpha)
	{
		double updateProb = 0;

		Random random = new Random(System.currentTimeMillis());
		double rand = random.nextDouble();
		
		double Py = jNode.gameTotalPayoff;
		double Px = iNode.gameTotalPayoff;
		double T,S;
		
		T = 1;
		S = 0;
		
		
		double k;
		
		if(iNode.getAllNeighbours().size() > jNode.getAllNeighbours().size()){
			k = iNode.getAllNeighbours().size();
		}else{
			k = jNode.getAllNeighbours().size();
		}
		
		
		updateProb = (Py-Px)/(k * beta);
		
		if(updateProb < 0){
			updateProb = 0;
		}
		
		updateProb = ((1 - alpha) * rand) +  alpha * updateProb;
		
		return updateProb;
	}
	
	/**
	 * This method uses an update strategy with a time lag from the payoff values of neighbours
	 * @param iNode
	 * @param jNode
	 * @param beta
	 * @param alpha
	 * @param lambda
	 * @return
	 */
	public double getUpdateStrategyProbability(Node iNode, Node jNode, double beta, double alpha, int lambda)
	{
		double updateProb = 0;

		Random random = new Random(System.currentTimeMillis());
		double rand = random.nextDouble();
		
		double Py = jNode.getPayoffBuffer()[lambda];
		double Px = iNode.getPayoffBuffer()[lambda];
		double T,S;
		
		T = 1;
		S = 0;
		
		
		//System.out.println("Px " +Px);
		//System.out.println("Py " +Py);
		double k;
		
		if(iNode.getAllNeighbours().size() > jNode.getAllNeighbours().size()){
			k = iNode.getAllNeighbours().size();
		}else{
			k = jNode.getAllNeighbours().size();
		}
		
		
		updateProb = (Py-Px)/(k * beta);
		
		if(updateProb < 0){
			updateProb = 0;
		}
		
		updateProb = ((1 - alpha) * rand) +  alpha * updateProb;
		
		return updateProb;
	}

} //End class
