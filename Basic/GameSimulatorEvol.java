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

public class GameSimulatorEvol

{

	/**
	 * @param args
	 * @throws IOException
	 * Written originally by Darshana around 14.04.2013
	 * Modified by Piraveen on 18.04.2013
	 * Integrated into Basic package of RC project by piraveen on 07.05.2013
	 */
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException
	{

		NetworkGenerator ng = new NetworkGenerator();
		final String root = "C:\\Documents and Settings\\piraveenan\\Desktop\\SOFTBASE1\\APASimulation\\Output\\GAME\\beta\\"; //Piraveens code to remember the output folder
		

		Network network = ng.generateFromFile("C:\\Documents and Settings\\piraveenan\\Desktop\\SOFTBASE1\\APASimulation\\Data\\hi\\Mod 1000 2000"); //P
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
			
			betaout2 = new BufferedWriter(new FileWriter(root+"evolution time Mod 1000 2000 beta 2-1.xls"));

		} 
		catch (IOException e1)
		{
			System.out.println("The file is not found");
			e1.printStackTrace();
		}


		double beta = 2.1;
		
		double alpha =  1.0;
		int lambda = 0;
			

				
				fromTime = System.currentTimeMillis();


				System.out.println(network.getSize() + " "
						+ network.getAllLinks().size()+" "+network.realMaxDegree());


				Random random = new Random(System.currentTimeMillis());
				
				double avg_CordCount = 0;
				double avg_NonCount = 0;
				
				double avg_CordVal = 0; //avg over multiple iterations
				double avg_NonVal = 0;
				
				//Declare variables to hold averages of averages over multiple runs
				
				int TIMESTEPS =  2500;
				double[] aa_CoCout = new double[TIMESTEPS];
				double[] aa_DefCout = new double[TIMESTEPS];
				double[] aa_CoDeg = new double[TIMESTEPS];
				double[] aa_DefDeg = new double[TIMESTEPS];
				double[] aa_MI = new double[TIMESTEPS];
				
				double INITS = 50;

				for (int q = 0; q <INITS; q++)
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


			
					
					for (int k = 0; k < TIMESTEPS; k++) //Timesteps
					{

						/****/

						if (k != 0) 
						{
							GameSimulatorEvol cg = new GameSimulatorEvol();

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

						//set up for calculating state-based network mutual information
						
						// also calculate average degrees
						
						double coAverageDegree = 0;
						double nonAverageDegree = 0;
						
						for (int t = 0; t < network.getSize(); t++)
						{
							Node n = network.nodes.elementAt(t);
							
							if(n.isCoordinator)
							{
								n.setReading(1.0);
								coAverageDegree =  coAverageDegree +  n.getNumberOfNeighbours();
							}
							else
							{
								n.setReading(0);
								nonAverageDegree =  nonAverageDegree + n.getNumberOfNeighbours();
							}
						}

						if(coordinatorCount!=0)
						{
							coAverageDegree =  coAverageDegree / (double)coordinatorCount;
						}
						else
						{
							coAverageDegree = 0;
						}
						if(defectorCount!=0)
						{
							nonAverageDegree = nonAverageDegree / (double)defectorCount;
						}
						else
						{
							nonAverageDegree = 0;
						}
						
						double MIstatbesed = 0;
					
						//Print the trend of evolution in some cases
						//if(q==0) //randomly choose a case
						{
							
							MIstatbesed = network.calculate_MI_straight_sensor();
							
							System.out.println("EVOL: "+defectorCount + "\t" + coordinatorCount+"\t" +MIstatbesed+"\t"+coAverageDegree + "\t" + nonAverageDegree+"\n");
							//betaout2.write(k+"\t"+defectorCount + "\t" + coordinatorCount+"\t"+MIstatbesed+"\t"+coAverageDegree + "\t" + nonAverageDegree);
							//betaout2.flush();
							//betaout2.write(System.getProperty(("line.separator")));
							//betaout2.flush();
						}
						
						aa_CoCout[k] = aa_CoCout[k] + (1.0/INITS)*coordinatorCount; //The division is to average over q
						aa_DefCout[k] = aa_DefCout[k] + (1.0/INITS)*defectorCount;
						aa_CoDeg[k] = aa_CoDeg[k] + (1.0/INITS)*coAverageDegree;
						aa_DefDeg[k] = aa_DefDeg[k] + (1.0/INITS)*nonAverageDegree;
						aa_MI[k] = aa_MI[k] + (1.0/INITS)*MIstatbesed;

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


					
					avg_CordCount =  avg_CordCount +  coordinatorCount;
					avg_NonCount = avg_NonCount + defectorCount;
					
					avg_CordVal =  avg_CordVal +  totalCoordinationValue;
					avg_NonVal = avg_NonVal + totalDefectorValue;
		




				} //End q (multiple initializations in a network) for-loop
				
				//Print Array
				
				for(int arr = 0; arr < TIMESTEPS; arr++)
				{
					
					
					
					System.out.println("ARR: "+"\t"+aa_DefCout[arr] + "\t" + aa_CoCout[arr]+"\t"+aa_MI[arr]+"\t"+aa_CoDeg[arr] + "\t" + aa_DefDeg[arr]+"\n");
					betaout2.write(arr+"\t"+aa_DefCout[arr] + "\t" + aa_CoCout[arr]+"\t"+aa_MI[arr]+"\t"+aa_CoDeg[arr] + "\t" + aa_DefDeg[arr]);
					betaout2.flush();
					betaout2.write(System.getProperty(("line.separator")));
					betaout2.flush();
				}
				
	
				
				avg_CordCount =  avg_CordCount / 200; //average over multiple ijitializations
				avg_NonCount =  avg_NonCount /200;
				
				avg_CordVal =  avg_CordVal / (200*avg_CordCount); //average over multiple initializations and all nodes, but cumulative payoff
				avg_NonVal =  avg_NonVal /(200*avg_NonCount);
				
			

				toTime = System.currentTimeMillis();	
				System.out.println("time take is : " + (toTime - fromTime));

			
			

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
