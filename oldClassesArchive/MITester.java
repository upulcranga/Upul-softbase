package oldClassesArchive;

import java.util.Random;
import java.util.Vector;


import Basic.Network;
import Basic.Node;
import Distributions.QDistribution_PowerLaw;

public class MITester {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	
		/* STAR NETWORK*/
		
		/* Create a star network */
		
		Network starNet = new Network();
		
		/* Star graph has 17 nodes */
		Node centreNode = new Node(1,16);
		starNet.addNode(centreNode);
		
		for (int i = 2; i <= 17; i++) 
		{
			Node node = new Node(i,1);
			starNet.addNode(node);
			try
			{
				starNet.addLink(centreNode, node);
			}
			catch(Exception e)
			{
				System.out.println(e.getMessage());
				System.out.println("Could not add link in star graph");
			}
		}
		
		/* Now star graph built - print out its properties */
		
		AssortativenessCalculator ac = new AssortativenessCalculator(starNet);
		System.out.println("STAR graph r=" + ac.calculate_r());
		QDistribution_PowerLaw starq = new QDistribution_PowerLaw(starNet);
		System.out.println("STAR graph <k>=" + starq.getAverageDegree());
		
		System.out.println("STAR graph H(q)=" + starNet.calculateHq());
		System.out.println("STAR graph Hc(q)=" + starNet.calculateHqc());
		System.out.println("STAR graph I(q)=" + starNet.calculateMI());
		
		/* NOTE: STAR GRAPH GIVES EXPECTED RESULTS NOW */
		
		/* ERDOS RENYI GRAPH */
		int ER_NODES = 300;
		int ER_LINKS = (int)(300*0.5*6.82);
		
		Network ERNet = new Network();
		
		Random rand = new Random();
		
		/*Generate Nodes */
		for (int i = 1; i <= ER_NODES; i++) 
		{
			/*Follow Uniform Degree Distribution:  2-12 degree */
			Node node = new Node(i,2+ (Math.abs(rand.nextInt())%11) );
			ERNet.addNode(node);
			
		}
		
		/*Link Nodes randomly */
		 
		Vector<Node>ernodes = ERNet.getAllNodes();
		for (int i = 1; i <= ER_LINKS; i++) 
		{
			Node node1 = ernodes.elementAt((Math.abs(rand.nextInt())%ernodes.size()) );
			Node node2 = ernodes.elementAt((Math.abs(rand.nextInt())%ernodes.size()) );
			try
			{
				if( (node1.getNumberOfLinks() <=11) && (node2.getNumberOfLinks() <=11)) /*We can add MAX 12 links */
				{
					ERNet.addLink(node1, node2);
				}
				else
				{
					ER_LINKS++;
				}
			}
			catch(Exception e)
			{
				System.out.println(e.getMessage());
				System.out.println("Could not add link in ER graph");
			}
		}
		
		/* Now ER graph built - print out its properties */
		
		ac = new AssortativenessCalculator(ERNet);
		System.out.println("ER graph r=" + ac.calculate_r());
		starq = new QDistribution_PowerLaw(ERNet);
		System.out.println("ER graph <k>=" + starq.getAverageDegree());
		
		System.out.println("ER graph H(q)=" + ERNet.calculateHq());
		System.out.println("ER graph Hc(q)=" + ERNet.calculateHqc());
		System.out.println("ER graph I(q)=" + ERNet.calculateMI());
		
		starq = new QDistribution_PowerLaw(2.8,2);

	}

}
