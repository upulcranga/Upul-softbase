/**
 * @author Piraveenan
 */
package Basic;

import Distributions.QDistribution_PowerLaw;

public class Link implements Comparable {
	
	private Node iNode; /*Node at one end of this link */
	private Node jNode; /* Node at the other end of this Link */
	private boolean MARKED;  /* means this Link is marked for deletion  -  needed because deleting links straightaway can lead to bugs */
	
	
	public Link(Node node1, Node node2)
	{
		iNode = node1;
		jNode = node2;
		MARKED = false;
	}
	
	public int iRemainingDegree()
	{
		return iNode.getAssignedRemainingDegree();
	}
	
	public int jRemainingDegree()
	{
		return jNode.getAssignedRemainingDegree();
	}
	
	public int iRealRemainingDegree()
	{
		return (iNode.getNumberOfLinks()-1);
	}
	
	public int jRealRemainingDegree()
	{
		return (jNode.getNumberOfLinks()-1);
	}
	
	public int iDegree()
	{
		return iNode.getNumberOfLinks();
	}
	
	public int jDegree()
	{
		return jNode.getNumberOfLinks();
	}
	
	public long iNodeID()
	{
		return iNode.getID();
	}
	
	public Node jNode()
	{
		return jNode;
	}
	
	public Node iNode()
	{
		return iNode;
	}
	
	public long jNodeID()
	{
		return jNode.getID();
	}
	
	public int iNodeLinks()
	{
		return iNode.getNumberOfLinks();
	}
	
	public int jNodeLinks()
	{
		return jNode.getNumberOfLinks();
	}
	
	public void markForDeletion()
	{
		MARKED = true;
	}
	
	public boolean isMarked()
	{
		return MARKED;
	}
	
	/*
	 * Note: this method compares two Links based on the sum of their degrees on each end.
	 * This implementation needs to be changed if any other way of comparing becomes necessary
	 */
	
	public int compareTo(Object ol)
	{
		Link otherLink = (Link)ol;
		int sumDegrees = iNodeLinks() + jNodeLinks();
		int otherSumDegrees = otherLink.iNodeLinks() + otherLink.jNodeLinks();
		if(sumDegrees < otherSumDegrees)
		{
			return 1;
		}
		return -1;
	}
	
	/*
	 * Note: this method compares two Links based on the DIFFERENCE of their degrees on each end.
	 * That is, the more assortative links go first.
	 * This implementation needs to be changed if any other way of comparing becomes necessary
	 */
	public int compareToDIFF(Object ol)
	{
		Link otherLink = (Link)ol;
		int difDegrees = Math.abs(iNodeLinks() -  jNodeLinks());
		int otherDiffDegrees = Math.abs(otherLink.iNodeLinks() - otherLink.jNodeLinks());
		if(difDegrees < otherDiffDegrees)
		{
			return 1;
		}
		else if(difDegrees > otherDiffDegrees)
		{
			return -1;
		}
		else //equal difference
		{
			if(iNodeLinks() > otherLink.iNodeLinks() )
			{
				return 1;
			}
			return -1;
		}
	}
	
	
	/*
	 * Note: this method compares two Links based on the contribution of assortativeness they make
	 * That is, the more assortative links go first.
	 * A links assortativeness is calculated as the SUM of local assortativeness contributions of each Node of taht link
	 * Local assortativeness of a node is calculated omitting the scaling terms, since scaling terms can be calculated only if know of the network
	 * and scaling terms are nnot needed for comparison.
	 * 
	 */
	public int compareToEJK(Object ol)
	{
		Link otherLink = (Link)ol;
		double sumLocalr = (calculate_localEJK(iNode())  +  calculate_localEJK(jNode()) );
		double otherSumLocalr = (calculate_localEJK(otherLink.iNode()) + calculate_localEJK(otherLink.jNode()));
		if(sumLocalr > otherSumLocalr)
		{
			return 1;
		}
		else 
		{
			return -1;
		}
		
	}
	
	public double calculate_localEJK( Node node)
	{
		if(node.getNumberOfLinks()==0)
		{
			//System.out.println("ERROR: Linkless node in network");
			return 0;
		}
	
		//Now calculate the average remaining degree of all neighbours!
		
	
		int sumRemainingDegrees = 0;
		for(int y=0; y<node.getLinks().size();y++)
		{
			if(node.getLinks().elementAt(y).iNodeID() != node.getID())
			{
				//This is a neghbour
				sumRemainingDegrees = sumRemainingDegrees + node.getLinks().elementAt(y).iRealRemainingDegree();
			}
			
			if(node.getLinks().elementAt(y).jNodeID() != node.getID())
			{
				//This is a neighbour
				sumRemainingDegrees = sumRemainingDegrees + node.getLinks().elementAt(y).jRealRemainingDegree();
			}
		}
		
		double avgRD = ((double)sumRemainingDegrees )/ ((double)node.getNumberOfLinks());
		
	
		double j = node.getNumberOfLinks()-1;
		double ejkNODE = (j*(j+1)*avgRD)/(2.0);
		
		return ejkNODE;
	}
	
	

}
