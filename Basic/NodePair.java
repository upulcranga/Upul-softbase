package Basic;

/**
 * This is a generic class, like a C Struct, which constitures two instances of Node.
 * This represents the concept of a 'pair' of nodes
 * Note that this is a 'transparent' struct, both constituent nodes are public
 * Just works as a construct to hold tow nodes together
 * @author piraveenan
 *
 */
public class NodePair {
	
	public Node xNode;
	public Node yNode;
	
	public NodePair(Node x, Node y)
	{
		xNode = x;
		yNode = y;
	}

}
