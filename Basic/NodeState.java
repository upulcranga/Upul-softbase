package Basic;

/**
 * BetweennessCentrality : Implementation in Java to calculate the
 * betweenness centrality index distribution of a given unweighted graph G=(V, E)
 * 
 * NodeState.java
 *
 * Date : 07/20/2008
 */

/**
 * The enumeration which specifies the three different types of states a 
 * node can be in
 *
 * @author  Namal Senarathne
 * @organization    Department of Computer Science & Engineering, University of Moratuwa
 */
public enum NodeState
{
    NOT_ENQUEUED,   /* The node is not yet put into the queue, so the distance is INF */
    QUEUED,         /* The node is already is in the queue */
    DEQUEUED        /* The node was removed earlier, thus has its shortest distance
                       is already calculated */
}
