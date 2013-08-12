package Basic;

/**
 * BetweennessCentrality : Implementation in Java to calculate the
 * betweenness centrality index distribution of a given unweighted graph G=(V, E)
 * 
 * BetweennessCentrality.java
 *
 * Date : 07/20/2008
 */


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.List;
import java.util.Properties;
import java.util.Random;

/**
 * This is the a modified version of Namal Senaratne's code
 * The types of Graph(Network) and Node have been changed to match Piraveenan's code.
 * modified: 05/26/2009
 *
 * @author  Namal Senarathne
 * @organization    Department of Computer Science & Engineering, University of Moratuwa
 */
public class ComponentAnalyzer
{ 
    public static void main(String[] args) throws InterruptedException, IOException
    {  
    	Random random = new Random(System.currentTimeMillis());
       
    	BufferedWriter out = null;
    	try{
    	
    	out = new BufferedWriter(new FileWriter(".\\LA-Extended\\attack 2.xls"));
    	}
    	
    	 catch (IOException e) {
 			System.out.println("cant write to file");
 			e.printStackTrace();
 		}

       NetworkGenerator gen = new NetworkGenerator();
     // String inputFilePath = properties.getProperty("file_path");
       
       Network graph = new Network();
       graph= gen.generateFromFile();
       
      
      ComponentAnalyzer bc = new ComponentAnalyzer();
       

  
      
  	graph.sortNodes();
    int gSize = graph.getSize();
	int[][] paths =  new int[gSize][gSize];
	
	
	for(int i = 0; i < graph.getSize(); i++)
	{
       
	      try {
			graph.calculateShortestPath((int)graph.getAllNodes().elementAt(i).getID());
			paths =  graph.returnPathArrays((int)graph.getAllNodes().elementAt(i).getID(), paths);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	     System.out.println("i="+i);
	     
	}  
	
	

	
	
	boolean x =  graph.isShortestPath(paths, 0,1);
	boolean y =  graph.isShortestPath(paths, 1,0);

       
     
       
       System.out.println("Size of graph AT START = " + graph.getSize());
       
       int BCSize = graph.getBiggestComponentSize_sp(paths);
       System.out.println("Biggest Component of graph AT START= " + BCSize);
    
      System.out.println("Biggest Component of graph = " + BCSize);
      System.out.println("Biggest Component percentage of graph = " + ((double)BCSize / (double)graph.getSize()));
       out.write("0"+"\t"+BCSize);
       out.flush();
       out.write( System.getProperty(("line.separator")));
       out.flush();
       

       //Simulate targeted attacks
       int initialSize = graph.getSize();
       for(int i=0; i<initialSize; i++ )
       {
    	   try {
    		int id = 0;
    		id = graph.getHighestDegreeNodeID();
    		//id = (int)Math.floor(random.nextDouble()*graph.getSize());
    		System.out.println("node to be removed degree = " + graph.getNode(id).getNumberOfLinks());
    		if(graph.getSize() > 1)
    		{
    			if(graph.getNode(id).getNumberOfLinks() <= 1)
    			{
    				int warning = 0;
    			}
				graph.removeNode(id);
				System.out.println("Size of graph = " + graph.getSize());
				
				//reset shortest paths 
				gSize = graph.getSize(); //the size is ONE LESS now!
		
	
			
					for(int p = 0; p < graph.getSize(); p++)
					{
				       
					      try {
							graph.calculateShortestPath((int)graph.getAllNodes().elementAt(p).getID());
							paths =  graph.returnPathArrays((int)graph.getAllNodes().elementAt(p).getID(), paths);
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
				    
				  
				       if((p%100) ==0)
				       {
				    	   System.out.println("p="+p);
				       }
					}  
					
				BCSize = graph.getBiggestComponentSize_sp(paths);
			    System.out.println("Biggest Component of graph now= " + BCSize);
			    //System.out.println("Biggest Component percentage of graph now= " + ((double)graph.getBiggestComponentSize() / (double)graph.getSize()));
			    out.write((initialSize - graph.getSize())+""+"\t"+BCSize);
				out.flush();
				out.write( System.getProperty(("line.separator")));
				out.flush();
				if(BCSize== 251)
				{
					graph.createTextFile("ecoli 251");
				}
				if(BCSize== 349)
				{
					//graph.createTextFile("ecoli 349");
				}
    		}
    		
    	
		       
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
       }
       
  
       System.out.println("betweeness tester terminates");
       
    }

    /*
     * Calculates Betweenness centrality distribution according to Brande's algorithm
     */
    void getBrandesCentrality(Network graph) throws InterruptedException, IOException
    {
    	File outFile = new File("Centralityoutput.txt");
    	if(outFile.exists())
    		outFile.delete();
    	BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
    	
    	/* size is incremented by 1 so that node's can start from 1 rather
          than 0, this eases the operation */
    	int gSize = graph.getSize() + 1;  

    	 float div = (float)(gSize - 2)*(gSize - 3);
    	 float[] BCentrality;
		try {
			BCentrality = graph.calcBetweenessCentrality();
		
    	 for(int v = 0; v < gSize-1; v++)
    	 {
           System.out.println(v + "\t" + (BCentrality[v]/div));
           writer.write(v + "\t" + (BCentrality[v]/div) + "\n");
    	 }
       
    	 writer.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

    /**
     * Another implementation of BC centrality by Namal Senaratne
     * @param graph
     * @throws InterruptedException
     * @throws IOException
     */
    void getNamalCentrality(Network graph) throws InterruptedException, IOException, Exception
    {
    	File outFile = new File("output.txt");
    	if(outFile.exists())
    		outFile.delete();
    	BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
    	
    	/* size is incremented by 1 so that node's can start from 1 rather
          than 0, this eases the operation */
    	int gSize = graph.getSize();
      
    	float[] BCDistribution = new float[gSize];
       
    	int[][] paths =  new int[gSize][gSize];
    	for(int i = 0; i < graph.getSize(); i++)
    	{
           
          graph.calculateShortestPath(i);
           graph.updateAllPairSPArrays(i, BCDistribution);
           paths =  graph.returnPathArrays(i, paths);
           
           System.out.println("i="+i);
    	}  
    	
    	
    	
    	boolean x =  graph.isShortestPath(paths, 0,1);
    	boolean y =  graph.isShortestPath(paths, 1,0);
      
       
    	float div = (float)(gSize - 2)*(gSize - 3); 
    	for(int v = 0; v < gSize; v++)
    	{
    		System.out.println(v + "\t" + (BCDistribution[v]/div));
    		writer.write(v + "\t" + (BCDistribution[v]/div) + "\n");
    	}
       
    	writer.close();
    }
    
}

