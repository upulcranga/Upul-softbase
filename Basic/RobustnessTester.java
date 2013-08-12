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
 * This programme simulates betweenness and degree centrality based attacks
 * and measures the Robustness of networks
 * 
 * Modified version of 'BetweenessTester.java' which was itself a modified version of Namal Senaratnes code
 * @date 29.03.2012
 */


public class RobustnessTester
{ 
	public static void main(String[] args) throws InterruptedException, IOException
	{ 
		

		RobustnessTester sim =  new RobustnessTester();
		
		sim.simulate(1);
		
	}

  
    
    public  void simulate(int callArg) throws InterruptedException, IOException
    { 
    	
    	
    	System.out.println("start:"+callArg);
    	Random random = new Random(System.currentTimeMillis());

    	String FDFileName = "Copy of Random ER 10-store"; //Sometimes the name might be in the format of "copy of xyz"
    	String attacktype= "BC";
    	BufferedWriter out = null;
    	try{
    	
    	out = new BufferedWriter(new FileWriter(".\\..\\Output\\ROBUSTNESS\\" + FDFileName+" "+attacktype + " .xls"));
    	}
    	
    	 catch (IOException e) {
 			System.out.println("cant write to file");
 			e.printStackTrace();
 		}

       NetworkGenerator gen = new NetworkGenerator();

       
       Network graph = new Network();
       graph= gen.generateFromFile(".\\..\\Data\\attacknets\\"+FDFileName);
       //graph = gen.growNetworkPA(500, 1.2);
       //graph = gen.growNetworkPA(500, 1.2);
       //graph =  gen.generateERRandom(500, 1200);
       //graph.createTextFile(".\\..\\Data\\attacknets\\ER1");
       graph.createTextFile(".\\..\\Data\\attacknets\\"+FDFileName);
      // graph.createTextFile(".\\..\\Data\\attacknets\\"+FDFileName+ "-store");
      graph= gen.generateFromFile(".\\..\\Data\\attacknets\\"+FDFileName); //Printing and re-reading to avoid ID problems

       System.out.println("network read");
       double N = graph.getSize();

       
       System.out.println("Size of graph AT START = " + graph.getSize());
       System.out.println("Size of graph LINKS = " + graph.getNoOfLinks());
    
  

  
     graph.populate_BC();
       //graph.populate_CC();
  	//graph.sortNodes(); //NB: the nodes themselves have the right 'compareto' mechanisms.
    int gSize = graph.getSize();
	boolean[][] paths =  new boolean[gSize][gSize];
	try {
		paths =  graph.returnIndirectAdjacencyArray(paths);
	} catch (Exception e1) {
		// TODO Auto-generated catch block
		e1.printStackTrace();
	}
	

	// for(int i=0; i<graph.getSize(); i++ )
     //{
		 //System.out.println("CC node= " + graph.getAllNodes().elementAt(i).getID()+ " is "+graph.getAllNodes().elementAt(i).getCCCentrality());
	//	System.out.println(""+graph.getAllNodes().elementAt(i).getBCCentrality());
	       
    // }

	

	

       
	double mi = graph.calculate_MI_straight();

       
       System.out.println("Size of graph AT START = " + graph.getSize());
       System.out.println("Size of graph LINKS = " + graph.getNoOfLinks());
       
       int BCSize = graph.getBiggestComponentSize_sp(paths);
       System.out.println("Biggest Component of graph AT START= " + BCSize);
    
      System.out.println("Biggest Component of graph = " + BCSize);
      System.out.println("Biggest Component percentage of graph = " + ((double)BCSize / (double)graph.getSize()));
      out.write("Index \t net size \t  removednodes \t Larg.comp.size \t MI");
      out.flush();
      out.write( System.getProperty(("line.separator")));
      out.flush();
      
       out.write("0"+"\t"+graph.getSize()+"\t"+"0"+"\t"+BCSize+"\t"+mi);
       out.flush();
       out.write( System.getProperty(("line.separator")));
       out.flush();
       
 

       //Simulate targeted attacks
       int initialSize = graph.getSize();
       double robustnes_sum = 0;
       for(int i=0; i<initialSize; i++ )
       {
    	   try {
    		int id = 0;
    
    		id = graph.getHighestDegreeNodeID();
    	
    		//graph.populate_rho();
    		graph.populate_BC();
    		//graph.populate_CC();
    		//graph.sortNodes();
    
    		if(attacktype.equals("BC"))
    		{
    			System.out.println("Attacktype BC");
    			id = graph.getHighestBCNodeID();
    		}
    		
    		if(attacktype.equals("CC"))
    		{
    			System.out.println("Attacktype CC");
    			id = graph.getHighestCCNodeID();
    		}
    		
    		if(attacktype.equals("R"))
    		{
    			System.out.println("Attacktype Random");
    			id = Math.abs(random.nextInt())%graph.getSize();
    		}
    		
    		//id = (int) graph.getAllNodes().elementAt(0).getID();
 
    		try
    		{
    			//id = (int)Math.floor(random.nextDouble()*paths.length);
    		}
    	      
	   		 catch (Exception e) {
	   			System.out.println("random ID deleted already");
	   			continue;
	   		}
	   		 
    		System.out.println("node to be removed degree = " + graph.getNode(id).getNumberOfLinks());
    		//System.out.println("node to be removed bc = " + graph.getNode(id).getBCCentrality());
    		System.out.println("node to be removed cc = " + graph.getNode(id).getCCCentrality());
    		
    		if(graph.getSize() > 1)
    		{
    			if(graph.getNode(id).getNumberOfLinks() <= 1)
    			{
    				int warning = 0;
    			}
				graph.removeNode(id);
    		
				
		
				System.out.println("Size of graph = " + graph.getSize());
				if(graph.getSize()<=1)
				{
					 //adjustment for remaing indices where component size = 1
					for(int k = (i+1); k <N; k++)
					{
						mi = graph.calculate_MI_straight();
						robustnes_sum = robustnes_sum + k;
						BCSize = graph.getSize();
						
						out.write((k)+"\t"+ graph.getSize()+"\t"+(initialSize - graph.getSize())+"\t"+BCSize+"\t"+mi);
						out.flush();
						out.write( System.getProperty(("line.separator")));
						out.flush();
					}
				      
					out.write((initialSize )+"\t"+(graph.getSize() )+"\t"+(initialSize )+"\t"+0+"\t"+0);
					out.flush();
					out.write( System.getProperty(("line.separator")));
					out.flush();
					
				    double rb = (6*robustnes_sum) / (N*(N+1)*(N-1)); //CHECK FORMULA!
				    
				    System.out.println("NETWORK ROBUSTNESS="+100*rb+" percent");
				    out.write("\t"+"rb"+""+"\t"+100*rb+" percent");
					out.flush();
					
				}
				graph.createTextFile(".\\..\\Data\\attacknets\\"+FDFileName);
				graph= gen.generateFromFile(".\\..\\Data\\attacknets\\"+FDFileName);
				mi = graph.calculate_MI_straight();
				
				//reset shortest paths 
				gSize = graph.getSize(); //the size is ONE LESS now!
		
	
			
					    	  paths =  graph.returnIndirectAdjacencyArray(paths);
							//graph.calculateShortestPath((int)graph.getAllNodes().elementAt(p).getID());
							//paths =  graph.returnPathArrays((int)graph.getAllNodes().elementAt(p).getID(), paths);
						
					}  
					
				BCSize = graph.getBiggestComponentSize_sp(paths);
			    System.out.println("Biggest Component of graph now= " + BCSize);
			    //System.out.println("Biggest Component percentage of graph now= " + ((double)graph.getBiggestComponentSize() / (double)graph.getSize()));
			    out.write((i+1)+"\t"+ graph.getSize()+"\t"+ (initialSize - graph.getSize())+""+"\t"+BCSize+"\t"+mi);
				out.flush();
				out.write( System.getProperty(("line.separator")));
				out.flush();
				
				robustnes_sum = robustnes_sum + ((i+1)*BCSize);
				
    		}
    		
    	
		       
		 catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.println("NEW BLOCK: NODE not found  -  populating tail");
			 //adjustment for remaing indices where component size = 1
			for(int k = (i+1); k <N; k++)
			{
				robustnes_sum = robustnes_sum + k;
				BCSize = graph.getSize();
				
				out.write((k)+"\t"+ graph.getSize()+"\t"+(initialSize - graph.getSize())+"\t"+"1"+"\t"+mi);
				out.flush();
				out.write( System.getProperty(("line.separator")));
				out.flush();
			}
			System.out.println("NEW BLOCK -  BREAKING");
			break;
		}
       }
       
      
       
       double rb = (6*robustnes_sum) / (N*(N+1)*(N-1)); //CHECK FORMULA!
       System.out.println("NETWORK ROBUSTNESS="+100*rb+ " percent");
       out.write("\t"+"rb"+""+"\t"+100*rb+" percent");
       out.flush();
       System.out.println("betweeness tester terminates");
       
    }

    
   
      
       
    
}

