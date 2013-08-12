/*
 * FileProcessor.java
 *
 * Created on January 9, 2008, 10:30 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package Basic;

import java.io.*;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Pattern;
import Basic.*;
/**
 *
 * @author kishan
 */
public class FileProcessor {
    
    File file = null;
    List nodes;
    
    int count = 0;
    Network net;
    
    public FileProcessor() {
        net = new Network();
        nodes = new ArrayList();
       
    }
    
    public void SetFile(String fname){
        file =  new File(fname);
    }
    
    public void Process(){
    	int lineCounter = 0;
        try {
            Scanner scanner = new Scanner(file);
            while ( scanner.hasNextLine() ){
            	
             	
            	//To measure performance of programme
            	//==========================================================
            	lineCounter ++;
            	//System.out.println("Parsing Line "+ lineCounter);
            	//==========================================================
            	
            	
            	if(lineCounter > 800000)
            	{
            		return;
            	}
            	
                processlines(scanner.nextLine());
            }
            scanner.close();
        } catch (IOException ex){
            
        }
        
    }
    
    public void Process_ForTimeNodes(){
    	int lineCounter = 0;
        try {
            Scanner scanner = new Scanner(file);
            while ( scanner.hasNextLine() ){
            	
             	
            	//To measure performance of programme
            	//==========================================================
            	lineCounter ++;
            	//System.out.println("Parsing Line "+ lineCounter);
            	//==========================================================
            	
            	
            	if(lineCounter > 800000)
            	{
            		return;
            	}
            	
                processlines_ForTimeNodes(scanner.nextLine());
            }
            scanner.close();
        } catch (IOException ex){
            
        }
        
    }
    
    public void ProcessVERBOSE(){
        try {
            Scanner scanner = new Scanner(file);
            int lineCounter = 0;
            int lineCounterk = 0;
            
            while ( scanner.hasNextLine() ){
            	lineCounter++;
            	if(lineCounter > 1000)
            	{
            		lineCounterk++;
            		lineCounter = 0;
            		System.out.println("nodes lines 1k " + lineCounterk);
            	}
           
                processlines(scanner.nextLine());
            }
            scanner.close();
        } catch (IOException ex){
            
        }
        
    }
    
    /**
     * Creates networks with 'TimeNodes', a sublass of nodes
     * @param line
     */
    public void processlines_ForTimeNodes(String line){
        Pattern p = Pattern.compile(" ");
        String[] items = p.split(line);
        //lmap[count][1] = items[0];
        //lmap[count][2] = items[1];
        int ii =0;
        for(int i = 0; i < items.length ; i++){
            if(ii > 1) break;
                    if(!items[i].equals("")){
        if(!nodes.contains(items[i]) ){
            nodes.add(count ,items[i] );
            TimeNode x = new TimeNode(count);
           //Node x = new Node(Integer.parseInt(items[i])); //INSERTING NODES WITH IDs given in the file
            //Rather than with new IDs, like Kishan's stupid code above does.
            net.addNode(x);
            count++;
            ii++;
            }
          }
        }
       
        
    }
    
    public void processlines(String line){
        Pattern p = Pattern.compile(" ");
        String[] items = p.split(line);
        //lmap[count][1] = items[0];
        //lmap[count][2] = items[1];
        int ii =0;
        for(int i = 0; i < items.length ; i++){
            if(ii > 1) break;
                    if(!items[i].equals("")){
        if(!nodes.contains(items[i]) ){
            nodes.add(count ,items[i] );
            Node x = new Node(count);
          // Node x = new Node(Integer.parseInt(items[i])); //INSERTING NODES WITH IDs given in the file
            //Rather than with new IDs, like Kishan's stupid code above does.
            net.addNode(x);
            count++;
            ii++;
            }
          }
        }
        /*
        for(int i = 0; i < items.length ; i++){
                    if(!items[i].isEmpty()){
        if(!nodes.contains(items[1]) ){
            nodes.add(count ,items[1] );
            Node x = new Node(count);
            net.addNode(x);
            count++;
        }
         }
        }
         */
        
    }
    
    public void  showall(){
        Object key = null, value;
        Iterator it = nodes.iterator();
        while (it.hasNext()) {
            
            key =  it.next();
            
            
        }
    }
    
    public int  getSize(){
        return nodes.size();
        
    }
    
    public List  getList(){
        return nodes;
        
    }
    
    public Network  getNetwork(){
        return net;
        
    }
    
    public int buildlinks(){
        int x,y;
        try {
            Scanner scanner = new Scanner(file);
           
            while ( scanner.hasNextLine() ){
            	
                Pattern p = Pattern.compile(" ");
                String[] items = p.split(scanner.nextLine());
                
                int ii = 0;
                String [] items1 = new String[items.length];
                        //0 3 5
                for(int i = 0; i < items.length ; i++){
                    if(!items[i].equals("")){
                       //System.out.println(i);
                       items1[ii] = items[i].trim();
                     /*                      
                       if(i == 0){
                            items1[ii] = items[i];
                       }
                       if(i > 1 ){
                            items1[ii] = items[i]; 
                       }
                      */
                    ii++;
                    
                    if(ii >1)
                        break;
                    }
                    
                   // System.out.println(items1.length);
                }
                
               //System.out.println(items1[0] + " " + items1[1]);
               x = nodes.indexOf(items1[0].trim());
               y = nodes.indexOf(items1[1].trim());
              //x = Integer.parseInt(items1[0].trim()) ; //Replacing Kishan's stupid code to match with file ID
             // y= Integer.parseInt(items1[1].trim()) ;
                //System.out.println(x + " " + y );
                try {
                	if(x!=y) //leave self links out: TODO Get rid of this condition
                	{
                		net.addLink(x,y);
                		
                		
                		
                		//Namal's code
                		//nodeList[nodeID1].addAdjacendyNode(nodeID2);
                		Node nodex = net.getNode(x);
                		nodex.addAdjacendyNode(y);
                		Node nodey = net.getNode(y);
                		nodey.addAdjacendyNode(x);
                		
                	}
                } catch (Exception ex) {
                ex.printStackTrace();
                return 0;
                }
                }
            scanner.close();
        } catch (IOException ex){
            return 0;
        }
        return 1;
        
        
      /*  
        for(int i=0;i <= lmap.length; i++){
            x = nodes.indexOf(lmap[i][0]);
            y = nodes.indexOf(lmap[i][1]);
            try {
                net.addLink(x,y);
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        
        */
    }
    
    public int buildlinksVERBOSE(){
        int x,y;
        try {
            Scanner scanner = new Scanner(file);
            int lineCounter = 0;
            int lineCounterk = 0;
            while ( scanner.hasNextLine() ){
            	lineCounter++;
            	if(lineCounter > 1000)
            	{
            		lineCounterk++;
            		lineCounter = 0;
            		System.out.println("lines 1k " + lineCounterk);
            	}
                Pattern p = Pattern.compile(" ");
                String[] items = p.split(scanner.nextLine());
                
                int ii = 0;
                String [] items1 = new String[items.length];
                        //0 3 5
                for(int i = 0; i < items.length ; i++){
                    if(!items[i].equals("")){
                       //System.out.println(i);
                       items1[ii] = items[i].trim();
                     /*                      
                       if(i == 0){
                            items1[ii] = items[i];
                       }
                       if(i > 1 ){
                            items1[ii] = items[i]; 
                       }
                      */
                    ii++;
                    
                    if(ii >1)
                        break;
                    }
                    
                   // System.out.println(items1.length);
                }
                
               //System.out.println(items1[0] + " " + items1[1]);
                x = nodes.indexOf(items1[0].trim());
                y = nodes.indexOf(items1[1].trim());
                //System.out.println(x + " " + y );
                try {
                	if(x!=y) //leav self links out: TODO Get rid of this condition
                	{
                		net.addLink(x,y);
                	}
                } catch (Exception ex) {
                ex.printStackTrace();
                return 0;
                }
                }
            scanner.close();
        } catch (IOException ex){
            return 0;
        }
        return 1;
        
        
      /*  
        for(int i=0;i <= lmap.length; i++){
            x = nodes.indexOf(lmap[i][0]);
            y = nodes.indexOf(lmap[i][1]);
            try {
                net.addLink(x,y);
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        
        */
    }
    
    /**
     * This method builds links only if an attribute of the link is beyond a THRESHOLD. otherwise the link is ignored. The attribute is given in the thrid column.
     * @param THRESHOLD
     * @return
     */
    
    public int buildlinks_conditional(double THRESHOLD){
        int x,y;
        int lineCounter = 0;
        double linkWeight;
        try {
            Scanner scanner = new Scanner(file);
           
            while ( scanner.hasNextLine() ){
            	
            	
            	
            	
            	//To measure performance of programme
            	//==========================================================
            	lineCounter ++;
            	System.out.println("Parsing again line "+ lineCounter);
            	//==========================================================
            	
            	

            	if(lineCounter > 800000)
            	{
            		return 1;
            	}
            	
            	
            	
            	
                Pattern p = Pattern.compile("\t");
                String[] items = p.split(scanner.nextLine());
                
                int ii = 0;
                String [] items1 = new String[items.length];
                        //0 3 5
                for(int i = 0; i < items.length ; i++){
                    if(!items[i].equals("")){
                       //System.out.println(i);
                       items1[ii] = items[i].trim();
                     /*                      
                       if(i == 0){
                            items1[ii] = items[i];
                       }
                       if(i > 1 ){
                            items1[ii] = items[i]; 
                       }
                      */
                    ii++;
                    
                    if(ii >1)
                        break;
                    }
                    
                   // System.out.println(items1.length);
                }
                
               //System.out.println(items1[0] + " " + items1[1]);
               // x = nodes.indexOf(items1[0].trim());
               // y = nodes.indexOf(items1[1].trim());
                x = Integer.parseInt(items1[0].trim()) ; //Replacing Kishan's stupid code to match with file ID
                y= Integer.parseInt(items1[1].trim()) ;
                //linkWeight = nodes.indexOf(items[2].trim());
                linkWeight = Integer.parseInt(items[2].trim());
                
               // x = Integer.parseInt(items1[0].trim()) ; //Replacing Kishan's stupid code to match with file ID
              // y= Integer.parseInt(items1[1].trim()) ;
                //System.out.println(x + " " + y );
                try {
                	if(x!=y && linkWeight > THRESHOLD) //leave self links out: 
                	{
                		net.addLink(x,y);
                		
                		
                		
                		//Namal's code
                		//nodeList[nodeID1].addAdjacendyNode(nodeID2);
                		Node nodex = net.getNode(x);
                		nodex.addAdjacendyNode(y);
                		Node nodey = net.getNode(y);
                		nodey.addAdjacendyNode(x);
                		
                	}
                } catch (Exception ex) {
                ex.printStackTrace();
                return 0;
                }
                }
            scanner.close();
        } catch (IOException ex){
            return 0;
        }
        return 1;
        
        
      /*  
        for(int i=0;i <= lmap.length; i++){
            x = nodes.indexOf(lmap[i][0]);
            y = nodes.indexOf(lmap[i][1]);
            try {
                net.addLink(x,y);
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        
        */
    }
}



