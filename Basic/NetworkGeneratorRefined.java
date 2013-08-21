package Basic;

import java.io.BufferedWriter;
import java.text.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.Vector;

import org.apache.commons.math3.stat.Frequency;

import Distributions.EDistribution;
import Distributions.EDistributionGenerator;
import Distributions.QDistribution;
import Distributions.QDistribution_PowerLaw;
import oldClassesArchive.AssortativenessCalculator;

/**
 * This class is used to generate a network with given topology
 * Networks are hard coded
 * @author Piraveenan
 *
 */
public class NetworkGeneratorRefined {
	
	public Network growNetwork(int size, double joiningp) throws Exception{

		Network net = new Network();
		Network authorNet = new Network();
		BufferedWriter out = null;

		Random rand = new Random();
		//First Node
		Node firstNode  = new Node(0);
		firstNode.assignWeight(generateRandom(0.0,20.0)); //weight is equivalent to the impact factor
		net.addNode(firstNode);
		firstNode.authorArray = new Integer[3];
		
		Node firstAuthor = new Node(0);//add a new author to authornet
		authorNet.addNode(firstAuthor);
		Node secondAuthor = new Node(1);
		authorNet.addNode(secondAuthor);
		authorNet.addLink(firstAuthor,secondAuthor);
		Node thirdAuthor = new Node(2);
		authorNet.addNode(thirdAuthor);
		authorNet.addLink(secondAuthor,thirdAuthor);
		authorNet.addLink(firstAuthor, thirdAuthor);
		Integer [] au = {0,1,2};
		firstNode.authorArray=au;
		firstAuthor.papers.add(firstNode.getID());
		secondAuthor.papers.add(firstNode.getID());
		thirdAuthor.papers.add(firstNode.getID());
		int newAuthors=0;
		
		BufferedWriter wrOut = null;
		wrOut = new BufferedWriter(new FileWriter("f://CNRes/authorpIndex.txt"));
		BufferedWriter wrhOut = null;
		wrhOut = new BufferedWriter(new FileWriter("f://CNRes/authorhIndex.txt"));
		BufferedWriter wrpOut = null;
		wrpOut = new BufferedWriter(new FileWriter("f://CNRes/authorpRank.txt"));
		

		//Add the rest of the nodes according to PA
		double p = joiningp; //Joining probability - a parameter
		for(int i = 1; i < size; i++)
		{
			if (i%200==0){
				System.out.println("here I go: "+i);
				Node node = new Node(i);
				node.assignWeight(generateRandom(0.0,5.0)); //weight is equivalent to the impact factor
				int[] authorsArray = generateAuthors(i,node.authorArray.length);//get the new authors and old authors count
				newAuthors=authorsArray[node.authorArray.length];
				//Integer[] randArray = new Integer [node.authorArray.length];
				node.authorArray[0]=0;
				for (int y=1;y<(node.authorArray.length-authorsArray[node.authorArray.length]);y++){//pick old authors randomly
					int gen = generateRandom(1,(authorNet.getSize()-1));
					Integer auth= authorNet.getAllNodes().elementAt(gen).ID;
					while(Arrays.asList(node.authorArray).contains(auth)){
						auth= authorNet.getNode(generateRandom(1,authorNet.getSize()-1)).getID();
					}
					node.authorArray[y]=auth;
					authorNet.getNode(auth).papers.add(node.getID());
				}
				
				for(int z=(node.authorArray.length-authorsArray[node.authorArray.length]);z<node.authorArray.length;z++){//create nodes for new authors
					node.authorArray[z] = authorNet.getSize();
					Node author = new Node(authorNet.getSize());
					authorNet.addNode(author);
					author.papers.add(node.getID());
				}
				
				
				//node.authorArray=randArray;
				
				//System.out.println("Iteration: "+i+" AuthorNet Size: "+authorNet.getSize()+" New Authors: " + newAuthors);
				
				for (int j=0;j<node.authorArray.length;j++){
					for (int k=j+1;k<node.authorArray.length;k++){
						//System.out.println(j+" and "+k);
						//System.out.println(authorNet.getNode(j).getID());
						//System.out.println(node.authorArray[1]);
						if(!(authorNet.getNode(node.authorArray[j]).isLinked(authorNet.getNode(node.authorArray[k])))){
							try 
							{
								//System.out.println("WOOT");
								authorNet.addLink(authorNet.getAllNodes().elementAt(node.authorArray[j]), authorNet.getAllNodes().elementAt(node.authorArray[k]));
							} catch (Exception e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
					}
				}
				
				int citationsToNetwork = (int)(citationToNetwork(i,node.numberOfCitations)*node.numberOfCitations);//citations to network as a function of time
				int links=0;
				for(int j = 0; j <authorNet.getNode(0).papers.size() && links<citationsToNetwork;j++)
				{
					Node destNode = net.getNode(authorNet.getNode(0).papers.get(j));
					
						try 
						{
							net.addLink(node, destNode);
							//System.out.println(node.getID()+" cited "+destNode.getID());
							links++;
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					
					//System.out.println("CitationsToNetwork = "+j);
				}

				net.addNode(node);	
				ArrayList <Double> pRank=calculatepRank(authorNet,net);//calculate page rank
				wrOut.write(i+" ");
				wrhOut.write(i+" ");
				wrpOut.write(i+" ");
				
				ArrayList <Integer> hIndex=calculatehIndex(authorNet,net);//calculate h-index
				for (int j = 0; j < authorNet.getSize(); j++) {
					wrOut.write(authorNet.getNode(j).percentile + " ");
					wrhOut.write(authorNet.getNode(j).hIndex + " ");
					wrpOut.write(pRank.get(j).toString() + " ");
					
				}
				
				wrOut.write("\n");
				wrhOut.write("\n");
				wrpOut.write("\n");
			}
			
			else {
				Node node = new Node(i);
				node.assignWeight(generateRandom(0.0,20.0)); //weight is equivalent to the impact factor
				int[] authorsArray = generateAuthors(i,node.authorArray.length);//get the new authors and old authors count
				newAuthors=authorsArray[node.authorArray.length];
				//Integer[] randArray = new Integer [node.authorArray.length];
				for (int y=0;y<(node.authorArray.length-authorsArray[node.authorArray.length]);y++){//pick old authors randomly
					int gen = generateRandom(0,(authorNet.getSize()-1));
					//System.out.println(gen);
					//System.out.println(authorNet.getNode(gen).getID());
					Integer auth= authorNet.getAllNodes().elementAt(gen).ID;
					//System.out.println("auth is: "+auth + " gen is: "+gen);
					while(Arrays.asList(node.authorArray).contains(auth)){
						auth= authorNet.getNode(generateRandom(1,authorNet.getSize()-1)).getID();
					}
					node.authorArray[y]=auth;
					authorNet.getNode(auth).papers.add(node.getID());
				}
				
				for(int z=(node.authorArray.length-authorsArray[node.authorArray.length]);z<node.authorArray.length;z++){//create nodes for new authors
					node.authorArray[z] = authorNet.getSize();
					Node author = new Node(authorNet.getSize());
					authorNet.addNode(author);
					author.papers.add(node.getID());
				}
				
				
				//node.authorArray=randArray;
				
				//System.out.println("Iteration: "+i+" AuthorNet Size: "+authorNet.getSize()+" New Authors: " + newAuthors);
				
				for (int j=0;j<node.authorArray.length;j++){
					for (int k=j+1;k<node.authorArray.length;k++){
						//System.out.println(j+" and "+k);
						//System.out.println(authorNet.getNode(j).getID());
						//System.out.println(node.authorArray[1]);
						if(!(authorNet.getNode(node.authorArray[j]).isLinked(authorNet.getNode(node.authorArray[k])))){
							try 
							{
								//System.out.println("WOOT");
								authorNet.addLink(authorNet.getAllNodes().elementAt(node.authorArray[j]), authorNet.getAllNodes().elementAt(node.authorArray[k]));
							} catch (Exception e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
					}
				}
				
				double probs[] =  new double[net.getSize()];
				int citationsToNetwork = (int)(citationToNetwork(i,node.numberOfCitations)*node.numberOfCitations);//citations to network as a function of time

				for(int j = 0; j <net.getSize();j++)
				{
					Node destNode = net.getAllNodes().elementAt(j);
					double prob = p;
					if(net.getNoOfLinks() > 0 && destNode.getNumberOfLinks() > 0)
					{
						//prob = p*((double)destNode.getNumberOfLinks()) / (double)net.getNoOfLinks();
						prob = (prob*destNode.getWeight()* (double)destNode.getNumberOfInLinks()) / (double)net.getWeightedDegreeSum(); //For weighted PA
					}
					else //No links made yet
					{
						//prob = p;
						prob = prob*destNode.getWeight(); //For weighted PA
					}
					probs[j] = prob;
				}
				int links=0;
				for(int j = 0; j <probs.length && links<citationsToNetwork;j++)
				{
					Node destNode = net.getAllNodes().elementAt(j);
					//Throw dice
					if(rand.nextDouble() <=  probs[j])
					{
						try 
						{
							net.addLink(node, destNode);
							//System.out.println(node.getID()+" cited "+destNode.getID());
							links++;
						} catch (Exception e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					//System.out.println("CitationsToNetwork = "+j);
				}

				net.addNode(node);
				//System.out.println("Number of Authors:" + node.authorArray.length);
				//System.out.println("Number of Citations:" + node.numberOfCitations + " Number of Citations to the Network:" + citationsToNetwork);

				//net.createTextFile("f://CNRes/paperNet_"+i+".txt");
				//authorNet.createTextFile("f://CNRes/authorNet_"+i+".txt");
				
				
				ArrayList <Double> pRank=calculatepRank(authorNet,net);//calculate page rank
				//System.out.println(pRank.toString());
				wrOut.write(i+" ");
				wrhOut.write(i+" ");
				wrpOut.write(i+" ");
				
				ArrayList <Integer> hIndex=calculatehIndex(authorNet,net);//calculate h-index
				//System.out.println(hIndex.toString());
				//net.writeIndices(authorNet, "f://CNRes/authorIndex_"+i+".txt");
				for (int j = 0; j < authorNet.getSize(); j++) {
					wrOut.write(authorNet.getNode(j).percentile + " ");
					wrhOut.write(authorNet.getNode(j).hIndex + " ");
					wrpOut.write(pRank.get(j).toString() + " ");
					
				}
				
				wrOut.write("\n");
				wrhOut.write("\n");
				wrpOut.write("\n");
			}
		}
		
		System.out.println("Papernet size nodes "+net.getSize()+ " links "+ net.getNoOfLinks());
		System.out.println("Authornet size nodes "+authorNet.getSize()+ " links "+ authorNet.getNoOfLinks());
		wrOut.flush();
		wrhOut.flush();
		wrpOut.flush();
		return net;

	}
	
	public ArrayList<Double> calculatepRank(Network authors, Network net) throws Exception{
        
		for(int i=0;i<net.getSize();i++){
			Node thisNode = net.getNode(i);
			thisNode.updatePageRank(0.5, net.getSize());
		}
		
		ArrayList <Double> pageRankSum = new ArrayList <Double>();
		
		for (int j = 0; j < authors.getSize(); j++) {
			ArrayList<Integer> papers = authors.getNode(j).papers;
			double pr = 0.0;
			for (int k = 0; k < papers.size(); k++) {
				pr = pr + net.getNode(papers.get(k)).getPageRank();//summing up the page rank values for each paper
			}
			pageRankSum.add(pr);
		}
		Frequency newf = new Frequency();
		ArrayList<Double> ret = new ArrayList<Double>();
		for(int x=0;x<pageRankSum.size();x++){
			newf.addValue(pageRankSum.get(x));
		}
		for(int x=0;x<pageRankSum.size();x++){
			int dec = (int)(newf.getCumPct(pageRankSum.get(x))*1000);
			double d = dec/10.0;
			ret.add(d);
			authors.getNode(x).percentile=d;
		}
		
		//System.out.println(pageRankSum.toString());
		
		return pageRankSum;
	}
	
	public ArrayList<Integer> calculatehIndex(Network authors, Network papernet) throws Exception{
		
		ArrayList <Integer> hIndex = new ArrayList<Integer>();
		
		for (int i=0;i<authors.getSize();i++){
			ArrayList <Integer> papers = authors.getNode(i).papers;
			ArrayList <Integer> nCite = new ArrayList <Integer>();
			for (int j=0;j<papers.size();j++){
				Node paperNode=papernet.getNode(papers.get(j));
				nCite.add(paperNode.getNumberOfInLinks());
			}
			Collections.sort(nCite);
			Collections.reverse(nCite);
			//System.out.println("Author: "+i+" "+nCite.toString());
			//int hInd=0;
			int k=0;
			for (k=0;k<nCite.size();k++){
				if (k==0 && nCite.get(k)==0) break;
				if (k >=nCite.get(k)) break;
			}
			authors.getNode(i).hIndex=k;
			hIndex.add(k);
		}
		return hIndex;
	}

	public int generateRandom(int min, int max){
		return min + (int)(Math.random() * ((max - min) + 1)); //return a random number of citations to be initiated.
	}
	
	public int generateRandome(int min, int max){
		return min + (int)(Math.random() * ((max - min))); //return a random number of citations to be initiated exclusive of max.
	}

	public double generateRandom(double min, double max){
		return min + (Math.random() * ((max - min) + 1)); //return a random number of citations to be initiated.
	}

	public double citationToNetwork(int i, int j){
		if (i<500){
			return 0.7/500 * Math.max(i, j);
		}
		else
			return 0.7;
	}

	public int[] generateAuthors(int i,int n){//n=number of authors for the paper, i=iteration variable to determine the percentage of old authors
		int [] authors = new int[n+1];
		int newA=0;
		int oldA=0;
		for(int j=0;j<n;j++){
			double r = generateRandom(0,1);
			if(i<500){
				if(r<(0.2*i/500)){
					oldA++;
					authors[j]=0;//pick an old author randomly
				}
				else
					newA++;
					authors[j]=1;//generate a new author
			}
			else if (i>=500){
				if(r<0.2){
					oldA++;
					authors[j]=0;//pick an old author randomly
				}
				else
					newA++;
					authors[j]=1;//generate a new author
			}
		}
			
		
		authors[n]=newA;
		return authors;
		
	}
	
}
