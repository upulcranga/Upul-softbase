package Basic;

/**
 * This class generates some random and Preferentiall attachemnt networks and prints them to files
 * @author pir011
 *
 */
public class GenericNetGenerator {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String location = ".\\LA-Extended\\simplenets\\";
		
		NetworkGenerator gen = new NetworkGenerator();
		
		Network net;
		/*
		net  = gen.generateERRandom(1000);
		net.createTextFile(location + "ra1");
		
		 
		
		net  = gen.generateERRandom(2000);
		net.createTextFile(location +"ra2");
		
		net  = gen.generateERRandom(3000);
		net.createTextFile(location +"ra3");
		
		net  = gen.generateERRandom(4000);
		net.createTextFile(location +"ra4");
		
		net  = gen.generateERRandom(5000);
		net.createTextFile(location +"ra5"); 
		
		net = gen.growNetworkPA(1000, 0.5);
		net.createTextFile(location +"pa1");
		

		net = gen.growNetworkPA(1000, 0.75);
		net.createTextFile(location +"pa2");
		

		net = gen.growNetworkPA(1000, 2);
		net.createTextFile(location +"pa3"); 
		

		net = gen.growNetworkPA(1000, 2.5);
		net.createTextFile(location +"pa4");
		
		net = gen.growNetworkPA(1000, 5);
		net.createTextFile(location +"pa5"); 
		
		net = gen.growNetworkAPA(1000, 40, 1.1, 0.8);
		net.createTextFile(location +"apa1");
		
		net = gen.growNetworkAPA(1000, 40, 1.1, 0.5);
		net.createTextFile(location +"apa2");
		
		net = gen.growNetworkAPA(1000, 40, 1.1, -0.8);
		net.createTextFile(location +"apa3");
		
		net = gen.growNetworkAPA(1000, 40, 1.1,- 0.5);
		net.createTextFile(location +"apa4");
		
		net = gen.growNetworkAPA(1000, 40, 1.1, -0.1);
		net.createTextFile(location +"apa5"); */
		
		//net = gen.growIGNet(1000);
		//net.createTextFile(location +"starIG2");
		
		net = gen.growNetworkPA(12, 0.5);
		net.createTextFile(location +"starPA");
		
		/*

		net = gen.growIGNet(1000);
		net.createTextFile(location +"ig2");

		net = gen.growIGNet(1000);
		net.createTextFile(location +"ig3");
		

		net = gen.growIGNet(1000);
		net.createTextFile(location +"ig4");
		

		net = gen.growIGNet(1000);
		net.createTextFile(location +"ig5"); */
		
		System.out.println("DONE!");
	}

}
