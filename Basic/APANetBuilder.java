package Basic;

import java.util.Random;

import javax.swing.*;

import java.awt.*;
import java.util.*;
import java.awt.event.*;
import java.io.*;

import oldClassesArchive.AssortativenessCalculator;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartPanel;
import java.awt.Color;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.event.*;
import org.jfree.chart.renderer.DefaultPolarItemRenderer;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.StandardCategoryItemLabelGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.CategoryItemRenderer;
import org.jfree.chart.renderer.category.LineAndShapeRenderer;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.ui.RectangleInsets;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.category.LineAndShapeRenderer;
import org.jfree.chart.axis.DateAxis;

import Distributions.EDistribution;
import Distributions.QDistribution_PowerLaw;

public class APANetBuilder {

	private ChartPanel cp;

	private DefaultCategoryDataset eDistData;

	private JFreeChart eDistChart;

	private String series1;

	private String series2; // r=1

	private String series3; // r=0

	private String series4; // r=-1

	private String netseries1; // E r

	private String netseries2; // E r net

	private String netseries3; // Q

	private String netseries4; // Q net

	private ChartPanel cp3;

	private DefaultCategoryDataset netEDistData;

	private JFreeChart netEDistChart;

	private ChartPanel cp4;

	private DefaultCategoryDataset netQDistData;

	private JFreeChart netQDistChart;

	private ChartPanel cp2;

	private DefaultCategoryDataset qDistData;

	private JFreeChart qDistChart;

	private String seriesq;

	private JTextField prField;

	private JTextField exrField;

	private JTextField rField;

	private JTextField maxMinusRField;

	private JTextField maxAvgKField;

	private JTextField exavgkField;

	private JTextField exnField;

	private JTextField exGammaField;

	private JTextField avgkField;

	private JTextField nField;

	private JTextField miField;

	private JTextField indField;

	private JButton buildNetButton;

	// private QDistribution_PowerLaw qDistribution;
	private double avg_k;

	private double gamma;

	private final int MAXLENGTH = 100;

	private BufferedWriter rwriter;

	private BufferedWriter iwriter;

	public static void main(String[] args) {

		/* Get an instance of the simulator */
		APANetBuilder sim = new APANetBuilder();

		// Get the Q distribution for the desired average degree
		double res = 0.01;
		double gamma_res = 0.025;
		double gamma = 0.5 - gamma_res;
		BufferedWriter out = null;
		try {
			out = new BufferedWriter(new FileWriter("gammas3.xls"));
		} catch (IOException e) {
			System.out.println("cant write");
			e.printStackTrace();
		}
		System.out.println("STARTING NEWMAN TESt");
		while (gamma < 0) {
			gamma = gamma + gamma_res;
			QDistribution_PowerLaw qDist = generateQ(gamma);

			// CHECK THE NEWMAN MINIMIZATION - THIS PART DOES NOT BELONG HERE
			// AND SHOULD BE MADE INTO ANOTHER CLASS
			// //////////////////////////////////////////////////////////////////////////
			System.out.print("STARTING NEWMAN TEST: gamma " + gamma + "  ");
			double minR = 0;

			double[] x = new double[qDist.getLength()];
			double[] q = qDist.getDist();
			Random rand = new Random();
			for (int ii = 0; ii < 1; ii++) // brute force change x distribution
			{
				if (ii % 1000000 == 0) {
					// System.out.println("looping "+ii);
				}
				double sum = 0;
				int cou = 0;
				double p = 0;

				double[] xrand2 = new double[qDist.getLength()];
				if (qDist.getLength() == 6) {
					for (double i1 = 0; i1 < 1.0; i1 = i1 + res) {
						xrand2[1] = i1;
						for (double i2 = 0; i2 < 1.0; i2 = i2 + res) {
							xrand2[2] = i2;
							for (double i3 = 0; i3 < 1.0; i3 = i3 + res) {
								xrand2[3] = i3;
								for (double i4 = 0; i4 < 1.0; i4 = i4 + res) {
									xrand2[4] = i4;
									for (double i5 = 0; i5 < 1.0; i5 = i5 + res) {
										xrand2[5] = i5;

										if ((xrand2[1] + xrand2[2] + xrand2[3] + xrand2[4]+ xrand2[5]) < 1 - (0.5 * res)) {

											continue;
										} else if ((xrand2[1] + xrand2[2]
												+ xrand2[3] + xrand2[4]+ xrand2[5]) > 1 + (0.5 * res)) {

											continue;
										} else {
											// System.out.print(true);
											cou++;

											EDistribution Ed = generateNewmanEdist(
													q, xrand2);
											boolean b = checkEDistValidity(Ed,
													qDist);
											if (b) {
												// System.out.println("VALID
												// X");
												double mq = weightedmean(q);
												double mx = weightedmean(xrand2);
												double sq2 = variance(q);
												// calculate rd
												double rd = -((mq - mx) * (mq - mx))
														/ sq2;
												// System.out.println("NEWMAN R
												// " + rd);
												if (rd < minR) {
													minR = rd;
												}

											} else {
												// System.out.println("INVALID
												// X");
											}
										}
									}
								}

							}
						}

					}
				} else {
					System.out.println("length not five");
				}
				if (qDist.getLength() == 4) {

				}

				// double xQuantum = (1.0 / x.length*1000); //In total, x must
				// have xlength*1000 quantim units
				// populate x;
				/*
				 * double[] xrand = new double[x.length]; double ssum = 0;
				 * xrand[0] = 0; for(int i = 1; i < (x.length ); i++) { xrand[i] =
				 * Math.abs(rand.nextDouble()) ; ssum = ssum + xrand[i];
				 *  } for(int i = 1; i < (x.length ); i++) { xrand[i] = xrand[i] /
				 * ssum; }
				 */

				/*
				 * int startInt = Math.abs(rand.nextInt()) % x.length; for(int i =
				 * startInt; i < (x.length + startInt); i++) //brute force
				 * change x distribution { if(i < (x.length + startInt-1)) p =
				 * (1.0-sum)*rand.nextDouble(); else //last element p= 1.0 -sum;
				 * x[i%x.length] = p;//(1.0 / x.length); sum = sum + p; }
				 */

			}

			System.out.println("NEWMAN R min " + minR);
			try {
				out.write("" + gamma + "\t" + minR + "\n");
				out.flush();
			} catch (IOException e) {
				System.out.println("no file");
				e.printStackTrace();
			}

			// System.out.println("COMPLETED NEWMAN TEST");

			// ////////////////////////////////////////////////////////////////////////////
		}
	}

	/**
	 * assuming x is a probability distribution, always between 0 and 1 and
	 * summing up to 1
	 * 
	 * @param x
	 * @return
	 */
	public static double weightedmean(double[] x) {
		double m = 0;
		for (int i = 0; i < x.length; i++) {
			m = m + (i * x[i]);
		}

		return m;
	}

	public static double variance(double[] x) {
		double m1 = 0;
		double m2 = 0;
		for (int i = 0; i < x.length; i++) {
			m1 = m1 + (i * i * x[i]);
			m2 = m2 + (i * x[i]);
		}

		return m1 - (m2 * m2);
	}

	public void build() {

		Random random = new Random();
		Network network;

		int SLIDINGSTEPS = 20;

		/* Read Parameters from Screen */
		double ra = Double.parseDouble(exrField.getText());
		int N = Integer.parseInt(exnField.getText());
		avg_k = Double.parseDouble(exavgkField.getText());
		gamma = Double.parseDouble(exGammaField.getText());

		/* Get the Q distribution for the desired average degree */
		QDistribution_PowerLaw qDist = generateQ();

		/*
		 * If specified assortativeness is minus, Calculate what is the maximum
		 * minus assortativeness possible for this Q distribution, and scale the
		 * user specified r accordingly.
		 */
		double maxMinusR = Math.round(1000 * maximumMinusR(qDist)) / 1000.0;
	
		maxMinusRField.setText("" + (Math.round(1000 * maxMinusR) / 1000.0)
				+ "");
		exrField.setBackground(Color.WHITE);
		double scaledR = ra;
		if (scaledR < 0) {
			if (scaledR < maxMinusR) {
				scaledR = maxMinusR;
				exrField.setText("" + scaledR);
				exrField.setBackground(Color.RED);
			}

			scaledR = scaledR * (1.0 / Math.abs(maxMinusR)); /*
																 * r is scaled
																 * up to work
																 * with E(j,k)
																 * generation
																 */

		}

		/* Get Edistribution for desired r and Q distribution */
		EDistribution eDist = generateEDistribution(qDist, scaledR,
				SLIDINGSTEPS);

		/* Update display of the given Q distribution and E distribution */
		update_EDistDataset(eDist, qDist);

		/*
		 * Build up a network of N nodes the APA way, for the given
		 * EDistribution and Q distribution
		 */

		// sim.indField.setText("BUILDING NETWORK");
		network = buildAPANetwork(qDist, eDist, N);
		// sim.indField.setText("FINISHED NETWORK");

		/* Update display of the grown Q dsitribution and E distribution */
		update_NETEDistDataset(eDist, network.getLinkDistribution());
		update_QDistDataset(qDist, network.getRemainingDegreeDistribution());

		/* Measure r for the grown network, and display it, along with error */
		AssortativenessCalculator ac = new AssortativenessCalculator(network);
		double calculated_r = Math.round(1000 * ac.calculate_r()) / 1000.0;
		double calculated_avg_k = Math.round(1000 * network.getAverageDegree()) / 1000.0;
		System.out.println("R: " + calculated_r);

		System.out.println("Links " + ((network.getAllLinks().size())));
		System.out.println("Nodes " + ((network.getAllNodes().size())));

		ra = Math.round(1000 * ra) / 1000.0;
		rField.setText("" + calculated_r);

		double percentage_error = Math.abs(100 * (ra - calculated_r)
				/ calculated_r);
		percentage_error = Math.round(1000 * percentage_error) / 1000.0;
		prField.setText("" + percentage_error + "%");
		nField.setText("" + network.getSize());
		avgkField.setText("" + calculated_avg_k);
		// maxAvgKField.setText(""+(Math.round(1000*qDist.calculateMaximumAverageDegree(gamma,
		// MAXLENGTH))/1000.0)+"");

		double mutualInformation = network.calculateMI();

		mutualInformation = Math.round(1000 * mutualInformation) / 1000.0;
		miField.setText("" + mutualInformation);
		System.out.println("I(Q): " + mutualInformation);

		try {

			rwriter.write("" + calculated_r);
			rwriter.flush();
			rwriter.newLine();
			iwriter.write("" + mutualInformation);
			iwriter.flush();
			iwriter.newLine();
		} catch (Exception e) {
			// Eat it!
		}
		/* Display The grown Network */
		// TODO: This needs to be implemented

	}

	/* Constructor */
	public APANetBuilder() {
		try {
			rwriter = new BufferedWriter(new FileWriter("apa.xls"));
			iwriter = new BufferedWriter(new FileWriter("apa2.xls"));

		} catch (Exception e) {

		}
		/* Set the default average degree and gamma values */
		avg_k = 2.6;
		gamma = 1;

		/* Set data series for JFreeChart */
		series1 = "Current E Distribution";
		series2 = "ED for r=1";
		series3 = "ED for r=-1";
		series4 = "ED for r=0";
		netseries1 = "Theoritic E Distribution";
		netseries2 = "Grown E distribution";
		netseries3 = "Theoritic Q distribution";
		netseries4 = "Grown q distribution";
		seriesq = "QD";
		eDistData = new DefaultCategoryDataset();
		qDistData = new DefaultCategoryDataset();
		netEDistData = new DefaultCategoryDataset();
		netQDistData = new DefaultCategoryDataset();
		JPanel displayPanel = new JPanel();

		/* Generate ResultPanel */
		JPanel rPanel = resultPanel();

		/* Set initial Values */
		exavgkField.setText("" + avg_k);
		exrField.setText("0.0");
		exnField.setText("9000");
		miField.setText("0");
		exGammaField.setText("" + gamma);
		maxMinusRField
				.setText(""
						+ (Math.round(1000 * maximumMinusR(generateQ(avg_k,
								gamma))) / 1000.0));
		maxAvgKField.setText(""
				+ generateQ(avg_k, gamma).calculateMaximumAverageDegree(gamma,
						MAXLENGTH));

		/* Create other panels, including those requiring a Q Distribution */
		displayPanel.add(eDistPanel());
		displayPanel.add(qDistPanel());
		displayPanel.add(rPanel);
		displayPanel.add(netEDistPanel());

		/* Create Frame */
		JFrame frame = new JFrame();
		frame.setSize(1200, 1000);
		frame.add(displayPanel);
		frame.show();

	}

	public QDistribution_PowerLaw generateQ(double avg_k, double gamma) {

		// return generateQ();

		QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(avg_k, gamma);
		if (maxAvgKField != null) {
			maxAvgKField.setText(""
					+ (Math.round(1000 * qDist.calculateMaximumAverageDegree(
							gamma, MAXLENGTH)) / 1000.0) + "");
		}
		return qDist;

	}

	public QDistribution_PowerLaw generateQ() {
		// return tent
		double[] qArr = new double[6];
		qArr[0] = 0;
		qArr[1] = (0.1);
		qArr[2] = (0.2);
		qArr[3] = (0.4);
		qArr[4] = (0.2);
		qArr[5] = (0.1);
		// QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(qArr);

		QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(avg_k, gamma);
		if (maxAvgKField != null) {
			maxAvgKField.setText(""
					+ (Math.round(1000 * qDist.calculateMaximumAverageDegree(
							gamma, MAXLENGTH)) / 1000.0) + "");
		}
		return qDist;

	}

	public static QDistribution_PowerLaw generateQ(double gamma) {

		/*
		 * double[] qArr = new double[5]; qArr[0] = 0; qArr[1] = (1.0/6.0);
		 * qArr[2] = (2.0/6.0); qArr[3] = (1.0/6.0); qArr[4] = (2.0/6.0);
		 * QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(qArr);
		 */

		QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(6, gamma, true);
		//QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(3, gamma);

		return qDist;

	}

	/**
	 * 'Scale' a Edistribution
	 */
	private double[] getScaledEDist(int joiningNodeQ, EDistribution eDist,
			QDistribution_PowerLaw qDist) {
		if (joiningNodeQ < 0) {
			System.out
					.println("A Node Tries to join with Unsuitable remaining Degree: Error");
		}
		double[] scaledDist = new double[eDist.getDistribution().length];

		/* Copy the values from main E distribution */
		double probSum = 0;
		for (int i = 0; i < scaledDist.length; i++) {
			
				scaledDist[i] = (eDist.getProbability(joiningNodeQ, i) / qDist
						.getDegreeProbability(i + 1));
			

			probSum = probSum + scaledDist[i];
		}

		/*
		 * Scale the copied Distribution so that summation of its probabilities
		 * become 1
		 */
		for (int i = 0; i < scaledDist.length; i++) {
			scaledDist[i] = scaledDist[i] / probSum;
		}

		return scaledDist;
	}

	/**
	 * This method implements the Assortative Preferential Attachment Method
	 * 
	 * @param qDist
	 * @param eDist
	 * @param networkSize
	 * @return
	 */
	public Network buildAPANetwork(QDistribution_PowerLaw qDist, EDistribution eDist,
			int networkSize) {
		Random random = new Random();

		/* start with m0 nodes. 1:2 ratio between m0 and N is optimal */
		int m0 = (int) (0.5 * networkSize);

		Network network = new Network();

		/* Add half of the nodes with no links */
		/* Gennerate otehr half of the nodes, which will be joining one by one */
		int counter = 0;
		Vector<Node> joiningNodes = new Vector();
		for (int y = 0; y < qDist.getLength(); y++) {

			int nodesWithThisq = (int) Math.round((m0 * qDist
					.getDegreeProbability(y + 1)));
			if (y == qDist.getLength() - 1)// *If the last one, get the balance
											// of nodes */
			{
				nodesWithThisq = m0 - network.getSize();
			}
			for (int k = 1; k <= nodesWithThisq; k++) {
				counter++;
				Node node = new Node(counter, (y + 1));
				network.addNode(node);
				/* Generate another node, with the same degree which will join */
				Node node2 = new Node(m0 + counter, (y + 1));
				joiningNodes.addElement(node2);
			}

		}

		/* Shuffle the set of nodes that are going to join */
		joiningNodes = network.shuffle(joiningNodes);

		/* Add this half of nodes in the 'APA' way */

		int addedLinksCount = 0;
		boolean linkLimit;

		for (int i = 0; i < (m0); i++) /*
										 * Note that m0 = N/2 gives the optimum
										 * conditions
										 */
		{
			linkLimit = false;
			addedLinksCount = 0;

			/* Add new node */
			// Node node = new Node(m0+i, qDist.getRandomDegreeDUMMY(i));
			// /*TODO: Hck random change back */
			Node node = joiningNodes.elementAt(i); // new Node(m0+i,
													// qDist.getRandomDegree());
			network.addNode(node);

			/* ADD LINK */
			/* scale the eDist for this instance of adding */
			double[] scaledDist = getScaledEDist(node.getAssignedRemainingDegree(),
					eDist, qDist);

			/*
			 * Go through all existing nodes in the network in Random order, and
			 * assign probabilities
			 */
			network.shuffle();
			for (int j = 0; j < network.getSize(); j++)
			// while(true)
			{
				Node networkNode = network.getAllNodes().elementAt(j);
				// Node networkNode = network.getAllNodes().elementAt(
				// (Math.abs(random.nextInt())%network.getSize()) ); //TODO:
				// Check for duplicate links
				double nodeProb = scaledDist[networkNode.getAssignedRemainingDegree()];
				// double nodeProb =
				// scaledDist2[node.getRemainingDegree()][networkNode.getRemainingDegree()];
				if (networkNode.getNumberOfLinks() >= networkNode.getAssignedDegree()) {
					/* all 'Stumps' have been exhausted */
					nodeProb = 0;
				}

				/*
				 * if(networkNode.isLinked(node)) { // This network node is
				 * already linked to joining node nodeProb = 0; }
				 */

				if (networkNode.getID() == node.getID()) {
					/* We dont want to connect the joining ndoe to itself! */
					nodeProb = 0;
				}

				/*
				 * Toss a coin to see whether this node should connect to the
				 * new node
				 */
				if (random.nextDouble() < nodeProb) {
					try {
						/* Add Link */
						network.addLink(node, networkNode);
						addedLinksCount++;
						if (addedLinksCount >= node.getAssignedDegree()) {
							/* Should not add any more links */
							linkLimit = true;
							addedLinksCount = 0;
						}
					} catch (Exception e) {
						// Eat it!
						System.out.println("Problem adding links");
						linkLimit = true;
						addedLinksCount = 0;
					}
				}

				if (linkLimit) /*
								 * No need to add any more links to this joining
								 * node
								 */
				{
					break;
				}

			}

			/* Measure r */
			/*
			 * if( (eDist.Expectation()<4.48) && (eDist.Expectation()>4.38) ){
			 * AssortativenessCalculator ac = new
			 * AssortativenessCalculator(network); double ra = ac.calculate_r();
			 * System.out.println("r now: "+ra); }
			 */

			/* update display */
			// update_NETEDistDataset(eDist, network.getLinkDistribution());
			// update_QDistDataset(qDist, network.getNodeDistribution());

		}
		// network.printArray(network.getLinkDistribution(),"link dist");
		return network;

	}

	/**
	 * This method implements symetric sliding
	 * 
	 * @param qDist
	 * @param r
	 * @return
	 */
	public EDistribution generateEDistribution(QDistribution_PowerLaw qDist, double r,
			int totslidingSteps) {
		Random random = new Random();
		EDistribution EDist1 = generateAssortativeEDistribution(qDist);
		EDistribution EDistM1 = generateDisassortativeEDistribution(qDist);
		EDistribution EDist = generateAssortativeEDistribution(qDist);
		EDistribution EDist0 = generateNonassortativeEDistribution(qDist);
		int totalOneSideSlidingSteps = totslidingSteps / 2;
		int slidingSteps;
		if (r >= 0) {

			slidingSteps = (int) Math
					.round(((double) totalOneSideSlidingSteps * (1.0 - r)));
		} else {
			slidingSteps = (int) Math
					.round(((double) totalOneSideSlidingSteps * (1.0 + r)));
		}

		/*
		 * slide down by 'slidingSteps' number of steps from the r=1
		 * eDistribution
		 */

		double[][] eArr;
		if (r >= 0) /* Sliding from r=1 to r=0 */
		{
			eArr = EDist.getDistribution();
			for (int i = 0; i < qDist.getLength(); i++) {

				double stepSize = (EDist1.getProbability(i, i) - EDist0
						.getProbability(i, i))
						/ (double) totalOneSideSlidingSteps;
				eArr[i][i] = eArr[i][i] - (stepSize * slidingSteps);

				for (int j = 0; j < eArr[i].length; j++) {

					if (i == j) {
						continue;
					}

					eArr[i][j] = eArr[i][j]
							+ (((double) slidingSteps / (double) totalOneSideSlidingSteps) * EDist0
									.getProbability(i, j));

				}
			}
		} else /* Sliding from r=0 to r=-1 */
		{

			// EDist = generateDisassortativeEDistribution(qDist);
			eArr = EDistM1.getDistribution();
			
			for (int i = 0; i < qDist.getLength(); i++) // Note the index
														// starting from 1
			{

				// double stepSize = (
				// EDist1.getProbability(i,i)-EDist0.getProbability(i,i)) /
				// (double)totalOneSideSlidingSteps;
				double stepSize = (EDistM1.getProbability(i,
						(qDist.getLength()-1 - i)) - EDist0.getProbability(i,
						(qDist.getLength()-1 - i)))
						/ (double) totalOneSideSlidingSteps;

				// eArr[i][qDist.getLength()-i] = eArr[i][qDist.getLength()-i] -
				// (stepSize*slidingSteps);

				for (int j = 0; j < eArr[i].length; j++) {

					if ((i + j+1) == qDist.getLength()) /*
														 * Dis assorttative
														 * peaks
														 */
					{
						double min = Math.min(EDist1.getProbability(i, i),
								EDist1.getProbability((qDist.getLength()-1 - i),
										(qDist.getLength()-1 - i)));
						stepSize = Math.abs((min - EDist0.getProbability(i,
								(qDist.getLength()-1 - i)))
								/ (double) totalOneSideSlidingSteps);
						eArr[i][j] = eArr[i][j] - (stepSize * slidingSteps);

					} else if (i == j) /* Assortative peaks */
					{
						double diff = 0;

						diff = EDistM1.getProbability(i, i)
								- EDist0.getProbability(i, i);

						stepSize = Math
								.abs((EDistM1.getProbability(i, i) - EDist0
										.getProbability(i, i))
										/ (double) totalOneSideSlidingSteps);
						// stepSize = Math.abs(( diff-EDist0.getProbability(i,i)
						// ) / (double)totalOneSideSlidingSteps);
						if (diff >= 0) {
							eArr[i][j] = eArr[i][j] - (stepSize * slidingSteps);
						} else {
							eArr[i][j] = eArr[i][j] + (stepSize * slidingSteps);
						}
					} else /* Not a peak at any point in time */
					{
						stepSize = Math
								.abs((EDistM1.getProbability(i, j) - EDist0
										.getProbability(i, j))
										/ (double) totalOneSideSlidingSteps);
						eArr[i][j] = eArr[i][j] + (stepSize * slidingSteps);
					}

				}
			}
		}

		EDistribution newEDist = new EDistribution(eArr);

		return newEDist;
	}

	/**
	 * This function generates a uniform degree distribution for a given average
	 * degree (<k>) average degree has to be an integer, otherwise a uniform
	 * distributio cannot be created.
	 */
	public static double[] getUniformDegreeDistribution(int average_k) {
		int distLength = 2 * (average_k - 1) + 1;
		double[] arr = new double[distLength];
		for (int i = 2; i < distLength; i++) {
			arr[i] = (1.0 / (double) (distLength - 2.0)); // Zeroth element,
															// first element not
															// needed
		}
		return arr;

	}

	/**
	 * NOTE: What is 'UNIFORM' is the degree distribution, NOT q distribution
	 * 
	 * @param average_k
	 * @return
	 */
	public static QDistribution_PowerLaw getUniformQDistribution(int average_k) {
		double[] pDist = getUniformDegreeDistribution(average_k);
		QDistribution_PowerLaw qDist = new QDistribution_PowerLaw(pDist, true);
		return qDist;
	}

	/**
	 * It is not possible for all Q Distributions to generate r=-1 Specifically,
	 * assymetric q distributions will not be able to do that. This method
	 * calculates the maximum minus r that can be generated for a given q (k)
	 * distribution
	 */
	public double maximumMinusR(QDistribution_PowerLaw qDist) {
		EDistribution EDist1 = generateAssortativeEDistribution(qDist);
		EDistribution EDistM1 = generateDisassortativeEDistribution(qDist);
		EDistribution EDist0 = generateNonassortativeEDistribution(qDist);

		double maxMinusR = Math.abs(EDistM1.Expectation()
				- EDist0.Expectation())
				/ Math.abs(EDist1.Expectation() - EDist0.Expectation());
		return -1 * maxMinusR;
	}

	/**
	 * This method generates the e distribution for the given q distribution and
	 * perfect assortativess (r=1)
	 * 
	 */
	public EDistribution generateAssortativeEDistribution(QDistribution_PowerLaw qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		for (int y = 0; y < qDist.getLength(); y++) {
			/* when both indices are the same, there is peak probability */
			eArr[y][y] = qDist.getRemainingDegreeProbability(y);
			/* Otherwise zero, but no need to explicitly code that as arrays are automatically initialized to zero */
		}

		EDistribution ed = new EDistribution(eArr);
		return ed;
	}

	/**
	 * This method generates the e distribution for the given q distribution and
	 * perfect dis-assortativess (r=-1) Note that if qDistribution is non
	 * symmetric, r=-1 is not possible. This method then generates a close
	 * approximation of best (minusmost) minus r TODO: This still does not
	 * generate the best minus r. Make it better!
	 * 
	 */
	public EDistribution generateDisassortativeEDistribution(QDistribution_PowerLaw qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		double[][] fakeEARR = new double[qDist.getLength()][qDist.getLength()];
		int lastIndex = qDist.getLength() - 1;
		double Nq = qDist.getLength();
		

		for (int y = 0; y < qDist.getLength(); y++) /*
													 * Note the indexes starting
													 * from 1 rather than 0.
													 * Need to stay this way
													 */
		{
			if (qDist.getRemainingDegreeProbability(y) <= qDist
					.getRemainingDegreeProbability(lastIndex - y)) {
				/*
				 * when both indices are the 'opposite', there is peak
				 * probability
				 */
				eArr[y][lastIndex - y] = qDist
						.getRemainingDegreeProbability(y);
				/*
				 * Other elements are zero, but no need to explicitly code that
				 * as arrays are automatically initialized to zero
				 */
			} else /* Need to maintain symmetry */
			{
				eArr[y][lastIndex - y] = qDist
						.getRemainingDegreeProbability(lastIndex - y);
				// eArr[y][y] = qDist.getRemainingDegreeProbability(y) -
				// qDist.getRemainingDegreeProbability(qDist.getLength() - y);

				double extrap = Math.abs(qDist.getRemainingDegreeProbability(y)- qDist.getRemainingDegreeProbability(lastIndex - y));
				double extraPPortion = extrap / ((double) (eArr.length - 2.0));

				for (int l = 0; l < eArr.length; l++) {
					double partnerExtraPortion;
					if (qDist.getRemainingDegreeProbability(l) <= qDist
							.getRemainingDegreeProbability(lastIndex - l)) {

						partnerExtraPortion = 0;
					} else {
						partnerExtraPortion = Math.abs(qDist.getRemainingDegreeProbability(l)- qDist.getRemainingDegreeProbability(lastIndex - l))
								/ ((double) (eArr.length - 2.0));
					}

					double minPortion = Math.min(extraPPortion,
							partnerExtraPortion);
					double diffPortion = Math.abs(extraPPortion
							- partnerExtraPortion);
					if (l != (lastIndex - y)) {
						if (y != l) {
							eArr[y][l] = minPortion;
							eArr[y][y] = eArr[y][y] + diffPortion;
						} else {
							eArr[y][y] = eArr[y][y] + extraPPortion; /*
																		 * Assortative
																		 * Peak,
																		 * just
																		 * take
																		 * your
																		 * share
																		 */
						}
					}
				}

			}
		}

		EDistribution ed = new EDistribution(eArr);
		// EDistribution ed = new EDistribution(fakeEARR);
		return ed;
	}

	/**
	 * Calculating the ejk for disassortativeness according to Neman formula 31
	 * in the PRE 67
	 * 
	 * @param q
	 *            q(k)
	 * @param x
	 *            x(k) -must have the same length as q.
	 * @return
	 */
	public static EDistribution generateNewmanEdist(double[] q, double[] x) {
		double[][] eArr = new double[q.length][q.length];
		for (int j = 0; j < eArr.length; j++) {
			for (int k = 0; k < eArr.length; k++) {
				eArr[j][k] = (q[j] * x[k]) + (q[k] * x[j]) - (x[j] * x[k]);
			}

		}
		EDistribution ed = new EDistribution(eArr);
		return ed;
	}

	/**
	 * This fucntion checks the validity of an Ejk distribution 1) All values
	 * are betwen 0 and +1.0 2) ejk = ekj 3) rowns and columns sum up to q(k)
	 * 
	 * @param ed
	 * @param q
	 * @return
	 */
	public static boolean checkEDistValidity(EDistribution ed, QDistribution_PowerLaw q) {
		double[][] eArr = ed.getDistribution();
		int counter = 0;
		boolean check = true;
		for (int j = 0; j < eArr.length; j++) {
			double rowsum = 0;
			for (int k = 0; k < eArr.length; k++) {
				if (eArr[j][k] < 0) {

					return false;

				}

				if (eArr[j][k] > 1) {
					System.out.println(" > 1");
					return false;
				}

				if (eArr[k][j] != eArr[j][k]) {
					System.out.println("non symmetry");
					return false;
				}

				rowsum = rowsum + eArr[j][k];

			}

			/*
			 * if(rowsum != q.getRemainingDegreeProbability(j)) { return false; }
			 */
		}

		return check;
	}

	/**
	 * This method generates the e distribution for the given q distribution and
	 * perfect dis-assortativess (r=-1) Note that if qDistribution is non
	 * symmetric, r=-1 is not possible. This method then generates an
	 * approximation of best (minusmost) minus r This methods' approximation is
	 * WORSE than the one before this.
	 * 
	 */
	public EDistribution generateDisassortativeEDistribution2(
			QDistribution_PowerLaw qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		for (int y = 1; y < qDist.getLength(); y++) /*
													 * Note the indexes starting
													 * from 1 rather than 0.
													 * Need to stay this way
													 */
		{
			if (qDist.getRemainingDegreeProbability(y) <= qDist
					.getRemainingDegreeProbability(qDist.getLength() - y)) {
				/*
				 * when both indices are the 'opposite', there is peak
				 * probability
				 */
				eArr[y][eArr.length - y] = qDist
						.getRemainingDegreeProbability(y);
				/*
				 * Other elements are zero, but no need to explicitly code that
				 * as arrays are automatically initialized to zero
				 */
			} else /* Need to maintain symmetry */
			{
				eArr[y][eArr.length - y] = qDist
						.getRemainingDegreeProbability(qDist.getLength() - y);
				eArr[y][y] = qDist.getRemainingDegreeProbability(y)
						- qDist.getRemainingDegreeProbability(qDist.getLength()
								- y);
				/* gain Other elements are zero */

			}
		}

		EDistribution ed = new EDistribution(eArr);
		return ed;
	}

	public EDistribution generateDisassortativeEDistribution3(
			QDistribution_PowerLaw qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];

		eArr[1][4] = (1.0 / 6.0);
		eArr[2][2] = (1.0 / 6.0);
		eArr[2][3] = (1.0 / 6.0);
		eArr[2][4] = (0 / 12.0);

		eArr[3][2] = (1.0 / 6.0);
		eArr[4][1] = (1.0 / 6.0);
		eArr[4][2] = (0 / 12.0);
		eArr[4][4] = (1.0 / 6.0);

		EDistribution ed = new EDistribution(eArr);
		return ed;
	}

	/**
	 * This method generates the e distribution for the given q distribution and
	 * complete non-assortativess (r=0)
	 * 
	 */
	public EDistribution generateNonassortativeEDistribution(QDistribution_PowerLaw qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		for (int j = 0; j < qDist.getLength(); j++) {
			for (int k = 0; k < qDist.getLength(); k++) {
				eArr[j][k] = qDist.getRemainingDegreeProbability(j)
						* qDist.getRemainingDegreeProbability(k);
			}
		}

		EDistribution ed = new EDistribution(eArr);
		return ed;
	}

	/**
	 * This method generates the e distribution for the given q distribution and
	 * complete non-assortativess (r=0). static version
	 * 
	 */
	public static EDistribution generateNonAssortativeEDistribution(
			QDistribution_PowerLaw qDist) {
		double[][] eArr = new double[qDist.getLength()][qDist.getLength()];
		for (int j = 0; j < qDist.getLength(); j++) {
			for (int k = 0; k < qDist.getLength(); k++) {
				eArr[j][k] = qDist.getRemainingDegreeProbability(j)
						* qDist.getRemainingDegreeProbability(k);
			}
		}

		EDistribution ed = new EDistribution(eArr);
		return ed;
	}

	JPanel resultPanel() {
		JPanel rPanel = new JPanel();
		// rPanel.setLayout(new GridLayout(1,3));
		JPanel centrePanel = new JPanel();
		JPanel leftPanel = new JPanel();
		JPanel rightPanel = new JPanel();

		prField = new JTextField(4);
		prField.setEditable(false);
		exrField = new JTextField(4);
		exrField.setEditable(true);
		maxMinusRField = new JTextField(4);
		maxMinusRField.setEditable(false);
		maxAvgKField = new JTextField(4);
		maxAvgKField.setEditable(false);
		rField = new JTextField(4);
		rField.setEditable(false);
		avgkField = new JTextField(4);
		avgkField.setEditable(false);
		nField = new JTextField(4);
		nField.setEditable(false);
		exavgkField = new JTextField(4);
		exavgkField.setEditable(true);
		exnField = new JTextField(4);
		exnField.setEditable(true);
		exGammaField = new JTextField(4);
		exGammaField.setEditable(true);
		buildNetButton = new JButton("Build Network");
		buildNetButton.addActionListener(new ButtonListener());
		miField = new JTextField(4);
		miField.setEditable(false);

		/* Input Fields */
		leftPanel.add(new JLabel("Expected Network Size"));
		leftPanel.add(exnField);
		leftPanel.add(new JLabel("Expected <k>"));
		leftPanel.add(exavgkField);
		leftPanel.add(new JLabel("expected r"));
		leftPanel.add(exrField);
		leftPanel.add(new JLabel("Expected Gamma"));
		leftPanel.add(exGammaField);

		centrePanel.add(buildNetButton);

		/* Display Fields */
		rightPanel.add(new JLabel("MAX MINUS R"));
		rightPanel.add(maxMinusRField);
		rightPanel.add(new JLabel("MAX AVG K"));
		rightPanel.add(maxAvgKField);
		rightPanel.add(new JLabel("r"));
		rightPanel.add(rField);
		rightPanel.add(new JLabel("Percentage Error"));
		rightPanel.add(prField);
		rightPanel.add(new JLabel("<k>"));
		rightPanel.add(avgkField);
		rightPanel.add(new JLabel("Network Size"));
		rightPanel.add(nField);
		rightPanel.add(new JLabel("I(q)"));
		rightPanel.add(miField);

		rPanel.add(leftPanel);
		rPanel.add(centrePanel);
		rPanel.add(rightPanel);

		rPanel.setPreferredSize(new java.awt.Dimension(1100, 100));
		return rPanel;
	}

	JPanel networkPanel() {
		JPanel networkPanel = new JPanel();

		indField = new JTextField(12);
		indField.setEditable(false);

		networkPanel.add(new JLabel("IND"));
		networkPanel.add(indField);

		networkPanel.setPreferredSize(new java.awt.Dimension(700, 200));
		return networkPanel;
	}

	/* JFREECHART FUNCTIONS */
	ChartPanel eDistPanel() {
		eDistData = create_update_EDistDataset(
				generateNonassortativeEDistribution(generateQ()), generateQ());
		eDistData = update_EDistDataset(
				generateNonassortativeEDistribution(generateQ()), generateQ());

		eDistChart = ChartFactory.createLineChart("E distribution",
				"(j,k) pair", "e(j,k)", eDistData, PlotOrientation.VERTICAL,
				true, // include legend
				true, // tooltips
				false // urls
				);

		CategoryPlot plot = (CategoryPlot) eDistChart.getPlot();

		eDistChart.setBackgroundPaint(Color.white);

		plot.setBackgroundPaint(Color.WHITE);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);
		plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));

		CategoryItemRenderer r = plot.getRenderer();
		if (r instanceof LineAndShapeRenderer) {
			LineAndShapeRenderer renderer = (LineAndShapeRenderer) r;
			renderer.setBaseShapesVisible(true);
			renderer.setBaseShapesFilled(true);
		}

		r.setSeriesPaint(3, Color.PINK);

		r.setSeriesStroke(0, new BasicStroke(4.0f, BasicStroke.CAP_ROUND,
				BasicStroke.JOIN_ROUND, 1.0f));

		r.setSeriesStroke(1,
				new BasicStroke(2.0f, BasicStroke.CAP_ROUND,
						BasicStroke.JOIN_ROUND, 1.0f,
						new float[] { 10.0f, 6.0f }, 0.0f));

		r
				.setSeriesStroke(2, new BasicStroke(2.0f,
						BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 1.0f,
						new float[] { 6.0f, 6.0f }, 0.0f));

		r
				.setSeriesStroke(3, new BasicStroke(2.0f,
						BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 1.0f,
						new float[] { 2.0f, 6.0f }, 0.0f));

		NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();

		cp = new ChartPanel(eDistChart);
		cp.setPreferredSize(new java.awt.Dimension(900, 400));

		return cp;
	}

	ChartPanel netEDistPanel() {
		netEDistData = create_update_NETEDistDataset(
				generateNonassortativeEDistribution(generateQ()), generateQ());

		netEDistChart = ChartFactory.createLineChart(
				"E distribution of grown network", "(j,k) pair", "e(j,k)",
				netEDistData, PlotOrientation.VERTICAL, true, // include
																// legend
				true, // tooltips
				false // urls
				);

		CategoryPlot plot = (CategoryPlot) netEDistChart.getPlot();

		netEDistChart.setBackgroundPaint(Color.white);

		plot.setBackgroundPaint(Color.WHITE);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);
		plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));

		CategoryItemRenderer r = plot.getRenderer();
		if (r instanceof LineAndShapeRenderer) {
			LineAndShapeRenderer renderer = (LineAndShapeRenderer) r;
			renderer.setBaseShapesVisible(true);
			renderer.setBaseShapesFilled(true);
		}

		cp3 = new ChartPanel(netEDistChart);
		cp3.setPreferredSize(new java.awt.Dimension(900, 200));

		return cp3;
	}

	/* JFREECHART FUNCTIONS */
	ChartPanel qDistPanel() {
		qDistData = create_update_QDistDataset(generateQ());

		qDistChart = ChartFactory.createLineChart("q distribution", "k",
				"q(k)", qDistData, PlotOrientation.VERTICAL, true, // include
																	// legend
				true, // tooltips
				false // urls
				);

		CategoryPlot plot = (CategoryPlot) qDistChart.getPlot();

		qDistChart.setBackgroundPaint(Color.white);

		plot.setBackgroundPaint(Color.WHITE);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);
		plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));

		CategoryItemRenderer r = plot.getRenderer();
		if (r instanceof LineAndShapeRenderer) {
			LineAndShapeRenderer renderer = (LineAndShapeRenderer) r;
			renderer.setBaseShapesVisible(true);
			renderer.setBaseShapesFilled(true);
		}

		r.setSeriesStroke(0, new BasicStroke(4.0f, BasicStroke.CAP_ROUND,
				BasicStroke.JOIN_ROUND, 1.0f));

		cp2 = new ChartPanel(qDistChart);
		cp2.setPreferredSize(new java.awt.Dimension(900, 200));

		return cp2;
	}

	ChartPanel netQDistPanel() {
		netQDistData = create_update_QDistDataset(generateQ());

		netQDistChart = ChartFactory.createLineChart("q distribution", "k",
				"q(k)", netQDistData, PlotOrientation.VERTICAL, true, // include
																		// legend
				true, // tooltips
				false // urls
				);

		CategoryPlot plot = (CategoryPlot) netQDistChart.getPlot();

		netQDistChart.setBackgroundPaint(Color.white);

		plot.setBackgroundPaint(Color.WHITE);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);
		plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));

		CategoryItemRenderer r = plot.getRenderer();
		if (r instanceof LineAndShapeRenderer) {
			LineAndShapeRenderer renderer = (LineAndShapeRenderer) r;
			renderer.setBaseShapesVisible(true);
			renderer.setBaseShapesFilled(true);
		}

		r.setSeriesStroke(0, new BasicStroke(4.0f, BasicStroke.CAP_ROUND,
				BasicStroke.JOIN_ROUND, 1.0f));

		cp4 = new ChartPanel(netQDistChart);
		cp4.setPreferredSize(new java.awt.Dimension(900, 150));

		return cp4;
	}

	private DefaultCategoryDataset create_update_EDistDataset(
			EDistribution eDist, QDistribution_PowerLaw qDist) {
		double[][] eArr = eDist.getDistribution();
		double[][] eArrA = generateAssortativeEDistribution(qDist)
				.getDistribution(); // ASS
		double[][] eArrN = generateNonassortativeEDistribution(qDist)
				.getDistribution(); // NON
		double[][] eArrD = generateDisassortativeEDistribution(qDist)
				.getDistribution(); // DIS

		for (int i = 0; i < eArr.length; i++) {
			for (int j = 0; j < eArr[i].length; j++) {
				eDistData.addValue(eArr[i][j], series1, "" + i + "," + j + "");
				eDistData.addValue(eArrA[i][j], series2, "" + i + "," + j + "");
				eDistData.addValue(eArrD[i][j], series3, "" + i + "," + j + "");
				eDistData.addValue(eArrN[i][j], series4, "" + i + "," + j + "");
			}
		}

		return eDistData;
	}

	private DefaultCategoryDataset create_update_NETEDistDataset(
			EDistribution eDist, QDistribution_PowerLaw qDist) {
		double[][] eArr = eDist.getDistribution();

		for (int i = 0; i < eArr.length; i++) {
			for (int j = 0; j < eArr[i].length; j++) {
				netEDistData.addValue(eArr[i][j], netseries1, "" + i + "," + j
						+ "");
				netEDistData.addValue(eArr[i][j], netseries2, "" + i + "," + j
						+ "");

			}
		}

		return netEDistData;
	}

	private DefaultCategoryDataset create_update_QDistDataset(
			QDistribution_PowerLaw qDist) {

		for (int i = 0; i < qDist.getLength(); i++) {

			qDistData.addValue(qDist.getRemainingDegreeProbability(i),
					netseries3, "" + i + "");
			qDistData.addValue(qDist.getRemainingDegreeProbability(i),
					netseries4, "" + i + "");

		}

		return qDistData;
	}

	private DefaultCategoryDataset update_EDistDataset(EDistribution eDist,
			QDistribution_PowerLaw qDist) {
		if (qDist.getLength() > 10) /*
									 * Impossible to visualize with this screen
									 * size and resolution
									 */
		{
			eDistData.clear();
			return eDistData;
		}
		double[][] eArr = eDist.getDistribution();
		double[][] eArrA = generateAssortativeEDistribution(qDist)
				.getDistribution(); // ASS
		double[][] eArrN = generateNonassortativeEDistribution(qDist)
				.getDistribution(); // NON
		double[][] eArrD = generateDisassortativeEDistribution(qDist)
				.getDistribution(); // DIS

		int presentCols = (int) Math.sqrt(eDistData.getColumnCount());

		int commonCols = Math.min(presentCols, eArr.length);
		int diffCols = presentCols - eArr.length; /*
													 * Minus if incoming data is
													 * biggger in size
													 */

		for (int i = 0; i < (commonCols + Math.abs(diffCols)); i++) {
			if (i < commonCols) {
				for (int j = 0; j < (commonCols + Math.abs(diffCols)); j++) {
					if (j < commonCols) {

						try {
							eDistData.removeValue(series1, "" + i + "," + j
									+ "");
						} catch (Exception e) {
							// If it cant be removed its not there, dont worry
							// :-)
							// e.printStackTrace();
						}
						eDistData.addValue(eArr[i][j], series1, "" + i + ","
								+ j + "");
						try {
							eDistData.removeValue(series2, "" + i + "," + j
									+ "");
						} catch (Exception e) {
							// If it cant be removed its not there, dont worry
							// :-)
						}
						eDistData.addValue(eArrA[i][j], series2, "" + i + ","
								+ j + "");
						try {
							eDistData.removeValue(series3, "" + i + "," + j
									+ "");
						} catch (Exception e) {
							// If it cant be removed its not there, dont worry
							// :-)
						}
						eDistData.addValue(eArrD[i][j], series3, "" + i + ","
								+ j + "");
						try {
							eDistData.removeValue(series4, "" + i + "," + j
									+ "");
						} catch (Exception e) {
							// If it cant be removed its not there, dont worry
							// :-)
						}
						eDistData.addValue(eArrN[i][j], series4, "" + i + ","
								+ j + "");

					}

					else {
						if (diffCols >= 0) /* Existing Data size is bigger */
						{
							eDistData.removeValue(series1, "" + i + "," + j
									+ "");
							eDistData.removeValue(series2, "" + i + "," + j
									+ "");
							eDistData.removeValue(series3, "" + i + "," + j
									+ "");
							eDistData.removeValue(series4, "" + i + "," + j
									+ "");
						} else /* Incoming Data size is bigger */
						{
							eDistData.addValue(eArr[i][j], series1, "" + i
									+ "," + j + "");
							eDistData.addValue(eArrA[i][j], series2, "" + i
									+ "," + j + "");
							eDistData.addValue(eArrD[i][j], series3, "" + i
									+ "," + j + "");
							eDistData.addValue(eArrN[i][j], series4, "" + i
									+ "," + j + "");
						}
					}

				}
			} else {
				if (diffCols >= 0) /* Existing Data size is bigger */
				{
					for (int j = 0; j < (commonCols + Math.abs(diffCols)); j++) {
						eDistData.removeValue(series1, "" + i + "," + j + "");
						eDistData.removeValue(series2, "" + i + "," + j + "");
						eDistData.removeValue(series3, "" + i + "," + j + "");
						eDistData.removeValue(series4, "" + i + "," + j + "");
					}
				} else /* Incoming Data size is bigger */
				{
					for (int j = 0; j < (commonCols + Math.abs(diffCols)); j++) {
						eDistData.addValue(eArr[i][j], series1, "" + i + ","
								+ j + "");
						eDistData.addValue(eArrA[i][j], series2, "" + i + ","
								+ j + "");
						eDistData.addValue(eArrD[i][j], series3, "" + i + ","
								+ j + "");
						eDistData.addValue(eArrN[i][j], series4, "" + i + ","
								+ j + "");
					}
				}
			}
		}

		/* Now, remove the values and add them again, to sort the order */
		presentCols = (int) Math.sqrt(eDistData.getColumnCount());
		for (int i = 0; i < presentCols; i++) {
			for (int j = 0; j < presentCols; j++) {
				eDistData.removeValue(series1, "" + i + "," + j + "");
				eDistData.removeValue(series2, "" + i + "," + j + "");
				eDistData.removeValue(series3, "" + i + "," + j + "");
				eDistData.removeValue(series4, "" + i + "," + j + "");
				eDistData.addValue(eArr[i][j], series1, "" + i + "," + j + "");
				eDistData.addValue(eArrA[i][j], series2, "" + i + "," + j + "");
				eDistData.addValue(eArrD[i][j], series3, "" + i + "," + j + "");
				eDistData.addValue(eArrN[i][j], series4, "" + i + "," + j + "");
			}
		}

		return eDistData;
	}

	private DefaultCategoryDataset update_NETEDistDataset(EDistribution eDist,
			double[][] netEDist) {
		if (netEDist.length > 10) /*
									 * Impossible to visualize with this screen
									 * size and resolution
									 */
		{
			netEDistData.clear();
			return netEDistData;
		}

		double[][] eArr = eDist.getDistribution();
		int presentCols = (int) Math.sqrt(netEDistData.getColumnCount());

		int commonCols = Math.min(presentCols, eArr.length);
		int diffCols = presentCols - eArr.length; /*
													 * Minus if incoming data is
													 * biggger in size
													 */

		for (int i = 0; i < (commonCols + Math.abs(diffCols)); i++) {
			if (i < commonCols) {
				for (int j = 0; j < (commonCols + Math.abs(diffCols)); j++) {
					if (j < commonCols) {

						try {
							netEDistData.removeValue(netseries1, "" + i + ","
									+ j + "");
						} catch (Exception e) {
							// If it cant be removed its not there, dont worry
							// :-)
							// e.printStackTrace();
						}
						netEDistData.addValue(eArr[i][j], netseries1, "" + i
								+ "," + j + "");
						try {
							netEDistData.removeValue(netseries2, "" + i + ","
									+ j + "");
						} catch (Exception e) {
							// If it cant be removed its not there, dont worry
							// :-)
						}
						netEDistData.addValue(netEDist[i][j], netseries2, ""
								+ i + "," + j + "");

					}

					else {
						if (diffCols >= 0) /* Existing Data size is bigger */
						{
							netEDistData.removeValue(netseries1, "" + i + ","
									+ j + "");
							netEDistData.removeValue(netseries2, "" + i + ","
									+ j + "");

						} else /* Incoming Data size is bigger */
						{
							netEDistData.addValue(eArr[i][j], netseries1, ""
									+ i + "," + j + "");
							netEDistData.addValue(netEDist[i][j], netseries2,
									"" + i + "," + j + "");

						}
					}

				}
			} else {
				if (diffCols >= 0) /* Existing Data size is bigger */
				{
					for (int j = 0; j < (commonCols + Math.abs(diffCols)); j++) {
						netEDistData.removeValue(netseries1, "" + i + "," + j
								+ "");
						netEDistData.removeValue(netseries2, "" + i + "," + j
								+ "");

					}
				} else /* Incoming Data size is bigger */
				{
					for (int j = 0; j < (commonCols + Math.abs(diffCols)); j++) {
						netEDistData.addValue(eArr[i][j], netseries1, "" + i
								+ "," + j + "");
						netEDistData.addValue(netEDist[i][j], netseries2, ""
								+ i + "," + j + "");

					}
				}
			}
		}

		/* Now, remove the values and add them again, to sort the order */
		presentCols = (int) Math.sqrt(netEDistData.getColumnCount());
		for (int i = 0; i < presentCols; i++) {
			for (int j = 0; j < presentCols; j++) {
				netEDistData.removeValue(netseries1, "" + i + "," + j + "");
				netEDistData.removeValue(netseries2, "" + i + "," + j + "");
				netEDistData.addValue(eArr[i][j], netseries1, "" + i + "," + j
						+ "");
				netEDistData.addValue(netEDist[i][j], netseries2, "" + i + ","
						+ j + "");
			}
		}

		return netEDistData;
	}

	private DefaultCategoryDataset update_QDistDataset(QDistribution_PowerLaw qDist,
			double[] netq) {

		/*
		 * if(qDist.getLength()>100) { qDistData.clear(); return qDistData; }
		 */

		int presentCols = qDistData.getColumnCount();

		int commonCols = Math.min(presentCols, netq.length);
		int diffCols = presentCols - netq.length; /*
													 * Minus if incoming data is
													 * biggger in size
													 */

		if (qDist.getLength() > 100) /*
										 * Impossible to visualize with this
										 * screen size and resolution, limit to
										 * 100 points
										 */
		{
			commonCols = Math.min(presentCols, 100);
			diffCols = presentCols - 100;
		}

		for (int i = 0; i < (commonCols + Math.abs(diffCols)); i++) {
			if (i < commonCols) {
				try {
					qDistData.removeValue(netseries3, "" + i + "");
				} catch (Exception e) {
					// e.printStackTrace();
				}
				qDistData.addValue(qDist.getRemainingDegreeProbability(i),
						netseries3, "" + i + "");
				try {
					qDistData.removeValue(netseries4, "" + i + "");
				} catch (Exception e) {
					// e.printStackTrace();
				}
				qDistData.addValue(netq[i], netseries4, "" + i + "");

			}

			else {
				if (diffCols >= 0) /* Existing Data size is bigger */
				{
					qDistData.removeValue(netseries3, "" + i + "");
					qDistData.removeValue(netseries4, "" + i + "");
				} else /* Incoming Data size is bigger */
				{
					qDistData.addValue(netq[i], netseries3, "" + i + "");
					qDistData.addValue(netq[i], netseries4, "" + i + "");
				}
			}

		}

		return qDistData;
	}

	/* PRIVATE CLASSES needed for the frame */

	private class ButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent Event) {

			JButton theItem = (JButton) Event.getSource();

			// Call the building function
			build();

		}
	}

	/* End of class */
}
