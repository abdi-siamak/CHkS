package uk.ac.rdg.chks.bb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.*;

import uk.ac.rdg.chks.bb.BnB_Search.SearchType;
import uk.ac.rdg.chks.graph.Edge;
import uk.ac.rdg.chks.graph.TermGraph;
import uk.ac.rdg.chks.graph.Vertex;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
/**
 * Main class for Connected Heaviest k-Subgraph (CHkS)<br>
 * Branch and bound search strategy: DFS, BFS, Hybrid
 * 
 * @author Giuseppe Di Fatta and Salwa Musa
 *
 */
public class CHkS {
	public static HashMap<Integer, String> dictionary = new HashMap<>();
	private static final String DATADIR = "data";
	private static final char SEP = '/';
	private static DecimalFormat df = new DecimalFormat("#.##");

	public static int numberOfEdges = 0;
	public static int maxID = 0;
	public static int k;
	
	/** 
	 * The graph file is parsed to find the vertex with the maximum ID.
	 * An empty matrix with the number of cells equal to this maximum ID is initialised.
	 */
	public static int getMaxID(String filename) {
		int id = 0;
		String line = "";
		BufferedReader br = null;
		StringTokenizer st;
		File f = null;
		int i = 0;

		try {
			f = new File(filename);
//			System.out.print(f);
		} catch (Exception e) {
			e.printStackTrace();
		}
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(f)));
		} catch (FileNotFoundException e) {
			System.err.println("Error opening file here! FileNotFound: "+f);
			System.exit(-1);;
		}

		try {
			line = br.readLine();
			while (line != null) {
//					System.out.println("line: "+line);
					if ((line != null) && (line.trim().length()>0)){
						numberOfEdges++;
						st = new StringTokenizer(line);
						i = Integer.parseInt(st.nextToken());
						if (i > id)
							id = i;
						i = Integer.parseInt(st.nextToken());
						if (i > id)
							id = i;
					}
					line = br.readLine();
			}
		} catch (IOException e) {
				System.err.println("Error reading line .");
				e.printStackTrace();
				System.exit(-1);;
		}

		maxID = id;
		return maxID;
	}
	
	/**
	 *  this methods will return a TermGraph object with the vertices and edges.
	 * @param filename
	 * @return TermGraph
	 */
	public static TermGraph loadGraph(String filename) {
		//System.out.println("Loading the graph from file <"+filename+">...");
		
		String line = "";
		BufferedReader br = null;
		StringTokenizer st;
		File f = null;
		int source = 0, destination = 0, weight = 0;
		Edge newEdge = null;
		int length = getMaxID(filename) + 1; // as the vertices id will start by 0 for this reason we have to add 1
		//System.out.println("length" + length);

		ArrayList<Vertex> vertices = new ArrayList<Vertex>(length); //[length];
		for (int i = 0; i < length; i++) {
//			vertices[i] = new Vertex(i);
			vertices.add(new Vertex(i));
		}
		
		//Only for debugging: remove
//		for (int i = 0; i < length; i++) {
//			if(vertices.get(i).getID() != i) {
//				System.out.println("ERROR: ID " + vertices.get(i).getID()+ "!= i " + i);
//				System.exit(-1);
//			}
//		}
		
		ArrayList<Edge> edges = new ArrayList<Edge>();

		try {
			f = new File(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(f)));
		} catch (FileNotFoundException e) {
			System.err.println("Error opening file for writing!");
		}

		// ArrayList<Edge> edges = new ArrayList<Edge>();
		try {
			line = br.readLine();
			while (line != null) {
					if (line.trim().length()>0) {
						st = new StringTokenizer(line);
						source = Integer.parseInt(st.nextToken());
						destination = Integer.parseInt(st.nextToken());
						weight = Integer.parseInt(st.nextToken());
						vertices.get(source).addNeighbour(vertices.get(destination), weight);
						vertices.get(destination).addNeighbour(vertices.get(source), weight);
	
						newEdge = new Edge(vertices.get(source), vertices.get(destination), weight);
						edges.add(newEdge);
						vertices.get(source).addEdge(newEdge);
						vertices.get(destination).addEdge(newEdge);
					}
					line = br.readLine();
			}
			br.close();
		} catch (IOException e) {
			System.err.println("Error reading line .");
			e.printStackTrace();
			System.exit(-1);;
		}


		// information about the Graph and the file 
//		System.out.println("\ninformation about the Graph and the file");
//		System.out.println("Graph Path :" + filename);
		System.out.println("number of vertices : " + vertices.size());
		//System.out.println("***********MaxID==" + getMaxID(filename));
		System.out.println("number of edges : " + edges.size());
//		System.out.println("********************************************************");
//		System.out.println("Finished loading file " + filename);
		
		return new TermGraph(vertices, edges);
	}
	public static void ExcelToHashMap(String path) {
		try {
			FileInputStream file = new FileInputStream(path);
			Workbook workbook = new XSSFWorkbook(file);
			Sheet sheet = workbook.getSheetAt(0); // assuming data is in the first sheet

			for (Row row : sheet) {
				Cell keyCell = row.getCell(0);
				Cell valueCell = row.getCell(1);
				Integer key = (int) keyCell.getNumericCellValue();
				String value = valueCell.getStringCellValue();
				dictionary.put(key, value);
			}
			//workbook.close();
			file.close();
			// use the hashmap as needed
			//System.out.println(dataMap);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	/**
	 * Starting method for CHkS:
	 * mandatory parameter: search method (DFS, BFS, Hybrid
	 * optional parameter: file name for the term graph. If no file is provided one of the test graphs in the dir "data" is used.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		if((args.length != 1) && (args.length != 2)) {
			System.out.println("usage: java CHkS searchMethod (DFS|BFS|HYBRID) [graphfile.txt]\nExit.");
			System.exit(-1);
		}
	    SearchType searchMethod = BnB_Search.searchMethod(args[0]);
		if(searchMethod==null) {
			System.out.println("Search method is null.");
			System.out.println("usage: java CHkS searchMethod (DFS|BFS|HYBRID) [graphfile.txt]\nExit.");
			System.exit(-1);
		}
		
		//Test graph
//		String filename="simplegraph.txt"; //a simple test graph
		String filename="graph.txt"; //a test graph
		ExcelToHashMap("dictionary.xlsx"); // load the dictionary file of the graph
		//System.out.println(dictionary);
		//Internet-like topologies (brite)
//		String filename="brite10K.txt";
//		String filename="brite0-1K-2K.txt";
		Top_Solutions.topSize = 10; // the number of top best solutions
		k = 6; //the number of nodes in the subgraph
		Top_Solutions.setSimilarity(0.5); // set the value of similarity measure
//		int maxedges = k * (k - 1) / 2;
//		int lowerBound = 0;

		//Random graphs
//		String filename="rand40K.txt";
//		String filename="rand30K.txt";
//		String filename="rand20K.txt";
//		String filename="rand10K.txt";
//
//		String filename="rand10K15.txt";
//		String filename="rand7K15.txt";
//		String filename="rand5K.txt";
//		String filename="rand3K.txt";  //3000 nodes, 45000 edges: DFS (k=5): 30 minutes

		//Twitter data		
//		String filename="TermGraph30.txt"; //DFS (k=5) <15 sec
//		String filename="TermGraph20.txt"; //DFS (k=5) <10 sec
//		String filename="TermGraph10.txt"; //DFS (k=5) <10 sec

		if((args.length == 2) && (args[1] != null) && (args[1].length() > 0)){
			filename = args[1]; //optional input filename overrides the above assignment to the filename
			File f = new File(DATADIR+SEP+filename);
			if(!f.exists()) {
				System.out.println("Error: file <"+f.getAbsolutePath()+"> not found.\nExit.");
				System.exit(-1);
			}
		}

		System.out.println("##### CHkS DFS (k=" + k + "): " + filename + " #####");

		long startTime = System.currentTimeMillis();
		long lastTime = startTime;
		TermGraph theTermGraph = loadGraph(DATADIR+SEP+filename);

		long currentTime = System.currentTimeMillis();
		long loadingFileTime = currentTime - lastTime;
		lastTime = currentTime;
//		double fileTime = (double) loadingFileTime / 1000;
		System.out.println("Time for loading the graph from file : " + df.format(loadingFileTime/1000.0) + " seconds");

//		System.out.println("number of Vertices: " + tg.getNumberOfVertices());
//		System.out.println("number of Edges: " + tg.getTotalEdges());
//		System.out.println("The max number of vertices is: " + k);

//		System.out.println("--------------------------------------------");
		System.out.println("Ordering the vertices by weighted degree...");
		//TODO: also try to sort them by the sum of the top k-1 edges, and compare performance.
		//tg.sort_Vertice();// using Map sorting
		theTermGraph.sortVertices(); // Using collection
//		System.out.println("Finished ordering the vertices:");
		theTermGraph.print_sortedVertices(20);
		// tg.Print_EV();
//		System.out.println("--------------------------------------------");
		System.out.println("Ordering the edges...");
		//tg.sortEdges();   // Using  Heap
		theTermGraph.sort_Edges(); // Using Collection
//		System.out.println("Finished ordering the edges:");
		theTermGraph.printEdges(20);
//		System.out.println("--------------------------------------------");
//		System.exit(0);
		
		
		if(BnB_Node.OPT3) { //some preprocessing required for opt3
			System.out.println("Computing nearest neighbours...");
			Vertex.compute_kNN_for_all_vertices(theTermGraph, k);
		}
//		tg.getVertex_byID(8).compute_kNN(k);
//		tg.getVertex_byID(8).print_kNN();
//		System.out.println("--------------------------------------------");
		
		currentTime = System.currentTimeMillis();
		double preprocTime = (currentTime - lastTime)/1000.0;
		lastTime = currentTime;
		System.out.println("Total preprocessing time (excl. loading file): "+ df.format(preprocTime)+ " seconds");
		
		
		// Branch and Bound
		System.out.println("Starting Branch and Bound... ("+searchMethod+")\n");
		long startAlgorithm = System.currentTimeMillis();

		// BnBNode n=branchAB.bbHeap(k, lowerBound, tg);
		//branchAB.bbHeap(k, lowerBound, tg);
		//branchAB.bbHeap(k, tg);

		//BaB_Search.SearchType method = BaB_Search.SearchType.BFS;
		BnB_Search search = new BnB_Search(theTermGraph, k);
		search.start_search(searchMethod);

		System.out.printf("\n");
		int index = 1;
		Set<Map.Entry<BnB_Node, List<BnB_Node>>> entrySet = Top_Solutions.getChildSolutions().entrySet();
		List<Map.Entry<BnB_Node, List<BnB_Node>>> list = new ArrayList<>(entrySet);
		ListIterator<Map.Entry<BnB_Node, List<BnB_Node>>> iterator = list.listIterator(list.size());
		while (iterator.hasPrevious()){
			System.out.println("##### Children of solution: " + index+ " #####");
			Map.Entry<BnB_Node, List<BnB_Node>> entry = iterator.previous();
			int child = 1;
			for (BnB_Node childSolution: entry.getValue()){
				System.out.println("child: " + child );
				StringBuffer buf = new StringBuffer();
				for (Vertex v : childSolution.getVertices()) {
					buf.append(CHkS.dictionary.get(v.getID())+" ");
				}
				System.out.println(buf + "; "+childSolution.getWeight());
				child++;
			}
			System.out.println("-----------------------------------------");
			index++;
		}

		long totalBnB_time = System.currentTimeMillis() - startAlgorithm;
		double algTime = totalBnB_time / 1000.0;
		System.out.println("\nAlgorithm total time : " + df.format(algTime) + " seconds");
		System.out.println("(in UB computations: opt1 = "+BnB_Node.opt1Cnt+", opt2 = "+BnB_Node.opt2Cnt +
				", opt3 = "+BnB_Node.opt3Cnt + ")"); //+", opt3b = "+BnBNode.opt3bCnt+")");
		System.out.println("Total number of pruned nodes - for UB: "+BnB_Search.prunedDueToUpperBound + 
				" and for k>="+k+": "+ BnB_Search.prunedDueToVertices);
		
		System.out.println();
		System.out.println("ver\tfilename\t|V|\t|E|\tk\tmethod\tbestCost\tbestSols\tpreprocTime(sec)\ttime(sec)\tNumOfNode\tNumVisited\topt1Cnt\topt2\topt2Cnt\topt3\topt3Cnt\tHsize");
		System.out.println(BnB_Search.ver + "\t" + filename + "\t" + theTermGraph.getNumberOfVertices() + "\t" + theTermGraph.getTotalEdges() + "\t" + k + 
				"\t" + searchMethod + "\t" + search.m_bestCost + "\t" + search.bestSolutions.size() +
				"\t" + preprocTime + "\t" + algTime + 
				"\t" + BnB_Node.getTotNumOfNode() + "\t" + BnB_Search.visitedNodes + 
				"\t" + BnB_Node.opt1Cnt + "\t" + BnB_Node.OPT2 + "\t" + BnB_Node.opt2Cnt + "\t" + BnB_Node.OPT3+ "\t" + BnB_Node.opt3Cnt +
				"\t" + BnB_Search.maxHeap_size); 
//				+ "\t"+BranchAndBound.betterUBcnt);
		
		System.out.println("\nEnd.");
	}
}
