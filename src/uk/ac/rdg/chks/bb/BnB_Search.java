package uk.ac.rdg.chks.bb;

import java.util.*;

import uk.ac.rdg.chks.graph.TermGraph;
import uk.ac.rdg.chks.graph.Vertex;
import uk.ac.rdg.chks.graph.VertexComparator_byWeight_reverse;

/**
 * This class is responsible to create Branch and Bound nodes to visit all possible solutions.
 * It implements the search strategy algorithms: DFS, BFS, HYBRID.
 * 
 * @author Giuseppe Di Fatta and Salwa Musa
 * 
 */
public class BnB_Search {
	private static boolean verb = false;
	public static String ver = "v.1.1";
	
	public enum SearchType {
	    DFS,
	    BFS,
	    HYBRID
	}

	public static int maxHeap_size=0; //for stats: the maximum size of the heap over the entire execution
	private static final int MAXHEAPSIZE = 1000; //for the hybrid approach: the parameter to switch between BFS and DFS
    
	TermGraph theTG; //the term graph: user input
	int m_maxVertices; //k: user input
	int m_maxEdges; //the theoretical max number of edges for a given k
	
	double m_bestCost ; // double bestLowerBound; // bestLowerBound contains the highest weight among the subgraphs found
	Set<BnB_Node> bestSolutions = new HashSet<BnB_Node>(); // bestSolution which is the best subgraph we have obtained, which has to have the heaviest weight
	static int visitedNodes=0;
//	static int betterUBcnt=0;
	
	static int prunedDueToUpperBound = 0;
	static int prunedDueToVertices = 0;

	
	public BnB_Search(TermGraph tg, int k) {
		theTG = tg;
		m_maxVertices = k;
	    m_bestCost = 0;
	    m_maxEdges = k * (k - 1) / 2;
	    if (m_maxEdges > tg.getTotalEdges())
		    m_maxEdges = tg.getTotalEdges();
	}

    public static SearchType searchMethod(String value){
        if(value.equalsIgnoreCase(SearchType.DFS.toString()))
            return SearchType.DFS;
        else if(value.equalsIgnoreCase(SearchType.BFS.toString()))
            return SearchType.BFS;
        else if(value.equalsIgnoreCase(SearchType.HYBRID.toString()))
            return SearchType.HYBRID;
        else
            return null;
    }

	/**
	 * start the BnB search for each root node (for each vertex in order of weighted degree)
	 */
	void start_search(SearchType type) { //int maxVertices, TermGraph termG) {
		BnB_Node rootNode;
		
		List<Vertex> vertices_List = theTG.getSortingVertices(); // List of Sorting Vertices
		//Set<Vertex> sx1 = new HashSet<Vertex>();

//		boolean visited[] = new boolean[vertices_List.size()];

		ArrayList<Integer> ncnt = new ArrayList<Integer>();
		int nodeCnt, prevNodeCnt=0;
		Set<Vertex> sr = new HashSet<Vertex>(); //set of excluded root vertices for the sibling trees
		for (int i = 0; i < vertices_List.size(); i++) {
			Vertex rootVertex = vertices_List.get(i);
//			if(verb)
//				System.out.println("Computing BnB search tree N."+i+" rooted at vertex <"+rootVertex.getID()+">...");

			rootNode = new BnB_Node(this, rootVertex); //(maxVertices, termG, rootVertex);
			
			//System.out.println(rootNode);
			//System.out.println(rootVertex.getID());
			
			rootNode.excludeVertices(sr);
			sr.add(rootVertex);
			
			//start visiting the search tree rooted at this root node
			if(verb)
				System.out.println(rootNode);

	    	Comparator<BnB_Node> bnBNodeComparator = new Comparator<BnB_Node>() {
	            @Override
	            public int compare(BnB_Node n1, BnB_Node n2) { //decreasing order of the upper bound
//	                return n2.upperBound - n1.upperBound;
	            	if (n1.upperBound > n2.upperBound)
	        			return -1;
	        		else if (n1.upperBound == n2.upperBound)
	        			return 0;
	        		else
	        			return 1;
	            }
	        };
	        
			switch (type) {
		    case DFS: //using recursive calls
		    	dfs_search(rootNode); 
		    	break;
		    case BFS: //using priority queue (heap)
//		    	BnBNodesHeap bfs_heap = new BnBNodesHeap(100);
		    	PriorityQueue<BnB_Node> bfs_heap = new PriorityQueue<BnB_Node>(bnBNodeComparator);
		    	bfs_heap.add(rootNode);
		    	bfs_search(bfs_heap);
		    	break;
		    case HYBRID:
		    	PriorityQueue<BnB_Node> hfs_heap = new PriorityQueue<BnB_Node>(bnBNodeComparator);
		    	hfs_heap.add(rootNode);
		    	hfs_search(hfs_heap);
				break;
		    default:
				System.out.println("BaB search method NA. Quit.");
				System.exit(-1);
			}
			 
			nodeCnt = BnB_Node.getTotNumOfNode() - prevNodeCnt;
			ncnt.add(nodeCnt);
			prevNodeCnt = BnB_Node.getTotNumOfNode();
			if(verb)
				System.out.println("End of search tree N."+i+" rooted at vertex <"+rootVertex.getID()+">: nodeCnt="+nodeCnt + " (total="+BnB_Node.getTotNumOfNode()+ ")");
			if((i<=1) || (i%100 == 0)){/*
				System.out.println("Completed BnB search tree N."+i+"...  \t(dfs_search_calls: " + visitedNodes
						+ " and nodes in this BnB search tree: " + nodeCnt + "   and total number of nodes: " + BnB_Node.getTotNumOfNode()+")");
						*/
			}

		}
		System.out.println("completed >>> dfs_search_calls: "+visitedNodes+" and BnB nodes: "+BnB_Node.getTotNumOfNode());

		
		System.out.println();
		//System.out.println("----- Number of BnB nodes in each BnB tree:");
		/*
		int max = (ncnt.size()<20)?ncnt.size():20;
		for(int i=0; i<max; i++) {
			System.out.println("nodeCnt["+i+"] = "+ncnt.get(i));
		}
		System.out.println("...");
		System.out.println("-------");
		System.out.println("total="+BnB_Node.getTotNumOfNode()+ "\n");
		*/
		System.out.println("\n##### Number of best solutions found is " + bestSolutions.size() + ": #####");
		int i = 1;
		for (BnB_Node nsb : bestSolutions) {
			System.out.println("----------- best solution N."+i+" -----------");
			System.out.println(nsb);
			System.out.println("cost: "+ nsb.getWeight());
//			System.out.println("Upper bound for the Best Solution: " + nsb.getUpperBound() + "\n");
//			System.out.println("Lower bound for the Best Solution: " + nsb.getWeight() + "\n");
//			System.out.println("The solution in the level: " + nsb.getTreeLevel() + "\n");
//			System.out.println("The Vertices in the best solution are: " + nsb.getVertices() + "\n");
//			System.out.println("The Edges in the best solution are: " + "\n");
//			nsb.getEdges();
			
//			System.out.println("-----------------------------------------------");
			i = i + 1;
		}
		////////////////////////////////////////////////////////////////////////////////////////////
		System.out.println("\n##### Top " + Top_Solutions.topSize + " solutions: #####");;
		int t = 1;
		for (BnB_Node nsb:Top_Solutions.getSolutions()) {
			System.out.println("----------- "+t+ " best solution ---------------");
			System.out.println(nsb);
			System.out.println("cost: "+ nsb.getWeight());
//			System.out.println("Upper bound for the Best Solution: " + nsb.getUpperBound() + "\n");
//			System.out.println("Lower bound for the Best Solution: " + nsb.getWeight() + "\n");
//			System.out.println("The solution in the level: " + nsb.getTreeLevel() + "\n");
//			System.out.println("The Vertices in the best solution are: " + nsb.getVertices() + "\n");
//			System.out.println("The Edges in the best solution are: " + "\n");
//			nsb.getEdges();
			t = t + 1;
			//System.out.println("-----------------------------------------------");
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * DFS using recursive implementation<br>
	 * Node p is the parent node from which to generate child BnB nodes (siblings).<br>
	 * 
	 * @param p
	 */
	public void dfs_search(BnB_Node p) {
		visitedNodes++;
		//System.out.println("<dfs_search> visitedNodes: " + visitedNodes);

// 		Set<Vertex> candidateSet = new HashSet<Vertex>(); // candidate vertices for spawning children
		List<Vertex> candidatesList = new ArrayList<Vertex>(); // candidate vertices for spawning children

		if(verb) 
			System.out.println("<dfs_search> The neighbours of the last vertex: <" + p.last_v + ">" + "are:" + p.last_v.getNeighbourVertices());

		//1. add the candidates inherited from the previous parent node:
		for (Vertex c: p.candidate_vertices) {
			if ((!(p.vertices.contains(c)) && (!(p.excluded_vertices.contains(c)))))
				if(!candidatesList.contains(c))
					candidatesList.add(c);
		}
	
		// 2. add new candidates from the neighbours of the last added vertex
		for (Vertex vn : p.last_v.getNeighbourVertices()) {
			if ((!(p.vertices.contains(vn)) && (!(p.excluded_vertices.contains(vn))))) {
//				candidateSet.add(vn);
				if(!candidatesList.contains(vn))
					candidatesList.add(vn);
			}
		}
		
		// Sorting the candidates set based on weighted degree
//		candidatesList.addAll(candidateSet);
//		Collections.sort(candidatesList, Collections.reverseOrder()); //??
		Collections.sort(candidatesList, new VertexComparator_byWeight_reverse());
//		if(verb) 
//		if(p.last_v.getID() == 10)
//			System.out.println("<dfs_search> The sorted candidates are: " + candidatesList+"\n");

		p.setCandidateVertices(candidatesList); //set the updated candidates list into this parent node before generating child nodes
//		if(verb) 
//			System.out.println("<dfs_search> The sorted candidates at the parent are: " + p.candidate_vertices+"\n");

		Set<Vertex> sc = new HashSet<Vertex>(); //set of excluded candidate vertices for the sibling nodes
		for (Vertex cv: candidatesList) { // generate a child node for each candidate vertex
			BnB_Node child = new BnB_Node(p); //generate a child node: a copy of the parent node
//			if(verb)
//				System.out.println("<dfs_search> child with candidate <"+cv.getID()+">:\n"+child);
			
			child.excludeVertices(sc); //BUG found: this was done after expand, which also used to compute the UB
			child.expandSolution(cv); // adds a new vertex and its edges to the current solution
			child.computeUpperBound();
			
//			if(verb) {
//				System.out.println("<dfs_search> child with candidate <"+cv.getID()+">:\n"+child);
//				System.out.println("<dfs_search> the parent candidates were: "+p.candidate_vertices);
//				System.out.println("<dfs_search> the child candidates are:   "+child.candidate_vertices);
//			}
			
//			double ub2 = child.computeUB(true, false);
//			double ub3 = child.computeUB(true, true);
//			if(BnBNode.OPT2 && (!BnBNode.OPT3)) {
//				if(ub2 != child.upperBound) {
//					//error:
//					System.out.println("Error: ub2 != child.upperBound    "+ ub2 + " != " + child.upperBound);
//					System.exit(-1);
//				}
//				if(ub2 != ub3) {
//					if(ub2 < ub3) {
//						//error:
//						System.out.println("Error: ub2 != ub3    "+ ub2 + " != " + ub3);
//						System.exit(-1);
//					} else { //ub3 is smaller because opt3 makes us skip more edges
//						betterUBcnt++;
//					}
//				}
//			}
//			if(BnBNode.OPT2 && BnBNode.OPT3) {
//				if(ub3 != child.upperBound) {
//					//error:
//					System.out.println("Error: ub3 != child.upperBound    "+ ub3 + " != " + child.upperBound);
//					System.exit(-1);
//				}
//				if(ub2 != ub3) {
//					if(ub2 < ub3) {
//						//error:
//						System.out.println("Error: ub2 != ub3    "+ ub2 + " != " + ub3);
//						System.exit(-1);
//					} else { //ub3 is smaller because opt3 makes us skip more edges
//						betterUBcnt++;
//					}
//				}
//			}
			
			
			sc.add(cv); // update sx with cv: to be used for the next child/sibling
			//System.out.println("DEBUG: candidates:"+candidatesList.size()+ " excluded: "+nc.excluded_vertices.size()+ " Sx: "+sx.size());

			if (child.getNumOfVertices() == CHkS.k){
				Top_Solutions.insertSolution(child); //finding Top-k solutions
			}
			if (child.getWeight() >= m_bestCost) { //test for best cost
				if (child.getWeight() > m_bestCost) {
					//System.out.println(child.getWeight());
					m_bestCost = child.getWeight();
					bestSolutions.clear();
				}
				bestSolutions.add(child);
				//System.out.println("##### new best solution (of "+bestSolutions.size()+") found #####");
				//System.out.println(child);
				//System.out.println();
			}
			
//			if(visitedNodes %100 == 0)
//			System.out.println("visitedNodes: "+visitedNodes );
//			System.out.println(p);
//			System.out.println();

			if ((child.upperBound >=  m_bestCost) && (child.vertices.size() < m_maxVertices)) { //test for expanding or else pruning
//				System.out.println(child);

				dfs_search(child); //########## DFS recursive call ##########

			} else { //child node is pruned
				//System.out.println(p);
				//this is the case: if ((upperBound <  m_bnb.m_bestCost) || (vertices.size() >= m_bnb.m_maxVertices))

				if ((child.upperBound <  m_bestCost) && (child.vertices.size() < m_maxVertices)) {
					prunedDueToUpperBound++;  //only due to UB: before we reach the max num of vertices
				}
				if (child.vertices.size() >= m_maxVertices) {
					prunedDueToVertices++; //due to max num of vertices (regardless of UB)
				}
			}
			
		} //end for each candidate vertex
	}
	
	/**
	 * BFS using heap implementation<br>
	 * BnBNodesHeap nodeHeap is the priority queue with the BFS search nodes.
	 * <br>
	 * 
	 * @param nodeHeap
	 */
	public void bfs_search(PriorityQueue<BnB_Node> nodeHeap) {
		
		while(!nodeHeap.isEmpty()) {
			BnB_Node p = nodeHeap.remove();
//			System.out.println("<bfs_search> visiting next top node:\n"+p);
			visitedNodes++;
			
	// 		Set<Vertex> candidateSet = new HashSet<Vertex>(); // candidate vertices for spawning children
			List<Vertex> candidatesList = new ArrayList<Vertex>(); // candidate vertices for spawning children
	
			if(verb) 
				System.out.println("<bfs_search> The neighbours of the last vertex: <" + p.last_v + ">" + "are:" + p.last_v.getNeighbourVertices());
	
			//1. add the candidates inherited from the previous parent node:
			for (Vertex c: p.candidate_vertices) {
				if ((!(p.vertices.contains(c)) && (!(p.excluded_vertices.contains(c)))))
					if(!candidatesList.contains(c))
						candidatesList.add(c);
			}
		
			// 2. add new candidates from the neighbours of the last added vertex
			for (Vertex vn : p.last_v.getNeighbourVertices()) {
				if ((!(p.vertices.contains(vn)) && (!(p.excluded_vertices.contains(vn))))) {
	//				candidateSet.add(vn);
					if(!candidatesList.contains(vn))
						candidatesList.add(vn);
				}
			}
			
			// Sorting the candidates set based on weighted degree
	//		candidatesList.addAll(candidateSet);
	//		Collections.sort(candidatesList, Collections.reverseOrder()); //??
			Collections.sort(candidatesList, new VertexComparator_byWeight_reverse());
	//		if(verb) 
	//		if(p.last_v.getID() == 10)
	//			System.out.println("<dfs_search> The sorted candidates are: " + candidatesList+"\n");
	
			p.setCandidateVertices(candidatesList); //set the updated candidates list into this parent node before generating child nodes
	//		if(verb) 
	//			System.out.println("<dfs_search> The sorted candidates at the parent are: " + p.candidate_vertices+"\n");
	
			Set<Vertex> sc = new HashSet<Vertex>(); //set of excluded candidate vertices for the sibling nodes
			for (Vertex cv: candidatesList) { // generate a child node for each candidate vertex
				BnB_Node child = new BnB_Node(p); //generate a child node: a copy of the parent node
	//			if(verb)
	//				System.out.println("<dfs_search> child with candidate <"+cv.getID()+">:\n"+child);
				
				child.excludeVertices(sc); //BUG found: this was done after expand, which also used to compute the UB
				child.expandSolution(cv); // adds a new vertex and its edges to the current solution
				child.computeUpperBound();
				
	//			if(verb) {
	//				System.out.println("<dfs_search> child with candidate <"+cv.getID()+">:\n"+child);
	//				System.out.println("<dfs_search> the parent candidates were: "+p.candidate_vertices);
	//				System.out.println("<dfs_search> the child candidates are:   "+child.candidate_vertices);
	//			}
				
	//			double ub2 = child.computeUB(true, false);
	//			double ub3 = child.computeUB(true, true);
	//			if(BnBNode.OPT2 && (!BnBNode.OPT3)) {
	//				if(ub2 != child.upperBound) {
	//					//error:
	//					System.out.println("Error: ub2 != child.upperBound    "+ ub2 + " != " + child.upperBound);
	//					System.exit(-1);
	//				}
	//				if(ub2 != ub3) {
	//					if(ub2 < ub3) {
	//						//error:
	//						System.out.println("Error: ub2 != ub3    "+ ub2 + " != " + ub3);
	//						System.exit(-1);
	//					} else { //ub3 is smaller because opt3 makes us skip more edges
	//						betterUBcnt++;
	//					}
	//				}
	//			}
	//			if(BnBNode.OPT2 && BnBNode.OPT3) {
	//				if(ub3 != child.upperBound) {
	//					//error:
	//					System.out.println("Error: ub3 != child.upperBound    "+ ub3 + " != " + child.upperBound);
	//					System.exit(-1);
	//				}
	//				if(ub2 != ub3) {
	//					if(ub2 < ub3) {
	//						//error:
	//						System.out.println("Error: ub2 != ub3    "+ ub2 + " != " + ub3);
	//						System.exit(-1);
	//					} else { //ub3 is smaller because opt3 makes us skip more edges
	//						betterUBcnt++;
	//					}
	//				}
	//			}
				
				
				sc.add(cv); // update sx with cv: to be used for the next child/sibling
				//System.out.println("DEBUG: candidates:"+candidatesList.size()+ " excluded: "+nc.excluded_vertices.size()+ " Sx: "+sx.size());

				if (child.getNumOfVertices() == CHkS.k){
					Top_Solutions.insertSolution(child); //finding Top-k solutions
				}
				if (child.getWeight() >= m_bestCost) { //test for best cost
					if (child.getWeight() > m_bestCost) {
						//System.out.println(child.getWeight());
						m_bestCost = child.getWeight();
						bestSolutions.clear();
					}
					bestSolutions.add(child);
					//System.out.println("##### new best solution (of "+bestSolutions.size()+") found #####");
					//System.out.println(child);
					//System.out.println();
				}

	//			if(nc.id %50000 == 0)
	//				System.out.println(nc);
	
				if ((child.upperBound >=  m_bestCost) && (child.vertices.size() < m_maxVertices)) { //test for expanding or else pruning
	//				System.out.println(nc);
	
	//				dfs_search(child); //########## DFS recursive call ##########
					nodeHeap.add(child); //########## add the new node to the BFS heap ##########
					if(nodeHeap.size() > maxHeap_size) maxHeap_size = nodeHeap.size();

//					System.out.println("\n\n--- HEAP --- \n");
//					for(BnBNode x: nodeHeap) {
//						System.out.println(x);
//					}
//					System.out.println("--- --- --- \n\n");
					
				} else { //child node is pruned
					//System.out.println(nc);
					//this is the case: if ((upperBound <  m_bnb.m_bestCost) || (vertices.size() >= m_bnb.m_maxVertices))
	
					if ((child.upperBound <  m_bestCost) && (child.vertices.size() < m_maxVertices)) {
						prunedDueToUpperBound++;  //only due to UB: before we reach the max num of vertices
					}
					if (child.vertices.size() >= m_maxVertices) {
						prunedDueToVertices++; //due to max num of vertices (regardless of UB)
					}
				}
				
			} //end for each candidate vertex
		} //end while heap is not empty
	}

	/**
	 * Hybrid BFS+DFS<br>
	 * <br>
	 * 
	 * @param nodeHeap
	 */
	public void hfs_search(PriorityQueue<BnB_Node> nodeHeap) {
		while(!nodeHeap.isEmpty()) {
			BnB_Node p = nodeHeap.remove();
			
			if(nodeHeap.size() > MAXHEAPSIZE) {
//				System.out.println("<hfs_search> visiting next top node in DFS mode:\n"+p);
				dfs_search(p);   //########## DFS recursive call ##########
			} else { //BFS

//			System.out.println("<hfs_search> visiting next top node in BFS mode:\n"+p);
			visitedNodes++;



	// 		Set<Vertex> candidateSet = new HashSet<Vertex>(); // candidate vertices for spawning children
			List<Vertex> candidatesList = new ArrayList<Vertex>(); // candidate vertices for spawning children

			if(verb)
				System.out.println("<bfs_search> The neighbours of the last vertex: <" + p.last_v + ">" + "are:" + p.last_v.getNeighbourVertices());

			//1. add the candidates inherited from the previous parent node:
			for (Vertex c: p.candidate_vertices) {
				if ((!(p.vertices.contains(c)) && (!(p.excluded_vertices.contains(c)))))
					if(!candidatesList.contains(c))
						candidatesList.add(c);
			}

			// 2. add new candidates from the neighbours of the last added vertex
			for (Vertex vn : p.last_v.getNeighbourVertices()) {
				if ((!(p.vertices.contains(vn)) && (!(p.excluded_vertices.contains(vn))))) {
	//				candidateSet.add(vn);
					if(!candidatesList.contains(vn))
						candidatesList.add(vn);
				}
			}

			// Sorting the candidates set based on weighted degree
	//		candidatesList.addAll(candidateSet);
	//		Collections.sort(candidatesList, Collections.reverseOrder()); //??
			Collections.sort(candidatesList, new VertexComparator_byWeight_reverse());
	//		if(verb)
	//		if(p.last_v.getID() == 10)
	//			System.out.println("<dfs_search> The sorted candidates are: " + candidatesList+"\n");

			p.setCandidateVertices(candidatesList); //set the updated candidates list into this parent node before generating child nodes
	//		if(verb)
	//			System.out.println("<dfs_search> The sorted candidates at the parent are: " + p.candidate_vertices+"\n");

			Set<Vertex> sc = new HashSet<Vertex>(); //set of excluded candidate vertices for the sibling nodes
			for (Vertex cv: candidatesList) { // generate a child node for each candidate vertex
				BnB_Node child = new BnB_Node(p); //generate a child node: a copy of the parent node
	//			if(verb)
	//				System.out.println("<dfs_search> child with candidate <"+cv.getID()+">:\n"+child);

				child.excludeVertices(sc); //BUG found: this was done after expand, which also used to compute the UB
				child.expandSolution(cv); // adds a new vertex and its edges to the current solution
				child.computeUpperBound();

	//			if(verb) {
	//				System.out.println("<dfs_search> child with candidate <"+cv.getID()+">:\n"+child);
	//				System.out.println("<dfs_search> the parent candidates were: "+p.candidate_vertices);
	//				System.out.println("<dfs_search> the child candidates are:   "+child.candidate_vertices);
	//			}

	//			double ub2 = child.computeUB(true, false);
	//			double ub3 = child.computeUB(true, true);
	//			if(BnBNode.OPT2 && (!BnBNode.OPT3)) {
	//				if(ub2 != child.upperBound) {
	//					//error:
	//					System.out.println("Error: ub2 != child.upperBound    "+ ub2 + " != " + child.upperBound);
	//					System.exit(-1);
	//				}
	//				if(ub2 != ub3) {
	//					if(ub2 < ub3) {
	//						//error:
	//						System.out.println("Error: ub2 != ub3    "+ ub2 + " != " + ub3);
	//						System.exit(-1);
	//					} else { //ub3 is smaller because opt3 makes us skip more edges
	//						betterUBcnt++;
	//					}
	//				}
	//			}
	//			if(BnBNode.OPT2 && BnBNode.OPT3) {
	//				if(ub3 != child.upperBound) {
	//					//error:
	//					System.out.println("Error: ub3 != child.upperBound    "+ ub3 + " != " + child.upperBound);
	//					System.exit(-1);
	//				}
	//				if(ub2 != ub3) {
	//					if(ub2 < ub3) {
	//						//error:
	//						System.out.println("Error: ub2 != ub3    "+ ub2 + " != " + ub3);
	//						System.exit(-1);
	//					} else { //ub3 is smaller because opt3 makes us skip more edges
	//						betterUBcnt++;
	//					}
	//				}
	//			}


				sc.add(cv); // update sx with cv: to be used for the next child/sibling
				//System.out.println("DEBUG: candidates:"+candidatesList.size()+ " excluded: "+nc.excluded_vertices.size()+ " Sx: "+sx.size());
				if (child.getNumOfVertices() == CHkS.k){
					Top_Solutions.insertSolution(child); //finding Top-k solutions
				}
				if ((child.getWeight() >= m_bestCost) ) { //test for best cost
					if (child.getWeight() > m_bestCost) {
						m_bestCost = child.getWeight();
						bestSolutions.clear();
						//System.out.println(child.getWeight());
					}
					bestSolutions.add(child);
					//System.out.println("##### new best solution (of "+bestSolutions.size()+") found #####");
					//System.out.println(child.getWeight());
					//System.out.println();
				}

	//			if(nc.id %50000 == 0)
	//				System.out.println(nc);
	
				if ((child.upperBound >=  m_bestCost) && (child.vertices.size() < m_maxVertices)) { //test for expanding or else pruning
	//				System.out.println(nc);
	
	//				
//					if(nodeHeap.size() < MAXHEAPSIZE) {
						nodeHeap.add(child); //########## add the new node to the BFS heap ##########
						if(nodeHeap.size() > maxHeap_size) maxHeap_size = nodeHeap.size();
//					} else { //DFS
//						dfs_search(child);   //########## DFS recursive call ##########
//					}
					
//					System.out.println("\n\n--- HEAP --- \n");
//					for(BnBNode x: nodeHeap) {
//						System.out.println(x);
//					}
//					System.out.println("--- --- --- \n\n");
					
				} else { //child node is pruned
					//System.out.println(nc);
					//this is the case: if ((upperBound <  m_bnb.m_bestCost) || (vertices.size() >= m_bnb.m_maxVertices))
	
					if ((child.upperBound <  m_bestCost) && (child.vertices.size() < m_maxVertices)) {
						prunedDueToUpperBound++;  //only due to UB: before we reach the max num of vertices
					}
					if (child.vertices.size() >= m_maxVertices) {
						prunedDueToVertices++; //due to max num of vertices (regardless of UB)
					}
				}
				
			} //end for each candidate vertex
			} //end of IF-BFS
		} //end while heap is not empty
	}
	
}// end class
