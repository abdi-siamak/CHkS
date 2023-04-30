package uk.ac.rdg.chks.bb;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

//import com.sun.jdi.IntegerValue;
import uk.ac.rdg.chks.graph.Edge;
import uk.ac.rdg.chks.graph.Vertex;

/**
 * The Branch and Bound node <br>
 * Each BnBNode instance represents a candidate subgraph solution G=(V',E'). 
 * The weight (cost) is the sum of the weights of edge in E'.
 * The upperBound is an estimation of the maximum cost for the nodes in the branch rooted at this node.
 * <br> 
 * 
 * @author Giuseppe Di Fatta and Salwa Musa
 */
public class BnB_Node { //implements Comparable<BnBNode> {
	private static boolean verb = false; //to switch ON/OFF a number of messages in stdout
	private static boolean warn = false; //to switch ON/OFF warning messages in stdout

	/** 
	 * The BnB search run this BnBNode belongs to.
	 */
	private BnB_Search m_bnb;
	
	/**
	 * This BnB node unique ID
	 */
	private int id; 
	private static int nextID = 0;

	/**
	 * The parent node in the search tree.
	 */
	private BnB_Node parent;

	/** 
	 * The weight (cost) of the solution. This is also a trivial lowerbound (LB).
	 */
	private double weight;

	/**
	 * the upperbound (UB) for the solutions in the branch rooted at this node
	 */
	double upperBound;

	/**
	 * the vertices of the subgraph 
	 */
	ArrayList<Vertex> vertices;
	
	/** 
	 * the edges of the subgraph: <br>
	 * a set of integers, indicating the ordinal position 
	 * of the edge in the ordered edge list of the entire term graph.
	 */
	private HashSet<Integer> edges;

	/** 
	 * This set contains the candidates vertices to be used to extend the current subgraph and generate child nodes.
	 */
//	Set<Vertex> candidate_verticess;// =new HashSet<Vertex>();
	ArrayList<Vertex> candidate_vertices;

	/** 
	 * This set contains the excluded vertices NOT to be added in the current subgraph
	 */
	Set<Vertex> excluded_vertices;
	
	/** 
	 * The last vertex added to the current subgraph
	 */
	Vertex last_v;

	/**
	 * the level of the search tree of this node
	 */
	private int treeLevel;
	
	
	//counting the number of optimisation savings in computeUpperBound
	static int opt1Cnt=0; 	//OPT1: exclude the edges already in the subgraph 
	static int opt2Cnt=0; 	//OPT2: exclude edge if one of its vertices is excluded
	static int opt3Cnt=0;	//OPT3: exclude edge if one of its vertices is too far from the current solution
//	static int opt3bCnt=0;	//OPT3: for debugging/info

	public static final boolean OPT2 = true;
	public static final boolean OPT3 = false;

	
	/**
	 * constructor for a root node of a Branch and Bound search tree: 
	 * this BnB node contains one vertex and no edges.
	 */
	public BnB_Node(BnB_Search bnb, Vertex v) { //int k, TermGraph tg,
		m_bnb = bnb;
		
		last_v = v;
		parent = null;
		treeLevel = 0;
		vertices = new ArrayList<Vertex>();
		vertices.add(v);

		edges = new HashSet<Integer>();

		excluded_vertices = new HashSet<Vertex>();
//		excluded_vertices.add(v);
		
		candidate_vertices = new ArrayList<Vertex>(); //new HashSet<Vertex>();

		weight = 0;
		upperBound = 0;

		// the upper bound is the sum of the weights of m_maxEdges=(k*k-1)/2 edges
		for (int i = 0; i < m_bnb.m_maxEdges; i++) {
			upperBound += m_bnb.theTG.getEdge(i).getWeight();
		}

//		if(verb)
//			System.out.println("Root node with vertex <"+v.getID()+">: UB = " + upperBound + "   (cost for top "+m_bnb.m_maxEdges+" edges)");
	}
	

	/**
	 * Constructor for a generic BnB node (not a root node).
	 * It creates a copy of the parent node, still to be extended with 
	 * a new vertex by using the method {@link #expandSolution(Vertex)}.
	 * 
	 * @param parentNode
	 */
	public BnB_Node(BnB_Node parentNode) {
		m_bnb = parentNode.m_bnb;
		last_v = null; //v;  //we do that in expandSolution
		parent = parentNode;
		id=++nextID;
		treeLevel = parentNode.treeLevel + 1;
 		this.weight = parentNode.weight;
		this.upperBound = 0;

		// this array has the all vertices in the subgraph
		vertices = new ArrayList<Vertex>();
		for (Vertex u : parentNode.vertices) {
			vertices.add(u);
		}
		//vertices.add(last_v); //we do that in expandSolution

		// this set has the all edges in the subgraph
		edges = new HashSet<Integer>();
 
		// this.edges_subgraph = parentNode.edges_subgraph;
		for (Integer e : parentNode.edges) {
			edges.add(e);
		}
	
		// this HashSet for the candidates
		candidate_vertices = new ArrayList<Vertex>(); //new HashSet<Vertex>();
		for (Vertex c : parentNode.candidate_vertices) {
//			if ((c != last_v))  //we remove it in expandSolution
			
//			if(!candidate_vertices.contains(c)) //TODO: redundant ==> to be removed but first test with it.
 				candidate_vertices.add(c);
		}
 
		// this HashSet for the excluded Vertices
		excluded_vertices = new HashSet<Vertex>();
		for (Vertex ve : parentNode.excluded_vertices) {
			excluded_vertices.add(ve);
		}
	}

	public int containsSomeVerticesFrom(int vIDs[]) {
		int cnt=0;
		for(int x: vIDs) {
			for(Vertex v: vertices) {
				if(x == v.getID()) {
					cnt++;
					break;
				}
			}
		}
		return cnt;
	}
	
	public int containsOnlyVerticesFrom(int vIDs[]) {
		int cnt=0;
		for(Vertex v: vertices) {
			boolean found = false;
			for(int x: vIDs) {
				if(x == v.getID()) {
					cnt++;
					found = true;
					break;
				}
			}
			if(!found) {
				cnt = 0;
				break;
			}
		}
		return cnt;
	}
	
	 public static int getTotNumOfNode() {
		 return nextID;
	 }
	 
	/** 
	 * Adding a new vertex (last_v) to the solution along with all its induced edges 
	 * (i.e., edges from this vertex to any vertex already in the solution.
	 */
	public void expandSolution(Vertex newVertex) {
		last_v = newVertex;
		if(candidate_vertices.contains(last_v))
			candidate_vertices.remove(last_v);
		
		for (Vertex in: vertices) {
			int index = last_v.indexOfNeighbour(in);
			if (index > -1) {
				// numberIncludedVertices++;
				weight += last_v.getWeight(index);
				edges.add(last_v.getOrder(index));

//				if(verb) 
//					System.out.println("\t\tadded edge weight: " + last_v.getWeight(index) + " for edge <" + last_v + "," + in + ">");
			}
		}
		
		vertices.add(last_v);
//		computeUpperBound(); //now brought forward to make it explicit/visible

		if(verb) 
			System.out.println("vertex "+ last_v + " and its induced edges have been added. The UB has been updated.");
	}
	
	/**
	 * Exclude a set of vertices from possible expansions in nodes of this BnB branch.
	 * @param excludedV
	 */
	public void excludeVertices(Set<Vertex> excludedV) {
		excluded_vertices.addAll(excludedV);
	}

//	public void addCandidateVertices(Set<Vertex> vset) { //excludedV) {
////		if (!(excludedV.contains(last_v))) { //BUG: what was that?
//		for(Vertex v: vset) {
//			if(!excluded_vertices.contains(v)) {
//				candidate_vertices.add(v); //excludedV);
//			}
//		}
//	}
	
//	public void setCandidateVertices(Set<Vertex> vset) {
//		vertices_Candidates.addAll(vset);
//	}

	public void setCandidateVertices(List<Vertex> vList) {
		candidate_vertices.clear();
		
//		candidate_vertices.addAll(vList);
		for (Vertex v : vList) {
			if (!excluded_vertices.contains(v)) {
				candidate_vertices.add(v);
			}
		}
	}
	


	//**********************************************************
		/**
		 * Calculating the upper Bound
		 * 
		 * next upper bound = the weight of subgraph + (the sum of the next
		 * edges(2 k)-(2 s) so we need the maxEdges numbers and add their
		 * weights to calculate the upper bound.
		 * Therefore, we have to start the loop with 
		 * that will give us the next index to the edge which will start to get the weigh
		 * as we already dealt with other edges. #then  will stop the loop when we reach to index maxEdges 
		 * where MaxEdges will give the max number of next edges to calculate UB and we will put in our account the the other edhges 
		 * in the list which have been used before.
		 * if the end points of the edges are in excluding vertices we will not consider this edge, 
		 * by doing thing will tight the UB 
		 */
	public void computeUpperBound(){
		//It is important to implement, test and compare the optimisation methods we discussed about the upperbound.
		//0. all edges can be considered: we just take the top maxNextEdges weights in the sorted edge list
		//1. edges already included cannot be used: OPT1
		//2. edges that refer to excluded vertices cannot be used: OPT2
		//3. edges that have any of its two vertices far away from the current subgraph (further than k-|V| hops) cannot be included: OPT3.
		//TODO: Anything else?
		
		upperBound = weight;
		int numExtraEdges = choose2(m_bnb.m_maxVertices) - choose2(vertices.size());
		
		//select the top edges that could be added in this branch:
		Iterator<Edge> it = m_bnb.theTG.getEdgeIterator();
		Edge e;
		int cnt = numExtraEdges;
//		int x=10;
		while(it.hasNext()) {
			if(cnt == 0) break;
			e = it.next();
			
			if(edges.contains(e.getOrder())) {      //OPT1: exclude the edges already in the subgraph 
				opt1Cnt++;
				continue;
			}

			if(OPT2) { //OPT2: exclude edge if one of its vertices is excluded
				if(excluded_vertices.contains(e.getVertex1()) || excluded_vertices.contains(e.getVertex2())) { //OPT2: exclude edge if one of its vertices is excluded
					opt2Cnt++;
					continue;
				}
			}
			
			if(OPT3) { //OPT3: exclude edge if one of its vertices is too far from the current solution
				int maxd = m_bnb.m_maxVertices - vertices.size();
				boolean withinRadius1 = false;
				boolean withinRadius2 = false;
				for(Vertex v: vertices) {
					if((!withinRadius1) && (v.isWithinDistance(e.getVertex1(), maxd))) withinRadius1 = true;
					if((!withinRadius2) && (v.isWithinDistance(e.getVertex2(), maxd))) withinRadius2 = true;
					if(withinRadius1 && withinRadius2) break; //no reason to continue
				}
				if(!(withinRadius1 && withinRadius2)) { // if((!withinRadius1) || (!withinRadius2)) {
					opt3Cnt++;
					continue;
				}
			}
			
			//using this edge for estimating the UB:
			upperBound += e.getWeight();
		    cnt--;
		}
		if(cnt > 0) {
			if(warn)
				System.out.println("Warning: run out of edges but still need " + cnt + 
						" (numNextEdges="+numExtraEdges+")");
		}
   }

	public double computeUB(boolean opt2, boolean opt3){
		double ub = weight;
		int numExtraEdges = choose2(m_bnb.m_maxVertices) - choose2(vertices.size());
		
		Iterator<Edge> it = m_bnb.theTG.getEdgeIterator();
		Edge e;
		int cnt = numExtraEdges;
		while(it.hasNext()) {
			if(cnt == 0) break;
			e = it.next();
			
			if(edges.contains(e.getOrder())) {      //OPT1: exclude the edges already in the subgraph 
//				opt1Cnt++;
				continue;
			}

			if(opt2) { //OPT2: exclude edge if one of its vertices is excluded
				if(excluded_vertices.contains(e.getVertex1()) || excluded_vertices.contains(e.getVertex2())) { //OPT2: exclude edge if one of its vertices is excluded
//					opt2Cnt++;
					continue;
				}
			}
			
			if(opt3) { //OPT3: exclude edge if one of its vertices is too far from the current solution
				int maxd = m_bnb.m_maxVertices - vertices.size();
				boolean withinRadius1 = false;
				boolean withinRadius2 = false;
				for(Vertex v: vertices) {
					if(v.isWithinDistance(e.getVertex1(), maxd)) withinRadius1 = true;
					if(v.isWithinDistance(e.getVertex2(), maxd)) withinRadius2 = true;
				}
				if((!withinRadius1) || (!withinRadius2)) {
//					opt3Cnt++;
					continue;
				}
			}
			
			
			//using this edge for estimating the UB:
			ub += e.getWeight();
		    cnt--;
		}
		if(cnt > 0) {
			if(warn)
				System.out.println("Warning: run out of edges but still need " + cnt + 
						" (numNextEdges="+numExtraEdges+")");
		}
		
		return ub;
   }

	public boolean isPruned() {
		if ((upperBound <  m_bnb.m_bestCost) || (vertices.size() >= m_bnb.m_maxVertices)) 
			return true;
		else 
			return false;
	}

	public double getWeight() {
		return weight;
	}

	public double getUpperBound() {
		return upperBound;
	}
	
	private int choose2(int i) {
		return i * (i - 1) / 2;
	}

	/**
	 * returns the number of vertices in the subgraph (current solution)
	 * @return
	 */
	public int getNumOfVertices() {
		return vertices.size();
	}

	public ArrayList<Vertex> getVertices() {
		return vertices;
	}
	public Set<Vertex> getVertices_excluded() {
		return excluded_vertices;
	}
	
//	public Set<Vertex> get_verticesCandidates() {
//		return candidate_vertices;
//	}

	public int getTreeLevel() {
		return treeLevel;
	}

//	//TODO: need to change this into a comparator in a separate class.
//	@Override
//	public int compareTo(BnBNode n) {
//	   
//       if (this.upperBound > n.upperBound)
//			return -1;
//		else if (this.upperBound == n.upperBound)
//			return 0;
//		else
//			return 1;
//	}

	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer();
		/*
		buf.append(  "----- BnB node#"+id+"  |V|="+vertices.size() + " |E|="+ edges.size() + " LB="+weight+" UB="+upperBound +" -----"
		           + "  (current bestCost="+m_bnb.m_bestCost+")\n");
		buf.append((parent!=null) ? "(parent=" + parent.id+ ") " : ("(root node)"));
		buf.append(((isPruned()) ? " (pruned)" : " (not pruned)") + " (tree level: " + treeLevel+")");
		buf.append("\n");
		
		buf.append("Vertices ("+vertices.size()+"):");
 		for (Vertex v : vertices) {
			buf.append(" <" + v+">");
			//buf.append(" <" + CHkS.dictionary.get(v) + ">");
		}
		buf.append("\n");
		 */
		////////////////////////////////////////////////////////////////////
		//buf.append("             ");
		for (Vertex v : vertices) {
			buf.append(CHkS.dictionary.get(v.getID()) + " ");
		}
		//buf.append("\n");
		/////////////////////////////////////////////////////////////////////
/*
		buf.append("Edges ("+edges.size()+"):");

		for (Integer e : edges) {
			buf.append(" " +m_bnb.theTG.getEdge(e).toString());
		}
		buf.append("\n");
		/////////////////////////////////////////////////////////////////////
		buf.append("          ");
		for (Integer e : edges) {
			buf.append(" (<" +CHkS.dictionary.get(m_bnb.theTG.getEdge(e).getVertex1().getID())+",");
			buf.append(" " +CHkS.dictionary.get(m_bnb.theTG.getEdge(e).getVertex2().getID())+">,");
			buf.append(" " +m_bnb.theTG.getEdge(e).getWeight()+")");
		}
		buf.append("\n");
		/////////////////////////////////////////////////////////////////////
		final int MAX=10;
*/
		/*
		if(verb) {
			buf.append("Candidate vertices ("+candidate_vertices.size()+"):");
			int max = MAX; //candidate_vertices.size()<MAX?candidate_vertices.size():MAX;
			for (Vertex v: candidate_vertices) {
				buf.append(" " + v.getID()+"");
				if(--max == 0) {
					buf.append(" ...");
					break;
				}
			}
			buf.append("\n");	
			
			buf.append("Excluded vertices ("+excluded_vertices.size()+"):");
			max = MAX; //excluded_vertices.size()<MAX?excluded_vertices.size():MAX;
			for(Vertex v: excluded_vertices) {
				buf.append(" " + v.getID()+"");
				if(--max == 0) {
					buf.append(" ...");
					break;
				}
			}
			buf.append("\n");
		}
		*/
		//buf.append("-----");
		return buf.toString();
	}
}
