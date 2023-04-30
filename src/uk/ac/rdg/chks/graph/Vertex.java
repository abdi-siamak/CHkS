package uk.ac.rdg.chks.graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/** 
 * The vertex class
 *
 */
public class Vertex implements Comparable<Vertex> {

	// The id of the vertex
	private int v_ID;
	// The number of neighbours of the vertex.
	private int degree;
	// The sum of weights of neighbours of the vertex.
	private double weighted_degree;
	// private int weighted_degree;
	// private ArrayList<Double> weightedDegree;

	private double top_weighted_degree; //weighted_degree for only top edges


	public ArrayList<Vertex> adjacentVertices;// adjacentvertex will be object of array list array as each element is
												// represent vertex and each element of this array is a list [0, 7, 2]
												// this list represent the nighbours of each vertex
	private ArrayList<Double> edgeWeights;

	private ArrayList<Edge> edges;// the edges to the neighbours
	
	private Map<Integer,Set<Vertex>> kNN;
	
	// constructor
	public Vertex(int ID) { // id of each vertex and in the file manger we have vertices[u]
		this.v_ID = ID;
		degree = 0;
		weighted_degree = 0.0;
		adjacentVertices = new ArrayList<Vertex>();
		edgeWeights = new ArrayList<Double>();
		edges = new ArrayList<Edge>();
		// weightedDegree =new ArrayList<Double>();
		// vertexEdges=new HashMap<>();
		
		kNN = new HashMap<Integer,Set<Vertex>>();
	}

	/**
	 * Compute the neareast neigbours up to a max number of hops (distance).
	 */
	public void compute_kNN(int maxDist) {
		Set<Vertex> nn = null;
		
		//d=0, the vertex itself
		nn= new TreeSet<Vertex>();
		nn.add(this);
		kNN.put(0,nn);
		
//		//d=1, the adjacent vertices
//		nn = new TreeSet<Vertex>();
//		for(Vertex v: adjacentVertices) {
//			nn.add(v);
//		}
//		kNN.put(1,nn);

		//d>0, the vertices that are distant d=1..maxdist from this vertex
		Set<Vertex> prev_nn = nn;
		for(int d=1; d<=maxDist; d++) {
			nn = new TreeSet<Vertex>();
			for(Vertex n: prev_nn) {
				for(Vertex v: n.adjacentVertices) {
					if(v!=this) {
						boolean found=false;
						for(int kk=(d-1); kk>0; kk--) {
							if(kNN.get(kk).contains(v)) {
								found=true;
								break;
							}
						}
						if(!found) {
							nn.add(v);
//							if(!nn.contains(v))
						}
					}
				}
			}
			kNN.put(d,nn);
			prev_nn = nn;
		}
	}
	
	/**
	 * return the nearest neighbours at a given distance
	 * 
	 * @param dist
	 * @return
	 */
	public Set<Vertex> get_NN(int dist){
		return kNN.get(dist);
	}
	
	public static void compute_kNN_for_all_vertices(TermGraph tg, int maxDist) {
		for (Vertex v: tg.getVertices()) {
			v.compute_kNN(maxDist);
		}
	}
	
	public void print_kNN() {
		System.out.println("--- kNN for vertex <" + this.v_ID+">");
		for(Integer k: kNN.keySet()) {
			Set<Vertex> nn = kNN.get(k);
			System.out.print("dist<"+k+"> size="+nn.size()+": ");
			for(Vertex v: nn) {
				System.out.print("<"+v.v_ID + "> ");
			}
			System.out.println();
		}
		System.out.println("---");
	}
	
//	public void compute_top_weighted_degree(int maxedges) {
////		System.out.println("<compute_top_weighted_degree>...");
//		Collections.sort(edges, new EdgeComparator_byWeight());
//		Edge e;
//		int min = (edges.size() < maxedges) ? edges.size() : maxedges;
//		for(int i=0; i<min; i++) {
//			e = edges.get(i);
////			System.out.println("<compute_top_weighted_degree> " + e);
//			top_weighted_degree += e.getWeight();
//		}
////		if(top_weighted_degree != weighted_degree) {
////			top_WD_opt_cnt++;
////			System.out.println("<compute_top_weighted_degree> top_weighted_degree opt cnt = "+
////					top_WD_opt_cnt + "      (" + top_weighted_degree + "     != " + weighted_degree+")");
////		}
//	}
	
	
	public int getID() {
		return v_ID;
	}

	public void setID(int ID) {
		this.v_ID = ID;
	}

	public void addEdge(Edge e) {
		edges.add(e);
	}

	public void addNeighbour(Vertex v, double weight) {

		adjacentVertices.add(v);
		edgeWeights.add(weight);
		// setWeight(v.v_ID,weight);
		degree++;// degree represent how many neighbours this vertex have
		weighted_degree = weighted_degree + weight;

	}

	
	public void resetWightD() {
		weighted_degree = 0.0;
	}

	public boolean isNeighbour_(Vertex v) {
		boolean flag = false;
		for (Vertex u : adjacentVertices) {
			if (u.v_ID == v.v_ID) {
				flag = true;
				break;
			}
		}
		return flag;
	}
	
	public boolean isWithinDistance(Vertex v, int maxd) {
		boolean found = false;
		for(int k=0; k<=maxd; k++) {
//		for(Integer k: kNN.keySet()) {
//			if(k > maxd) break;
			Set<Vertex> nn = kNN.get(k);
			if((nn != null) && (nn.contains(v))) {
				found = true;
				break;
			}
		}
		return found;
	}


	public void printAdjNeighbours() {
		System.out.println("At vertex <" + this.v_ID + "> adjacentVertices are ("+adjacentVertices.size()+"): ");
		for (Vertex v : adjacentVertices) {
			System.out.print(" " + v.v_ID);
		}
		System.out.println();
	}

	public void printNeighbours(int maxd) {
		System.out.println("Vertex " + this.v_ID);
		for(int k=0; k<=maxd; k++) {
			Set<Vertex> nn = kNN.get(k);
			if(nn != null) {
				System.out.print("Neighbours at distance "+k+ "("+nn.size()+"):");
				int max=10;
				for (Vertex v : nn) {
					System.out.print(" " + v.v_ID);
					if(--max == 0) {
						System.out.print(" ...");
						break;
					}
				}
				System.out.println();
			} else {
				System.out.println("Neighbours at distance "+k+ ": NONE");
			}
		}
	}
	
	public int getDegree() {
		return degree;
	}

	public double get_WeightedDegree() {
		// weightedDegree.add(weighted_degree);// this duplicated we have made this in
		// set degree but once we wont call the setweighted degree we added here but
		// that should be removed
		return weighted_degree;
	}

	public double get_topWeightedDegree() {
		return top_weighted_degree;
	}
	
	public ArrayList<Vertex> getNeighbourVertices() {
		return adjacentVertices;
	}

	public ArrayList<Double> getEdgeWeights() {
		return edgeWeights;
	}

	

	public ArrayList<Edge> getEdgeList() {
		return edges;
	}

	public double getWeight(int i) {
		return edgeWeights.get(i);
	}

	public int getOrder(int index) {
		return edges.get(index).getOrder();
	}


	public void printEdges(int M) {// M is TotalEdges
		System.out.println("Vertex : " + v_ID);
		for (int i = 0; i < M && i < adjacentVertices.size(); i++) {
			System.out.println(v_ID + " " + adjacentVertices.get(i) + " : " + edgeWeights.get(i));
		}
	}

	public void print_adjacentVertices() {
		System.out.println("Vertex " + v_ID);
		if (adjacentVertices.size() > 0) {
			for (int i = 0; i < adjacentVertices.size(); i++) {

				System.out.println("vertex : " + v_ID + " is connected to " + adjacentVertices.get(i) + " with weight "
						+ edgeWeights.get(i));
			}

		} else {
			System.out.println("No neighbours");
		}
	}

	public void printEdgesweights() {

		System.out.println("Vertex : " + v_ID);
		for (int i = 0; i < edgeWeights.size(); i++) {
			System.out.println("the weight are==> " + edgeWeights.get(i));
		}
	}

	/* This Method for removing Vertices and Edges */
	public Edge removeVertex(Vertex target) {

		int i = adjacentVertices.indexOf(target);
		degree--;
		weighted_degree -= edgeWeights.get(i);
		edgeWeights.remove(i);
		adjacentVertices.remove(i);

		// Removing Edges
		Edge rem = null;
		for (Edge e : edges) {
			if (e.isIncident(target) && e.isIncident(this)) {
				// System.out.println("This is " + this.v_ID);
				rem = e;
				// break;
			}
		}
		// target.edges.remove(rem);
		edges.remove(rem);
		System.out.println("Removing Edge " + rem);
		return rem;
	}

	/*
	 * @Override
	 * 
	 * public int compareTo(Vertex v) { if(v !=null){ double
	 * compare=((Vertex)v).weighted_degree; //For Ascending order // return
	 * this.weighted_degree-compare;
	 * 
	 * For Descending order return (int) (compare-this.weighted_degree);
	 * 
	 * }
	 */

	// this method for sorting vertices based on weighted degree in ASC
//	public int compareTo(Vertex v) {
//
//		int value = 0;
//
//		/*
//		 * if (v != null) {
//		 * 
//		 * if (this.get_WeightedDegree() < v.get_WeightedDegree()) value = 1; else if
//		 * (this.get_WeightedDegree() == v.get_WeightedDegree()) value = 0;
//		 * 
//		 * else value =-1; }
//		 */
//
//		// return value;
//
//		// for ordering dec
//
//		if (v != null) {
//
//			if (this.get_WeightedDegree() > v.get_WeightedDegree())
//				value = 1;
//			else if (this.get_WeightedDegree() == v.get_WeightedDegree())
//				value = 0;
//
//			else
//				value = -1;
//		}
//
//		return value;
//	}

	public int compareTo(Vertex v) {
		int value = 0;
		if (v != null) {
			if (this.v_ID > v.v_ID)
				value = 1;
			else if (this.v_ID == v.v_ID)
				value = 0;
			else
				value = -1;
		}
		return value;
	}
	
	public int indexOfNeighbour(Vertex v) {
		return adjacentVertices.indexOf(v);
	}

	public void setWeight(int index, double w) {
		edgeWeights.set(index, w);
	}
	
	public String toString_wAdj() { //with adj
		StringBuffer b = new StringBuffer();
		b.append("<"+v_ID+"> adj: ");
		for (Vertex v : adjacentVertices) {
			b.append("<" + v.v_ID + "> ");
		}
		return b.toString();
	}

	public String toStringObj() {
		return ""+((Object)this).toString();
	}
	
	public String toString() {
		return ""+v_ID;
	}
}
