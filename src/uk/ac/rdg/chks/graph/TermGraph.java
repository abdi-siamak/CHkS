package uk.ac.rdg.chks.graph;

import static java.util.stream.Collectors.toMap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** 
 * The term graph class
 */
public class TermGraph {
	// Contains the vertices of the graph
//	private Vertex[] vertices;
	private ArrayList<Vertex> vertices;

	// The edges of the graph
	private ArrayList<Edge> edges;
	
	//***********************************************************
 	// Using for sorting the vertices 
	//TODO: can we eliminate this additional lists? Can we just use the two above?
	public List<Vertex> vertices_List = new ArrayList<Vertex>();
//	public List<Edge> edges_List = new ArrayList<Edge>();

	Map<Integer, Double> wDM = new HashMap<Integer, Double>();// store all terms with unique ID
   	ArrayList<Double> wDVl = new ArrayList<Double>();
   	
	//***********************************************************
   	
   	/**
   	 * Constructor
   	 */
   	public TermGraph (ArrayList<Vertex> theVertices, ArrayList<Edge> theEdges) {
   		vertices = theVertices;
   		edges = theEdges;
   	}
   	   	
	public int getNumVertices() {
		return vertices.size();
//		int cnt=0;
//		for (int i = 0; i < vertices.length; i++) {
//			if(vertices[i]!=null){
//				cnt++;
//			}
//		}
//		return cnt;
	}
	
	
//	public void computeTopWD(int max) {
//		for(Vertex v: vertices) {
//			v.compute_top_weighted_degree(max);
//		}
//	}

	public void setVertices(ArrayList<Vertex> v) {
		vertices = v;
	}

	public ArrayList<Vertex> getVertices() {
		return vertices;
	}
	
	// Returns a vertex with a specific index
	public Vertex getVertex(int index) {
		return vertices.get(index);
	}

	public Vertex getVertex_byID(int id) {
		for(Vertex v: vertices)
			if(v.getID() == id)
				return v;
		return null;
	}
	
	// Returns the total number of vertices
	public int getNumberOfVertices() {
		return getNumVertices();// calling the function which return number of vertices
	}

	
	public List<Vertex> getSortingVertices() {
		return vertices_List;
	}
	public void printAllVertices() {

		System.out.println("The vertices are:");
 		for (int i = 0; i < vertices.size(); i++) {
//			if (vertices.get(i) != null) {
			if((vertices.get(i)!=null) && (vertices.get(i).getID()>=0)){
				//System.out.println("V[" + i + "]=" + vertices[i] + ",WD= " + vertices[i].getWeightedDegree()); // in case we use the arraylist of wighted degree we can use this method
				System.out.println("V[" + i + "]=" + vertices.get(i) + ",WD= " + vertices.get(i).get_WeightedDegree());// using double value 
			}
		}
	}

	
	
	// let us print all the edgs available in the graph
	public void edgesPrint() {
		
//		Collections.sort(edges);//  using collection is calling of compare to method in Edge class the main thing to sort the edges
		System.out.println("The Edges with their weights  are:");
		for (Edge ee : edges) { // to print the array list of objects of edges
			System.out.println("Edg = " + ee);// using toString to print it;
		}
	}
	
//	public  void printEdges(int M) {
//		//System.out.print(M);
//		for (int i = 0; i < M; i++) {
//			System.out.println(edges.get(i));
//		}
//	}
	
	public  void printEdges(int num) {
		if((edges==null) || (edges.size()==0)) {
			System.out.println("There are no edges: please make sure to have loaded a graph first.");
			return;
		} else {
			//System.out.println("The edges in the term graph are:");
			int min = (num < edges.size()) ? num : edges.size();
			Edge e;
	 		for (int i = 0; i < min; i++) {
				e = edges.get(i);
				if (e != null) {
	 				//System.out.println("edge[" + i + "]: " + e);
				} else {
	 				//System.out.println("edge[" + i + "]: NULL");
				}
			}
	 		if(num < edges.size())
				System.out.println("...");
	 		else
	 			System.out.println();
		}
	}
	
	public  void printEdges() {
		for (int i = 0; i < edges.size(); i++) {
			System.out.println(i+"] "+edges.get(i));
		}
	}

	public void maxMinWeights() {

		// get user with Max Edge.
		Edge ee = edges.stream().max(Comparator.comparing(Edge::getWeight)) //
				.get();
		System.out.println("the Max edge is" + ee);
	}

	// Returns an edge with a specific index
	public Edge getEdge(int i) {
      
		if (edges.get(i)!= null){
     		return edges.get(i);
		}
     		else
    			throw new IllegalStateException();
     }

	// Returns all edge in the graph.
	public ArrayList<Edge> getEdges() {
		return edges;
	}

	 
	public void setEdges(ArrayList<Edge> newEdges) {
		edges = newEdges;
	}

	public int getTotalEdges() {
		return edges.size();
	}

	public Iterator<Edge> getEdgeIterator() {
		return edges.iterator();
	}
	
	/** Weighted core decomposition method:<br>
	 * at each iteration step, remove the node with the minimum weighted degree
	 * Removing vertices and edge 
	 */
	public void wCoreDecomp(int k) {
		
		// k is the number of vertices, which the graph has to have after the deleting processing

		int numVertices = vertices.size();
		
		Edge remEdge;
		
		while(numVertices > k){ // the number of vertices will be decreased after the deleting
			
			double minWdegree = Double.MAX_VALUE; // to associate the first value as min then will compare with the rest of values,will continue till the number of vertices will be k
			Vertex minWdegreeV = null;
			int minWDVidx = -1;
            
			//find the minWdegreeV
			for (int i = 0; i < vertices.size(); i++) {
				if ((vertices.get(i)!=null) && (vertices.get(i).getID()>=0) && (vertices.get(i).get_WeightedDegree() < minWdegree)) { //this is the new current minWdegreeV
					minWdegreeV = vertices.get(i); 
					minWDVidx = i;
					minWdegree = vertices.get(i).get_WeightedDegree();	
				}
			}
			System.out.println("removing V[" + minWDVidx + "]=" + vertices.get(minWDVidx) + "   " + vertices.get(minWDVidx).get_WeightedDegree());
//			vertices[minWDVidx] = null; // As the matrix of vertices is static so that mean we cannot delete the elements from it, we will replace the value we deleting with null
			vertices.remove(minWDVidx);
				
			for(Vertex v: minWdegreeV.getNeighbourVertices()){ // loop during the neighbours of the vertex
				
				remEdge= v.removeVertex(minWdegreeV);
		        
				 //Removing Edges
				 edges.remove(remEdge);
				 //System.out.println("Removing Edge " + remEdge);
			}

			numVertices--;
			
		}
				
	}
	
	public void print_wDegrees() {

		// we have used nullsFirst(Comparator.naturalOrder() coz vertices array
		// has null values
		// as it is static array and when we delete the elements we replace them
		// with null

		// Arrays.sort(vertices,Comparator.nullsFirst(Comparator.naturalOrder()));
//		Arrays.sort(vertices, Comparator.nullsLast(Comparator.naturalOrder()));

		for (int i = 0; i < vertices.size(); i++) {

//			if (vertices.get(i) != null) {
			if((vertices.get(i)!=null) && (vertices.get(i).getID()>=0)){
				// vertices[i].getID();
				System.out.println("vertex<" + i + "> getWeightedDegree = " + vertices.get(i).get_WeightedDegree());
			}
		}
	}
	
	
     // Sorting Edges Using Collection
	public void sort_Edges() {
//		Collections.sort(edges); //nop
		Collections.sort(edges, new EdgeComparator_byWeight());
		int order = 0;
		for (Edge e : edges) {
			e.setOrder(order);
			order++;
		}
	}

	 //using  the heap edges to sort the edges 
//	
//	public void sortEdges() {
//		EdgeHeap ehSort = new EdgeHeap(edges.size());
//		while (edges.size() > 0) {
//			ehSort.insert(edges.remove(edges.size() - 1));// we have to insert edge.e)// removing the last element
//		}
//
//		Edge e = null;
//
//		System.out.println("Sorting list started");
//		int i = 0;
//		while (!ehSort.empty()) {
//			e = ehSort.remove();
//			 
//			    edges.add(e);
//			    e.setOrder(i);
//		    
//			i++;
//		}
//	}

	
	 // Using the heap vertices to sort the vertices
	
	public void sortVertices() {
		 // Using Collections.addAll() method 
//         Collections.addAll(vertices_List, vertices);
         for(Vertex v: vertices)
        	 vertices_List.add(v);
         
         /* Collections.sort method is sorting the 
        	elements of ArrayList in ascending order. */
       
         /* Collections.reverseOrder method is Sorting an ArrayList in descending order */
//         Collections.sort(vertices_List, new VertexComparator_byWeight()); //Collections.reverseOrder());
         Collections.sort(vertices_List, new VertexComparator_byWeight_reverse()); //Collections.reverseOrder());
        
        //******Using Heap for sorting Vertices********//
		/*VerticesHeap ehSort = new VerticesHeap(vertices.length);
		while (vertices_List.size() > 0) {
 			ehSort.insert(vertices_List.remove(vertices_List.size()-1));
		}

		Vertex v = null;

		System.out.println("Sorting list started");
		int i = 0;
		while (!ehSort.empty()) {
			v = ehSort.remove();
			vertices_List.add(v);
			i++;
		}
		printAllVertices(); // to print the vertices before  ording **/
	}
	
	// sorting vertices by using Map, we have saved thr vertices along with their weighted degree in the map then we have sort this map
	public void sort_Vertice() {

		for (Vertex v : vertices) {
			wDM.put(v.getID(), v.get_WeightedDegree());
		}

		Map<Integer, Double> sorted = wDM
				.entrySet()
				.stream()
				.sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
				.collect(
						toMap(Map.Entry::getKey, Map.Entry::getValue,
								(e1, e2) -> e2, LinkedHashMap::new));

		wDM = sorted;
		printAllVertices(); // to print the vertices before  ording 
 	}
	
	public void print_sortedVertices(int num) {
		if((vertices_List==null) || (vertices_List.size()==0)) {
			System.out.println("There are no vertices: please make sure to have loaded a graph first.");
			return;
		} else {
			//System.out.println("The sorted vertices in the term graph are:");
			int min = (num < vertices_List.size()) ? num : vertices_List.size();
			Vertex v;
	 		for (int i = 0; i < min; i++) {
				v = vertices_List.get(i);
				if ((v != null) && (v.getID()>=0)) {
	 				//System.out.println("V[" + i + "] = " + v + ", WD = " + v.get_WeightedDegree()+", degree = "+v.getDegree());
				} else {
	 				//System.out.println("V[" + i + "] = NULL");
				}
			}
	 		if(num < vertices_List.size())
				System.out.println("...");
	 		else
	 			System.out.println();
		}
	}

	public void print_sortedVertices() {
		//we have sorted the vertices using heap and we have create vertices_List to deal with vertices. so what we have done we sorted Verties aeeay and we have sorted then inside arraylist so will be easier to deal with
		System.out.println("The sorted vertices in the term graph are:");
 		for (int i = 0; i < vertices_List.size(); i++) {
 			System.out.println("V[" + i + "] = " + vertices_List.get(i)+ ", WD = " + vertices_List.get(i).get_WeightedDegree());// using double value 
		}
	}
	
	public void printVertices() { //not sorted
		System.out.println("The vertices in the term graph are:");
 		for (int i = 0; i < vertices.size(); i++) {
//			if (vertices.get(i) != null) {
			if((vertices.get(i)!=null) && (vertices.get(i).getID()>=0)){
 				System.out.println("V[" + i + "] = " + vertices.get(i));
			}
		}
		
		/*System.out.println("The nighbours  are:");
		for (int i = 0; i < vertices_List.size(); i++) {

 				System.out.println("V[" + i + "]=" + vertices_List.get(i)+ ",Neighbours= " + vertices_List.get(i).getNeighbourVertices()+", "+" EdgesList=" + vertices_List.get(i).getEdgeList());// using double value 

		}*/
	}

	public void Print_EV() {

		for (int i = 0; i < vertices.size(); i++) {
			if((vertices.get(i)!=null)&& (vertices.get(i).getID()>=0)){
				System.out.println("getEdgeList====" + vertices.get(i).getEdgeList());
				System.out.println("getEdgesWeight===" + vertices.get(i).getEdgeWeights());
				System.out.println("getNeighbourVertices===" + vertices.get(i).getNeighbourVertices());
				//System.out.println("getDegree===" + vertices[i].getDegree());
				System.out.println("getWeightedDegree===" + vertices.get(i).get_WeightedDegree());
				//System.out.println("getOrder===" + vertices[i].getOrder(i));
	
				System.out.println("*****************************************");
			}
		}
 
 		System.out.println("sort all the Edges in the Graph");
 		//Collections.sort(edges);//  using collection is calling of compare to method in Edge class
		edgesPrint();
		 //System.out.println("get ordering Edges in the Graph" + ee.getOrder());

		maxMinWeights();

	}
	
}
	


