package uk.ac.rdg.chks.graph;

/**
 * The edge class
 *
 */
public class Edge implements Comparable<Edge>{
	private Vertex source,destination;
	private int weight;
	/*The order is a value to keep track of their index in the sorted edges list*/
	private int order;
	
	public Edge(){
		this.source = null;
		this.destination = null;
		this.weight = 0;
		this.order = -1;
	}
	public Edge(Vertex source, Vertex destination, int weight){
		this.source = source;
		this.destination = destination;
		this.weight = weight;
		this.order = -1;
	}
	
	
	public Vertex getVertex1(){
		return source;
	}
	public Vertex getVertex2(){
		return destination;
	}
	public  int getWeight(){
		return weight;
	}
	
	public int getOrder(){
		return order;
	}
	
	public void setOrder(int order){
		this.order = order;
	}
	
		
	public boolean isIncident( Vertex v){
		if(this.source == v){
			return true;
		}
		if(this.destination == v){
			return true;
		}
		return false;
	}
	
	/*@Override
	public int compareTo(Edge e) {
		if( this.weight > e.weight ) 
			return 1;
		else if ( this.weight == e.weight)
			return 0;
		else
			return -1;
	}*/
	
	@Override
	public int compareTo(Edge e) {
		int res = this.getVertex1().getID() - e.getVertex1().getID();
		if(res == 0)
			res = this.getVertex2().getID() - e.getVertex2().getID();
		return res;
	}
  
	//nop!
//	public int compareTo(Edge e) {
//		int compareEdges_Wiegts = ((Edge) e).getWeight();
//		// For Ascending order
//		// return this.weight-compareage;
//
//		// For Descending order
//		return compareEdges_Wiegts - this.weight;
//	}
	
	
//	public static double inducedWeight(ArrayList<Vertex> subgraphvertices){
//		ArrayList<Vertex> neighbour;
//		
//		double totalWeight = 0;
//		
//		for (int i = 0; i < subgraphvertices.size(); i++) {
//			neighbour = subgraphvertices.get(i).getNeighbourVertices();
//			for (int j = 0; j < neighbour.size(); j++) {
//				if (subgraphvertices.contains(neighbour.get(j))) {
//					totalWeight += subgraphvertices.get(i).getWeight(j);
//				}
//			}
//		}
//		/**
//		 * as the graph is undirect for this reason edge such as (7,3) 
//		 * it will be computed twice once with vertex 7 as 3 its neighbour and 
//		 * once again with vertex 3 as 7 its neighbour  
//		 * */
//		return totalWeight/2;
//	}
//	
//	public static int inducedEdges(ArrayList<Vertex> myList) {
//		ArrayList<Vertex> neigh;
//		int total = 0;
//		for (int i = 0; i < myList.size(); i++) {
//			neigh = myList.get(i).getNeighbourVertices();
//			for (int j = 0; j < neigh.size(); j++) {
//				if (myList.contains(neigh.get(j))) {
//					total++;
//				}
//			}
//		}
//		return total/2;
//	}
	
	
	public String toString(){
		// to print Edges
		return "(<" + source + "," + destination + ">,w=" + weight+")" ;
	}
	
}
