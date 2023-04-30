package uk.ac.rdg.chks.graph;

import java.util.Comparator;

public class EdgeComparator_byWeight implements Comparator<Edge> {

	public int compare(Edge e1, Edge e2) {
		return e2.getWeight() - e1.getWeight();			// For Descending order
		// return e1.getWeight() - e2.getWeight();		// For Ascending order
	}
	
}
