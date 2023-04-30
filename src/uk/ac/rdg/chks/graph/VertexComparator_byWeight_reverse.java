package uk.ac.rdg.chks.graph;

import java.util.Comparator;

public class VertexComparator_byWeight_reverse implements Comparator<Vertex> {

	private int m_opt; //0: sort by WeightedDegre (all edges); 1: sort by WeightedDegre for only top edges
	
	/**
	 * Constructor
	 */
	public VertexComparator_byWeight_reverse(){
		m_opt = 0;
	}

	public VertexComparator_byWeight_reverse(int edgesNum){
		m_opt = edgesNum;
	}

	
	@Override
	public int compare(Vertex v1, Vertex v2) {
		int res = 0;
		double val1=0, val2=0;
		if(m_opt == 0) {
			val1 = v1.get_WeightedDegree();
			val2 = v2.get_WeightedDegree();
		} else if (m_opt == 1) {
			val1 = v1.get_topWeightedDegree();
			val2 = v2.get_topWeightedDegree();			
		} else {
			//error: option not available
		}
		
		if ((v1 != null) && (v2 != null)) {

			if (val1 > val2)
				res = -1;
			else if (val1 == val2)
				res = 0;
			else
				res = 1;
		}

		return res;
	}

}
