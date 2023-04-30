package uk.ac.rdg.chks.graph;

import java.util.Comparator;

public class VertexComparator_byWeight implements Comparator<Vertex> {

	@Override
	public int compare(Vertex v0, Vertex v1) {
		int value = 0;

		if ((v0 != null) && (v1 != null)) {

			if (v0.get_WeightedDegree() > v1.get_WeightedDegree())
				value = 1;
			else if (v0.get_WeightedDegree() == v1.get_WeightedDegree())
				value = 0;
			else
				value = -1;
		}

		return value;
	}

}
