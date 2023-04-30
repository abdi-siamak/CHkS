package uk.ac.rdg.chks.bb;

import uk.ac.rdg.chks.graph.Vertex;

import java.util.*;

public class Top_Solutions {
    public static int topSize; // the size of top K solutions
    private static double similarity;
    private static ArrayList<BnB_Node> solutionList= new ArrayList<>();
    private static LinkedHashMap<BnB_Node, List<BnB_Node>> childSolution= new LinkedHashMap<>();
    public static void insertSolution (BnB_Node newSolution){
        /*
        for (BnB_Node s:solutionList){
            System.out.println(s.getWeight());
        }
            System.out.println(----------------------);
         */

        if(!solutionList.isEmpty()) {
            for(BnB_Node solution:solutionList){
                if(newSolution.getWeight() >= solution.getWeight()){ // a new solution should be inserted
                    if (similarity(newSolution.getVertices(), solution.getVertices()) >= similarity){ // old solution should be moved into a child list of this slot
                        //System.out.println("similar");
                        //show(newSolution.getVertices(), solution.getVertices());
                        //System.out.println("new: " + newSolution.getWeight());
                        //System.out.println("old: " + solution.getWeight()+"\n");
                        if (!childSolution.containsKey(solution)){
                            ArrayList<BnB_Node> oldSolution = new ArrayList<>();
                            oldSolution.add(solution);
                            childSolution.put(newSolution, oldSolution); // track the child solutions with the object newSolution everytime
                        }else {
                            ArrayList<BnB_Node> oldSolution = new ArrayList<>();
                            oldSolution.add(solution);
                            oldSolution.addAll(childSolution.get(solution));
                            childSolution.remove(solution);
                            childSolution.put(newSolution, oldSolution);
                        }
                        solutionList.set(solutionList.indexOf(solution), newSolution);
                    }else if(similarity(newSolution.getVertices(), solution.getVertices()) < similarity){ // the new solution should be added and the old one should be shifted to the bottom
                        //System.out.println("not similar");
                        solutionList.add(solutionList.indexOf(solution), newSolution); // this command shift the old one to the bottom
                        if (solutionList.size()>topSize){ // if the solution list exceeds its size, remove from the bottom
                            solutionList.remove(topSize);
                        }
                    }
                    break;
                } else if (newSolution.getWeight() < solution.getWeight()) {
                    break;
                }
            }
        }else{
            solutionList.add(newSolution);
        }
    }

    private static void show (ArrayList<Vertex> vertices_1, ArrayList<Vertex> vertices_2) {
        Set<String> set1 = new HashSet<>();
        Set<String> set2 = new HashSet<>();
        for (Vertex v : vertices_1) {
            set1.add(CHkS.dictionary.get(v.getID()));
        }
        for (Vertex v : vertices_2) {
            set2.add(CHkS.dictionary.get(v.getID()));
        }
        System.out.println("new solutions: " + set1);
        System.out.println("old solutions: " + set2);
    }

    private static float similarity(ArrayList<Vertex> vertices_1, ArrayList<Vertex> vertices_2) {
        float union;
        float intersection;
        Set<Integer> set1 = new HashSet<>();
        Set<Integer> set2 = new HashSet<>();
        for (Vertex v : vertices_1) {
            set1.add(v.getID());
        }
        for (Vertex v : vertices_2) {
            set2.add(v.getID());
        }
        set1.retainAll(set2);
        intersection = set1.size();
        Set<Integer> un = new HashSet<>(set1);
        un.addAll(set2);
        union = un.size();
        //System.out.println(set1);
        //System.out.println(set2);
        //System.out.println(intersection/union);
        return intersection/union;
    }
    public static ArrayList<BnB_Node> getSolutions (){
        return solutionList;
    }
    public static LinkedHashMap<BnB_Node, List<BnB_Node>> getChildSolutions (){
        return childSolution;
    }
    public static void setSimilarity(double sim){
        similarity = sim;
    }
}
