/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graphtheory;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;

/**
 *
 * @author mk
 */
public class GraphProperties extends JPanel {
    static int time;

    public int[][] adjacencyMatrix;
    public int[][] distanceMatrix;
    public Vector<VertexPair> vpList;
    public List<List<Integer>> bridges = new ArrayList();
    private static GraphProperties instance;
    private Map<Integer, Double> degreeDistribution; // Stores degree frequency
    
    
    public int[][] generateAdjacencyMatrix(Vector<Vertex> vList, Vector<Edge> eList) {
        adjacencyMatrix = new int[vList.size()][vList.size()];

        for (int a = 0; a < vList.size(); a++)//initialize
        {
            for (int b = 0; b < vList.size(); b++) {
                adjacencyMatrix[a][b] = 0;
            }
        }

        for (int i = 0; i < eList.size(); i++) {
            adjacencyMatrix[vList.indexOf(eList.get(i).vertex1)][vList.indexOf(eList.get(i).vertex2)] = 1;
            adjacencyMatrix[vList.indexOf(eList.get(i).vertex2)][vList.indexOf(eList.get(i).vertex1)] = 1;
        }

        // Populate adjacency matrix with edge weights
        for (Edge edge : eList) {
            int index1 = vList.indexOf(edge.vertex1);
            int index2 = vList.indexOf(edge.vertex2);
            if (edge.weight != null) { // If the edge has a weight
                adjacencyMatrix[index1][index2] = edge.weight;
                adjacencyMatrix[index2][index1] = edge.weight; // Undirected graph
            } else {
                adjacencyMatrix[index1][index2] = 1; // Default weight for unweighted edges
                adjacencyMatrix[index2][index1] = 1; // Undirected graph
            }
        }

        getAllCutpoints(adjacencyMatrix, vList.size());
        computeDegreeDistribution(vList); // Compute degree distribution
        return adjacencyMatrix;
    }

   
    private void computeDegreeDistribution(Vector<Vertex> vList) {
        Map<Integer, Integer> degreeCount = new HashMap<>();
        int totalNodes = vList.size();

        // Count degrees
        for (int i = 0; i < totalNodes; i++) {
            int degree = 0;
            for (int j = 0; j < totalNodes; j++) {
                if (adjacencyMatrix[i][j] > 0) {
                    degree++;
                }
            }
            degreeCount.put(degree, degreeCount.getOrDefault(degree, 0) + 1);
        }

        // Convert to fraction
        degreeDistribution = new TreeMap<>();
        for (Map.Entry<Integer, Integer> entry : degreeCount.entrySet()) {
            degreeDistribution.put(entry.getKey(), entry.getValue() / (double) totalNodes);
        }

        // Display the bar plot
        displayDegreeDistribution();
    }

    private void displayDegreeDistribution() {
        JFrame frame = new JFrame("Degree Distribution");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.setSize(600, 400);
        frame.add(this);
        frame.setVisible(true);
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);

        if (degreeDistribution == null) return;

        Graphics2D g2d = (Graphics2D) g;
        int width = getWidth();
        int height = getHeight();
        int margin = 50;
        DecimalFormat df = new DecimalFormat("0.00"); // Format for two decimal places

        // Define scale
        int maxDegree = degreeDistribution.keySet().stream().max(Integer::compare).orElse(1);
        double maxFraction = degreeDistribution.values().stream().max(Double::compare).orElse(1.0);

        // Draw axes
        g2d.setColor(Color.BLACK);
        g2d.drawLine(margin, height - margin, width - margin, height - margin); // X-axis
        g2d.drawLine(margin, margin, margin, height - margin); // Y-axis

        // Draw bars
        int barWidth = (width - 2 * margin) / (maxDegree + 1);
        int i = 0;

        for (Map.Entry<Integer, Double> entry : degreeDistribution.entrySet()) {
            int degree = entry.getKey();
            double fraction = entry.getValue();

            int barHeight = (int) ((fraction / maxFraction) * (height - 2 * margin));
            int x = margin + i * barWidth;
            int y = height - margin - barHeight;

            g2d.setColor(Color.BLUE);
            g2d.fillRect(x, y, barWidth - 5, barHeight);
            g2d.setColor(Color.BLACK);
            g2d.drawString(String.valueOf(degree), x + barWidth / 4, height - margin + 15); // X-axis labels

            // Draw fraction labels on top of bars
            g2d.drawString(df.format(fraction), x + barWidth / 4, y - 5);

            i++;
        }

        // Label X-axis
        g2d.drawString("Degree", width / 2, height - 10);

        // Rotate and Label Y-axis (Fraction of Nodes)
        AffineTransform originalTransform = g2d.getTransform();
        g2d.rotate(-Math.PI / 2); // Rotate 90 degrees counterclockwise
        g2d.drawString("Fraction of Nodes", -height / 2, margin - 10);
        g2d.setTransform(originalTransform); // Restore original transformation

        // Y-axis labels for fraction values
        int numLabels = 5;
        for (int j = 0; j <= numLabels; j++) {
            double fractionLabel = maxFraction * j / numLabels;
            int yLabel = height - margin - (int) ((fractionLabel / maxFraction) * (height - 2 * margin));
            g2d.drawString(df.format(fractionLabel), 10, yLabel);
        }
    }


    static boolean[] getAllCutpoints(int[][] adjacencyMatrix, int vSize){

        if (vSize < 3){
            boolean[] temp = new boolean[vSize];
            for (int i = 0; i < vSize; i++){
                temp[i] = false;
            }
            return temp;
        }

        ArrayList<ArrayList<Integer> > adj = new ArrayList<ArrayList<Integer> >(vSize);

        for (int i = 0; i < vSize; i++) {
            adj.add(new ArrayList<Integer>());
        }


        for (int i = 1; i < vSize; i++){
            for (int j = 0; j < i; j++){
                if (adjacencyMatrix[i][j] > 0){
                    addEdge(adj, i, j);
                }
            }
        }

        int[] visited = new int[vSize];
        int[] insertion_time = new int[vSize];
        int[] minimum_insertion = new int[vSize];

        int count=1;
        for(int i = 0; i < vSize;i++){
            if(visited[i]==0){
                getBridge(i,-1,visited,insertion_time,minimum_insertion,count,adj);
            }
        }

        System.out.println("Articulation points of the graph");

        boolean[] isAP = AP(adj, vSize);

        return isAP;
    }

    static void addEdge(ArrayList<ArrayList<Integer>> adj, int u, int v)
    {
        adj.get(u).add(v);
        adj.get(v).add(u);
    }

    static void getBridge(int node,int parent,int[] visited,int[] insertion_time,int[] minimum_insertion,int count,ArrayList<ArrayList<Integer>> graph){

        visited[node]=1;
        insertion_time[node] = minimum_insertion[node] = count++;
        for(int nbr:graph.get(node)){
            if(nbr==parent) continue;

            if (visited[nbr] == 0) {
                getBridge(nbr, node, visited, insertion_time, minimum_insertion, count, graph);
                minimum_insertion[node] = Math.min(minimum_insertion[node], minimum_insertion[nbr]);

                if (minimum_insertion[nbr] > insertion_time[node]) {
                    System.out.println("Bridge edge is between Node " + nbr + " and Node " + node);
                    // Input nbr and node in a list that contains the pair
                    List<Integer> bridge = new ArrayList<>();
                    bridge.add(Math.min(nbr, node)); // Store smaller node first for consistency
                    bridge.add(Math.max(nbr, node)); // Store larger node second
                    GraphProperties gp = GraphProperties.getInstance();
                    if (!bridgeExists(gp.bridges, bridge)) { // Use the helper function
                        gp.bridges.add(bridge);
                    }

                }
            } else {
                minimum_insertion[node] = Math.min(minimum_insertion[node], insertion_time[nbr]);
            }

            /*
            if(visited[nbr]==0){
                getBridge(nbr,node,visited,insertion_time,minimum_insertion,count,graph);
                minimum_insertion[node]=Math.min(minimum_insertion[node],minimum_insertion[nbr]);
                if(minimum_insertion[nbr]>insertion_time[node]){
                    System.out.println("Bridge edge is between Node "+nbr+" and Node "+node);

                }
            }
            else{
                minimum_insertion[node]=Math.min(minimum_insertion[node],insertion_time[nbr]);
            }

             */
        }
    }

    private static boolean bridgeExists(List<List<Integer>> bridges, List<Integer> bridgeToCheck) {

        for (List<Integer> existingBridge : bridges) {
            if (existingBridge.equals(bridgeToCheck)) {
                return true; // Bridge already exists
            }
        }
        return false; // Bridge does not exist
    }

    static void APUtil(ArrayList<ArrayList<Integer> > adj, int u,
                       boolean visited[], int disc[], int low[],
                       int parent, boolean isAP[])
    {
        // Count of children in DFS Tree
        int children = 0;

        // Mark the current node as visited
        visited[u] = true;

        // Initialize discovery time and low value
        disc[u] = low[u] = ++time;

        // Go through all vertices adjacent to this
        for (Integer v : adj.get(u)) {
            // If v is not visited yet, then make it a child of u
            // in DFS tree and recur for it
            if (!visited[v]) {
                children++;
                APUtil(adj, v, visited, disc, low, u, isAP);

                // Check if the subtree rooted with v has
                // a connection to one of the ancestors of u
                low[u] = Math.min(low[u], low[v]);

                // If u is not root and low value of one of
                // its child is more than discovery value of u.
                if (parent != -1 && low[v] >= disc[u])
                    isAP[u] = true;
            }

            // Update low value of u for parent function calls.
            else if (v != parent)
                low[u] = Math.min(low[u], disc[v]);
        }

        // If u is root of DFS tree and has two or more children.
        if (parent == -1 && children > 1)
            isAP[u] = true;
    }

    static boolean[] AP(ArrayList<ArrayList<Integer> > adj, int V)
    {
        boolean[] visited = new boolean[V];
        int[] disc = new int[V];
        int[] low = new int[V];
        boolean[] isAP = new boolean[V];
        int time = 0, par = -1;

        // Adding this loop so that the
        // code works even if we are given
        // disconnected graph
        for (int u = 0; u < V; u++)
            if (visited[u] == false)
                APUtil(adj, u, visited, disc, low, par, isAP);

        for (int u = 0; u < V; u++)
            if (isAP[u] == true)
                System.out.print(u + " ");
        System.out.println();

        return isAP;
    }

    public int[][] generateDistanceMatrix(Vector<Vertex> vList) {
        distanceMatrix = new int[vList.size()][vList.size()];

        // Initialize distance matrix with adjacency matrix values
        for (int a = 0; a < vList.size(); a++) {
            for (int b = 0; b < vList.size(); b++) {
                if (a == b) {
                    distanceMatrix[a][b] = 0; // Distance to self is 0
                } else if (adjacencyMatrix[a][b] != 0) {
                    distanceMatrix[a][b] = adjacencyMatrix[a][b]; // Direct edge weight
                } else {
                    distanceMatrix[a][b] = Integer.MAX_VALUE; // No direct edge
                }
            }
        }

        // Apply Floyd-Warshall algorithm to compute shortest paths
        for (int k = 0; k < vList.size(); k++) {
            for (int i = 0; i < vList.size(); i++) {
                for (int j = 0; j < vList.size(); j++) {
                    if (distanceMatrix[i][k] != Integer.MAX_VALUE && distanceMatrix[k][j] != Integer.MAX_VALUE) {
                        distanceMatrix[i][j] = Math.min(distanceMatrix[i][j], distanceMatrix[i][k] + distanceMatrix[k][j]);
                    }
                }
            }
        }
        return distanceMatrix;
    }

    public void displayContainers(Vector<Vertex> vList) {
        vpList = new Vector<VertexPair>();
        int[] kWideGraph = new int[10];
        for (int i = 0; i < kWideGraph.length; i++) {
            kWideGraph[i] = -1;
        }

        VertexPair vp;

        for (int a = 0; a < vList.size(); a++) {    // assign vertex pairs
            for (int b = a + 1; b < vList.size(); b++) {
                vp = new VertexPair(vList.get(a), vList.get(b));
                vpList.add(vp);
                int longestWidth = 0;
                System.out.println(">Vertex Pair " + vList.get(a).name + "-" + vList.get(b).name + "\n All Paths:");
                vp.generateVertexDisjointPaths();
                for (int i = 0; i < vp.VertexDisjointContainer.size(); i++) {//for every container of the vertex pair
                    int width = vp.VertexDisjointContainer.get(i).size();
                    Collections.sort(vp.VertexDisjointContainer.get(i), new descendingWidthComparator());
                    int longestLength = vp.VertexDisjointContainer.get(i).firstElement().size();
                    longestWidth = Math.max(longestWidth, width);
                    System.out.println("\tContainer " + i + " - " + "Width=" + width + " - Length=" + longestLength);

                    for (int j = 0; j < vp.VertexDisjointContainer.get(i).size(); j++) //for every path in the container
                    {
                        System.out.print("\t\tPath " + j + "\n\t\t\t");
                        for (int k = 0; k < vp.VertexDisjointContainer.get(i).get(j).size(); k++) {
                            System.out.print("-" + vp.VertexDisjointContainer.get(i).get(j).get(k).name);
                        }
                        System.out.println();
                    }

                }
                //d-wide for vertexPair
                for (int k = 1; k <= longestWidth; k++) { // 1-wide, 2-wide, 3-wide...
                    int minLength = 999;
                    for (int m = 0; m < vp.VertexDisjointContainer.size(); m++) // for each container with k-wide select shortest length
                    {
                        minLength = Math.min(minLength, vp.VertexDisjointContainer.get(m).size());
                    }
                    if (minLength != 999) {
                        System.out.println(k + "-wide for vertexpair(" + vp.vertex1.name + "-" + vp.vertex2.name + ")=" + minLength);
                        kWideGraph[k] = Math.max(kWideGraph[k], minLength);
                    }
                }
            }
        }

        for (int i = 0; i < kWideGraph.length; i++) {
            if (kWideGraph[i] != -1) {
                System.out.println("D" + i + "(G)=" + kWideGraph[i]);
            }
        }


    }

    public static GraphProperties getInstance() {
        if (instance == null) {
            instance = new GraphProperties();
        }
        return instance;
    }

    public void drawAdjacencyMatrix(Graphics g, Vector<Vertex> vList, int x, int y) {
        int cSize = 20;
        g.setColor(Color.LIGHT_GRAY);
        g.fillRect(x, y - 30, vList.size() * cSize + cSize, vList.size() * cSize + cSize);
        g.setColor(Color.black);
        g.drawString("AdjacencyMatrix", x, y - cSize);
        for (int i = 0; i < vList.size(); i++) {
            g.setColor(Color.RED);
            g.drawString(vList.get(i).name, x + cSize + i * cSize, y);
            g.drawString(vList.get(i).name, x, cSize + i * cSize + y);
            g.setColor(Color.black);
            for (int j = 0; j < vList.size(); j++) {
                g.drawString("" + adjacencyMatrix[i][j], x + cSize * (j + 1), y + cSize * (i + 1));
            }
        }

        // Add here the showing of the cutpoints and bridges
        boolean[] isAP = getAllCutpoints(adjacencyMatrix, vList.size());

        g.drawString("Cutpoints", x,  2*cSize + y + (vList.size() * cSize));

        int spacer = 0;
        for (int u = 0; u < vList.size(); u++) {
            if (isAP[u] == true) {
                g.drawString(" " + u, 10 + x + spacer, 2 * cSize + y / 2 + y + (vList.size() * cSize));
                spacer = spacer + 15;
            }
        }

        g.drawString("Bridges", x,  2*cSize + 2*y + (vList.size() * cSize));

        GraphProperties gp = GraphProperties.getInstance();

        List<String> bridgeStrings = new ArrayList<>();
        for (List<Integer> bridge : gp.bridges) {
            String bridgeString = String.format("(%d-%d)", bridge.get(0), bridge.get(1));
            bridgeStrings.add(bridgeString);
        }
        int space = 0;
        for (int i = 0; i < bridgeStrings.size(); i++){
            g.drawString(bridgeStrings.get(i),10 + x + space,  2*cSize + 2*y + y/2 + (vList.size() * cSize));
            space = space + 35;
        }

    }

    public void drawDistanceMatrix(Graphics g, Vector<Vertex> vList, int x, int y) {
        int cSize = 20;
        g.setColor(Color.LIGHT_GRAY);
        g.fillRect(x, y-30, vList.size() * cSize+cSize, vList.size() * cSize+cSize);
        g.setColor(Color.black);
        g.drawString("ShortestPathMatrix", x, y - cSize);
        for (int i = 0; i < vList.size(); i++) {
            g.setColor(Color.RED);
            g.drawString(vList.get(i).name, x + cSize + i * cSize, y);
            g.drawString(vList.get(i).name, x, cSize + i * cSize + y);
            g.setColor(Color.black);
            for (int j = 0; j < vList.size(); j++) {
                g.drawString("" + distanceMatrix[i][j], x + cSize * (j + 1), y + cSize * (i + 1));
            }
        }
    }

    public Vector<Vertex> vertexConnectivity(Vector<Vertex> vList) {
        Vector<Vertex> origList = new Vector<Vertex>();
        Vector<Vertex> tempList = new Vector<Vertex>();
        Vector<Vertex> toBeRemoved = new Vector<Vertex>();
        Vertex victim;


        origList.setSize(vList.size());
        Collections.copy(origList, vList);

        int maxPossibleRemove = 0;
        while (graphConnectivity(origList)) {
            Collections.sort(origList, new ascendingDegreeComparator());
            maxPossibleRemove = origList.firstElement().getDegree();

            for (Vertex v : origList) {
                if (v.getDegree() == maxPossibleRemove) {
                    for (Vertex z : v.connectedVertices) {
                        if (!tempList.contains(z)) {
                            tempList.add(z);
                        }
                    }
                }
            }

            while (graphConnectivity(origList) && tempList.size() > 0) {
                Collections.sort(tempList, new descendingDegreeComparator());
                victim = tempList.firstElement();
                tempList.removeElementAt(0);
                origList.remove(victim);
                for (Vertex x : origList) {
                    x.connectedVertices.remove(victim);
                }
                toBeRemoved.add(victim);
            }
            tempList.removeAllElements();
        }

        return toBeRemoved;
    }

    private boolean graphConnectivity(Vector<Vertex> vList) {

        Vector<Vertex> visitedList = new Vector<Vertex>();

        recurseGraphConnectivity(vList.firstElement().connectedVertices, visitedList); //recursive function
        if (visitedList.size() != vList.size()) {
            return false;
        } else {
            return true;
        }
    }

    private void recurseGraphConnectivity(Vector<Vertex> vList, Vector<Vertex> visitedList) {
        for (Vertex v : vList) {
            {
                if (!visitedList.contains(v)) {
                    visitedList.add(v);
                    recurseGraphConnectivity(v.connectedVertices, visitedList);
                }
            }
        }
    }

    private class ascendingDegreeComparator implements Comparator {

        public int compare(Object v1, Object v2) {

            if (((Vertex) v1).getDegree() > ((Vertex) v2).getDegree()) {
                return 1;
            } else if (((Vertex) v1).getDegree() > ((Vertex) v2).getDegree()) {
                return -1;
            } else {
                return 0;
            }
        }
    }

    private class descendingDegreeComparator implements Comparator {

        public int compare(Object v1, Object v2) {

            if (((Vertex) v1).getDegree() > ((Vertex) v2).getDegree()) {
                return -1;
            } else if (((Vertex) v1).getDegree() > ((Vertex) v2).getDegree()) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    private class descendingWidthComparator implements Comparator {

        public int compare(Object v1, Object v2) {

            if (((Vector<Vertex>) v1).size() > (((Vector<Vertex>) v2).size())) {
                return -1;
            } else if (((Vector<Vertex>) v1).size() < (((Vector<Vertex>) v2).size())) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}