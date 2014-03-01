package ai;

import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Scanner;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class tsp {

	protected static ArrayList<char[]> map = new ArrayList<char[]>();;
	protected static TreeSet<Character> allCheckPoints = new TreeSet<Character>();
	private static final Pattern pattern = Pattern.compile("[A-Za-z]");
	public static int ROW_MAX = 0;
	public static int COL_MAX = 0;
	private static final Comparator<Node> my_total_order = new myComp();
	private static final Comparator<State> state_comparator = new stateComp();

	private static class stateComp implements Comparator<State> {
		public int compare(State s, State t) {
			int[] s_xy = new int[2];
			getIndex(s.current_checkpoint, s_xy);
			int[] t_xy = new int[2];
			getIndex(t.current_checkpoint, t_xy);
			if (s.f == t.f) {
				int i = 2;
				while (s_xy[1] == t_xy[1] && s_xy[0] == t_xy[0]) {
					if (s.visited_checkpoints.get(s.visited_checkpoints.size()
							- i) >= 0
							&& t.visited_checkpoints.get(t.visited_checkpoints
									.size() - i) > 0) {
						getIndex(
								s.visited_checkpoints.get(s.visited_checkpoints
										.size() - i), s_xy);
						getIndex(
								t.visited_checkpoints.get(t.visited_checkpoints
										.size() - i), t_xy);
					} else {
						break;
					}
					i++;
				}

				if (s_xy[0] == t_xy[0]) {
					return ((Integer) s_xy[1]).compareTo(t_xy[1]);
				} else {
					return ((Integer) s_xy[0]).compareTo(t_xy[0]);
				}

			} else {
				if (s.f < t.f) {
					return -1;
				} else {
					return 1;
				}
				// return ((Integer) s.f).compareTo(t.f);
			}

		}
	}

	private static class myComp implements Comparator<Node> {

		public int compare(Node s, Node t) {
			/*
			 * return (s.searchCost == t.searchCost) ? ((s.x == t.x) ?
			 * ((Integer) s.y) .compareTo(t.y) : ((Integer) s.x).compareTo(t.x))
			 * : (((Integer) s.searchCost).compareTo(t.searchCost));
			 */

			if (s.searchCost == t.searchCost) {
				if (s.x == t.x) {
					return ((Integer) s.y).compareTo(t.y);
				} else {
					return ((Integer) s.x).compareTo(t.x);
				}

			} else {
				if (s.searchCost < t.searchCost) {
					return -1;
				} else {
					return 1;
				}
			}

		}

	}

	private static class Node {
		char node;
		double searchCost;
		int x;
		int y;
		int pathCost;
		int[] parent = new int[2];

		Node(char node, double searchCost, int x, int y, int pathCost,
				int[] parent) {
			this.node = node;
			this.searchCost = searchCost;
			this.x = x;
			this.y = y;
			this.pathCost = pathCost;
			this.parent[0] = parent[0];
			this.parent[1] = parent[1];
		}
	}

	public static void main(String[] args) throws IOException {

		int taskNumber = 0;
		String inputFile = null;

		String outputFile = null;
		String outputLog = null;

		boolean firstT = true;

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-t") && firstT) {
				taskNumber = Integer.parseInt(args[i + 1]);
				firstT = false;
			}
			if (args[i].equals("-i")) {
				inputFile = args[i + 1];
			}
			if (args[i].equals("-op")) {
				outputFile = args[i + 1];
			} else if (args[i].equals("-ol")) {
				outputLog = args[i + 1];
			}
		}
		String str = null;
		Scanner s = null;
		FileReader fr = null;
		try {
			fr = new FileReader(inputFile);
			s = new Scanner(fr);
			while (s.hasNextLine()) {
				str = s.nextLine();
				getCheckPoints(str);
				map.add(str.toCharArray());
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			s.close();
			try {
				fr.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		ROW_MAX = map.size() - 1;
		COL_MAX = map.get(0).length - 1;

		if (taskNumber == 1) {
			shortestPath(outputFile, outputLog, false);
		} else if (taskNumber == 2) {
			tsp1(outputFile, outputLog);
		}

	}

	private static void getCheckPoints(String str) {

		Matcher matcher = pattern.matcher(str);
		while (matcher.find()) {
			allCheckPoints.add(str.charAt(matcher.start()));
		}

	}

	public static HashMap<String, Double> shortestPath(String outputFileName,
			String outputLogName, boolean fromtsp) throws IOException {

		HashMap<String, Double> allEdges = new HashMap<String, Double>();
		Character[] checkpoints = allCheckPoints
				.toArray(new Character[allCheckPoints.size()]);
		Arrays.sort(checkpoints);
		double edgeCost = 0.0;
		BufferedWriter outputFile = null;
		BufferedWriter outputLog = null;

		outputFile = new BufferedWriter(new FileWriter(outputFileName));
		outputLog = new BufferedWriter(new FileWriter(outputLogName));
		for (int i = 0; i < checkpoints.length; i++) {
			for (int j = i + 1; j < checkpoints.length; j++) {
				// System.out.print(checkpoints[i] + " and " + checkpoints[j]);
				edgeCost = AStarsearch(checkpoints[i], checkpoints[j],
						outputLog);
				// System.out.print(" " + edgeCost);
				// System.out.println();
				allEdges.put(checkpoints[i] + "," + checkpoints[j], edgeCost);
				if (fromtsp == false) {
					outputFile.write(checkpoints[i] + "," + checkpoints[j]
							+ "," + edgeCost);
					outputFile.newLine();
				}
			}
		}
		outputFile.close();
		outputLog.close();
		return allEdges;
	}

	private static class State {
		char current_checkpoint;
		ArrayList<Character> visited_checkpoints;
		double g;
		double h;
		double f;

		State(char current_checkpoint,
				ArrayList<Character> visited_checkpoints, double g, double h) {
			this.current_checkpoint = current_checkpoint;
			this.visited_checkpoints = visited_checkpoints;
			this.g = g;
			this.h = h;
			this.f = (int) (g + h);
		}
	}

	private static HashMap<String, Double> getSubgraph(
			HashMap<String, Double> allEdges, ArrayList<Character> visited) {

		HashMap<String, Double> subgraph = new HashMap<String, Double>();
		String[] parts = null;
		// adding source to subgraph
		for (Entry<String, Double> edge : allEdges.entrySet()) {
			parts = edge.getKey().split(",");
			if ((!visited.contains(parts[0].charAt(0)) && !visited
					.contains(parts[1].charAt(0)))
					|| (parts[0].charAt(0) == 'A' && !visited.contains(parts[1]
							.charAt(0)))
					|| (parts[1].charAt(0) == 'A' && !visited.contains(parts[0]
							.charAt(0)))) {
				subgraph.put(edge.getKey(), edge.getValue());
			}
		}
		return subgraph;
	}

	public static void tsp1(String outputFileName, String outputLogName)
			throws IOException {

		HashMap<String, Double> allEdges = shortestPath(outputFileName,
				outputLogName, true);
		HashMap<String, Double> subgraph = new HashMap<String, Double>();

		char currentCheckPoint = 'A';
		double h = 0, g = 0;
		double parentPathCost = 0;
		char current;
		BufferedWriter outputFile = new BufferedWriter(new FileWriter(
				outputFileName));
		BufferedWriter outputLog = new BufferedWriter(new FileWriter(
				outputLogName));

		PriorityQueue<State> pq = new PriorityQueue<State>(17, state_comparator);

		ArrayList<Character> temp = new ArrayList<Character>();
		temp.add('A');

		State child;
		h = MST(allEdges);
		State s = new State(currentCheckPoint, temp, g, h);
		pq.add(s);
		parentPathCost = 0;
		while (!pq.isEmpty()) {
			s = pq.remove();
			/*
			 * if (visited.contains(s.current_checkpoint)) { continue; }else{
			 * visited.add(s.current_checkpoint); }
			 */
			 String str = new String();
			 for(char ch : s.visited_checkpoints){
			 	str += Character.toString(ch);
			 }
			outputLog.write(str + "," + s.g + "," + s.h + ","
					+ s.f);
			outputLog.newLine();
			parentPathCost = s.g;
			current = s.current_checkpoint;
			if (isEqual(s.visited_checkpoints, allCheckPoints)) {
				s.visited_checkpoints.add('A');
				child = new State('A', s.visited_checkpoints, s.f, 0);
				str = new String();
				 for(char ch : s.visited_checkpoints){
				 	str += Character.toString(ch);
				 }
				outputLog.write(str + "," + child.g + ","
						+ child.h + "," + child.f);
				outputLog.newLine();
				for (char ch : s.visited_checkpoints) {
					outputFile.write(ch);
					outputFile.newLine();
				}
				outputFile.write("Total Tour Cost: " + child.f);
				outputFile.close();
				outputLog.close();
				return; // solution found
			}
			// visited.add(current);
			// explore neighbors

			subgraph = getSubgraph(allEdges, s.visited_checkpoints);
			h = MST(subgraph);
			for (char innerC : allCheckPoints) {
				if (s.visited_checkpoints != null
						&& s.visited_checkpoints.contains(innerC)) {
					continue;
				}
				if (current < innerC) {
					g = parentPathCost + allEdges.get(current + "," + innerC);
				} else {
					g = parentPathCost + allEdges.get(innerC + "," + current);
				}
				ArrayList<Character> current_visited = new ArrayList<Character>();
				if (s.visited_checkpoints != null) {
					current_visited.addAll(s.visited_checkpoints);
				}
				current_visited.add(innerC);
				child = new State(innerC, current_visited, g, h);
				pq.add(child);
			}

		}
		outputFile.close();
		outputLog.close();
	}

	public static double AStarsearch(char startNode, char endNode,
			BufferedWriter outputLog) throws IOException {
		int[] s = new int[2];
		getIndex(startNode, s);
		int[] t = new int[2];
		getIndex(endNode, t);
		int[] intermediate = new int[2];
		int[] c = { 0, 0 };
		char nextNode;
		Node from = null;
		boolean[][] visited_coordinates = new boolean[ROW_MAX + 1][COL_MAX + 1];
		PriorityQueue<Node> pq = new PriorityQueue<Node>(17, my_total_order);
		pq.add(new Node(startNode, getManhattanDist(s, t), s[0], s[1], 0, c));
		/*
		 * from 'A' to 'B' -----------------------------------------------
		 */
		outputLog.write("from '" + startNode + "' to '" + endNode + "'");
		outputLog.newLine();
		outputLog.write("-----------------------------------------------");
		outputLog.newLine();
		outputLog.write("x,y,g,h,f");
		outputLog.newLine();
		while (!pq.isEmpty()) {
			from = pq.remove();
			c[0] = from.x;
			c[1] = from.y;
			if (visited_coordinates[c[0]][c[1]] == true) {
				continue;
			} else {
				visited_coordinates[c[0]][c[1]] = true;
			}
			outputLog
					.write(from.y + "," + from.x + "," + from.pathCost + ","
							+ (from.searchCost - from.pathCost) + ","
							+ from.searchCost);
			outputLog.newLine();

			// System.out.println("("+from.x+","+from.y+") [parent("+from.parent[0]+","+from.parent[1]+")]");
			if (from.node == endNode) {
				outputLog
						.write("-----------------------------------------------");
				outputLog.newLine();
				return from.searchCost; // solution found
			}
			// U
			if (!(c[0] - 1 < 0 || c[0] - 1 > ROW_MAX)
					&& !(c[1] < 0 || c[1] > COL_MAX)) {
				nextNode = map.get(c[0] - 1)[c[1]];
				if (nextNode != '*') {
					intermediate[0] = c[0] - 1;
					intermediate[1] = c[1];
					pq.add(new Node(nextNode, 1 + from.pathCost
							+ getManhattanDist(intermediate, t),
							intermediate[0], intermediate[1],
							1 + from.pathCost, c));
				}
			}
			// R
			if (!(c[0] < 0 || c[0] > ROW_MAX)
					&& !(c[1] + 1 < 0 || c[1] + 1 > COL_MAX)) {
				nextNode = map.get(c[0])[c[1] + 1];
				if (nextNode != '*') {
					intermediate[0] = c[0];
					intermediate[1] = c[1] + 1;
					pq.add(new Node(nextNode, 1 + from.pathCost
							+ getManhattanDist(intermediate, t),
							intermediate[0], intermediate[1],
							1 + from.pathCost, c));
				}
			}
			// D
			if (!(c[0] + 1 < 0 || c[0] + 1 > ROW_MAX)
					&& !(c[1] < 0 || c[1] > COL_MAX)) {
				nextNode = map.get((c[0] + 1))[c[1]];
				if (nextNode != '*') {
					intermediate[0] = c[0] + 1;
					intermediate[1] = c[1];
					pq.add(new Node(nextNode, 1 + from.pathCost
							+ getManhattanDist(intermediate, t),
							intermediate[0], intermediate[1],
							1 + from.pathCost, c));
				}
			}
			// L
			if (!(c[0] < 0 || c[0] > ROW_MAX)
					&& !(c[1] - 1 < 0 || c[1] - 1 > COL_MAX)) {
				nextNode = map.get(c[0])[c[1] - 1];
				if (nextNode != '*') {
					intermediate[0] = c[0];
					intermediate[1] = c[1] - 1;
					pq.add(new Node(nextNode, 1 + from.pathCost
							+ getManhattanDist(intermediate, t),
							intermediate[0], intermediate[1],
							1 + from.pathCost, c));

				}
			}
		}
		return 0;
	}

	private static int getManhattanDist(int[] s, int[] t) {

		int diff_x = s[0] - t[0];
		int diff_y = s[1] - t[1];
		if (diff_x < 0) {
			diff_x *= -1;
		}
		if (diff_y < 0) {
			diff_y *= -1;
		}
		return diff_x + diff_y;
	}

	public static void getIndex(char node, int[] coordinates) {

		int i = 0, j = 0;
		for (; i <= ROW_MAX; i++) {
			for (j = 0; j <= COL_MAX; j++) {
				if (node == map.get(i)[j]) {
					coordinates[0] = i;
					coordinates[1] = j;
					return;
				}
			}
		}
	}

	private static HashMap<String, Double> getNeighbours(
			HashMap<String, Double> graph, char currentCheckPoint) {
		HashMap<String, Double> neighbours = new HashMap<String, Double>();
		String[] str = null;
		for (String c : graph.keySet()) {
			str = c.split(",");
			if (str[0].charAt(0) == currentCheckPoint) {
				neighbours.put(c, graph.get(c));
			} else if (str[1].charAt(0) == currentCheckPoint) {
				neighbours.put(str[1] + "," + str[0], graph.get(c));
			}
		}
		return neighbours;
	}

	private static class nbr {
		String str;
		double cost;

		nbr(String str, double cost) {
			this.str = str;
			this.cost = cost;
		}
	}

	private static final Comparator<nbr> nbr_comparator = new nbrComp();

	private static class nbrComp implements Comparator<nbr> {
		public int compare(nbr s, nbr t) {
			if (s.cost < t.cost) {
				return -1;
			} else if (s.cost > t.cost) {
				return 1;
			} else {
				int[] s_xy = new int[2];
				getIndex(s.str.charAt(2), s_xy);
				int[] t_xy = new int[2];
				getIndex(t.str.charAt(2), t_xy);
				if (s_xy[0] == t_xy[0]) {
					return ((Integer) s_xy[1]).compareTo(t_xy[1]);
				} else {
					return ((Integer) s_xy[0]).compareTo(t_xy[0]);
				}
			}
		}
	}

	public static double MST(HashMap<String, Double> graph) {

		ArrayList<String> mst = new ArrayList<String>();
		HashSet<Character> Vnew = new HashSet<Character>();
		double mstcost = 0.0;
		PriorityQueue<nbr> pq = new PriorityQueue<nbr>(11, nbr_comparator);
		nbr popped = null;
		TreeSet<Character> allPointsInSubgraph = new TreeSet<Character>();
		for (String entry : graph.keySet()) {
			allPointsInSubgraph.add(entry.charAt(0));
			allPointsInSubgraph.add(entry.charAt(2));
		}
		char current = 'A';
		Vnew.add(current);
		while (!isEqual(Vnew, allPointsInSubgraph)) {

			for (Entry<String, Double> v : getNeighbours(graph, current)
					.entrySet()) {
				if (!Vnew.contains(v.getKey().charAt(2))) {
					pq.add(new nbr(v.getKey(), v.getValue()));
				}
			}
			popped = pq.remove();
			while (Vnew.contains(popped.str.charAt(0))
					&& Vnew.contains(popped.str.charAt(2))) {
				popped = pq.remove();

			}
			mstcost += popped.cost;
			mst.add(popped.str + "," + popped.cost);
			current = popped.str.charAt(2);
			Vnew.add(current);
		}

		return mstcost;
	}

	private static boolean isEqual(HashSet<Character> Vnew, TreeSet<Character> V) {

		boolean isEqual = true;
		for (char ch : V) {
			if (!Vnew.contains(ch)) {
				isEqual = false;
				break;
			}
		}
		return isEqual;
	}

	private static boolean isEqual(ArrayList<Character> Vnew,
			TreeSet<Character> V) {

		boolean isEqual = true;
		for (char ch : V) {
			if (!Vnew.contains(ch)) {
				isEqual = false;
				break;
			}
		}
		return isEqual;
	}
}
