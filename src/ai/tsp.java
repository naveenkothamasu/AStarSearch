package ai;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Scanner;
import java.util.SortedMap;
import java.util.TreeMap;
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
				return ((Integer) s.searchCost).compareTo(t.searchCost);
			}

		}

	}

	private static class Node {
		char node;
		int searchCost;
		int x;
		int y;
		int pathCost;
		int[] parent = new int[2];

		Node(char node, int searchCost, int x, int y, int pathCost, int[] parent) {
			this.node = node;
			this.searchCost = searchCost;
			this.x = x;
			this.y = y;
			this.pathCost = pathCost;
			this.parent[0] = parent[0];
			this.parent[1] = parent[1];
		}
	}

	public static void main(String[] args) {

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
			shortestPath();
		} else if (taskNumber == 2) {
			tsp();
		}

	}

	private static void getCheckPoints(String str) {

		Matcher matcher = pattern.matcher(str);
		while (matcher.find()) {
			allCheckPoints.add(str.charAt(matcher.start()));
		}

	}

	public static ArrayList<String> shortestPath() {

		ArrayList<String> allEdges = new ArrayList<String>();
		char[] checkpoints = allCheckPoints.toString().toCharArray();
		Arrays.sort(checkpoints);
		double edgeCost = 0.0;
		/*
		 * for (int i = 0; i < checkpoints.length; i++) { for (int j = i + 1; j
		 * < checkpoints.length; j++) { System.out.print(checkpoints[i] +
		 * " and " + checkpoints[j]); edgeCost = AStarsearch(checkpoints[i],
		 * checkpoints[j]); System.out.print(" " + edgeCost);
		 * System.out.println(); allEdges.add(checkpoints[i] + "," +
		 * checkpoints[j] + "," + edgeCost); } }
		 */
		System.out.println(AStarsearch('E', 'F'));
		return allEdges;
	}

	private static class State {
		char current_checkpoint;
		List<Character> visited_checkpoints;
		double g;
		double h;

		State(char current_checkpoint, List<Character> visited_checkpoints,
				double g, double h) {
			this.current_checkpoint = current_checkpoint;
			this.visited_checkpoints = visited_checkpoints;
			this.g = g;
			this.h = h;
		}
	}

	public static void tsp() {
		ArrayList<String> allEdges = shortestPath();
		ArrayList<Character> visited = new ArrayList<Character>();
		allCheckPoints.remove(allCheckPoints.first());
		for (char current : allCheckPoints) {

		}
	}

	public static int AStarsearch(char startNode, char endNode) {
		int[] s = new int[2];
		getIndex(startNode, s);
		int[] t = new int[2];
		getIndex(endNode, t);
		int[] intermediate = new int[2];
		int[] c = { 0, 0 };
		char nextNode;
		Node from = null;
		boolean[][] visited_coordinates = new boolean[ROW_MAX + 1][COL_MAX + 1];
		// System.out.print(" (" +s[0] + ","+ s[1]+") to (");
		// System.out.println(t[0] + ","+ t[1] +")");
		PriorityQueue<Node> pq = new PriorityQueue<Node>(17, my_total_order);
		pq.add(new Node(startNode, getManhattanDist(s, t), s[0], s[1], 0, c));
		while (!pq.isEmpty()) {
			from = pq.remove();
			c[0] = from.x;
			c[1] = from.y;
			if(visited_coordinates[c[0]][c[1]] == true){
				continue;
			}else{
				visited_coordinates[c[0]][c[1]] = true;
			}
			System.out.println(from.y + "," + from.x + "," + from.pathCost
					+ "," + (from.searchCost - from.pathCost) + ","
					+ from.searchCost);
			
			//System.out.println("("+from.x+","+from.y+") [parent("+from.parent[0]+","+from.parent[1]+")]");
			if (from.node == endNode) {
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

	public ArrayList<String> MST(ArrayList<String> graph) {

		ArrayList<String> mst = new ArrayList<String>();
		SortedMap<Double, Character[]> sortedEdges = new TreeMap<Double, Character[]>();
		Character[] nodes = null;
		HashSet<Character> Vnew = new HashSet<Character>();
		HashSet<Character> V = new HashSet<Character>();

		for (String entry : graph) {
			String[] str = entry.split(",");
			nodes = new Character[2];
			nodes[0] = str[1].toCharArray()[0];
			nodes[1] = str[2].toCharArray()[0];
			V.add(nodes[0]);
			V.add(nodes[1]);
			sortedEdges.put(Double.parseDouble(str[0]), nodes);
		}
		Vnew.add(nodes[1]);
		Character[] currentNodes = { '1', '2' };
		double edgeCost = 0.0;
		while (!isEqual(Vnew, V)) {
			while (Vnew.contains(currentNodes[0])
					&& !V.contains(currentNodes[1])) {
				edgeCost = sortedEdges.firstKey();
				currentNodes = sortedEdges.remove(edgeCost);
			}
			Vnew.add(currentNodes[1]);
			mst.add(currentNodes[0] + "," + currentNodes[1] + "," + edgeCost);
		}
		return mst;
	}

	private boolean isEqual(HashSet<Character> Vnew, HashSet<Character> V) {

		boolean isEqual = true;
		for (char ch : Vnew) {
			if (!V.contains(ch)) {
				isEqual = false;
				break;
			}
		}
		return isEqual;
	}
}
