package ai;

public class Node {
	char node;
	int searchCost;
	int x;
	int y;
	int pathCost;
	int[] parent;

	Node(char node, int searchCost, int x, int y, int pathCost, int[] parent) {
		this.node = node;
		this.searchCost = searchCost;
		this.x = x;
		this.y = y;
		this.pathCost = pathCost;
		this.parent = parent;
	}

}
