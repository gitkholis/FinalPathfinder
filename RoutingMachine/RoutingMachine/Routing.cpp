#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
#include <random>
using namespace std;

const int INF = 999999999;

double readInputData(vector<vector<pair<int, int>>>& graph, 
				   vector<tuple<int, int, int>>& banned_turns, 
				   int& vertex_count, int& arc_count,
				   vector<pair<long double, long double>>& vertex_coordinates,
				   bool data_with_coordinates) {
	clock_t start_time = clock();
	int banned_turns_count;
	char source_file[] = "nyc_data_meters.txt";
	ifstream fin;
	fin.open(source_file);
	fin >> vertex_count >> arc_count;
	graph.resize(vertex_count + 1);
	for (int i = 1; i <= arc_count; ++i) {
		int from_vertex, to_vertex, arc_weight;
		fin >> from_vertex >> to_vertex >> arc_weight;
		graph[from_vertex].push_back({ to_vertex, arc_weight });
	}
	fin >> banned_turns_count;
	for (int i = 0; i < banned_turns_count; ++i) {
		int from_vertex, by_vertex, to_vertex;
		fin >> from_vertex >>by_vertex >> to_vertex;
		banned_turns.push_back({ from_vertex, by_vertex , to_vertex});
	}
	fin.close();
	if (data_with_coordinates) {
		char source_file_coord[] = "nyc_coordinates.txt";
		fin.open(source_file_coord);
		fin >> vertex_count;
		vertex_coordinates.push_back({ 0 , 0 });
		for (int i = 1; i <= vertex_count; ++i) {
			int vertex_id, latitude, longitude;
			fin >> vertex_id >> latitude >> longitude;
			vertex_coordinates.push_back({ (long double)latitude/1000000, (long double)longitude/1000000 });
		}
		fin.close();
	}
	clock_t end_time = clock();
	return (end_time - start_time);
}

pair<int, int> findArcsWeight(vector<vector<pair<int, int>>>& source_graph,
							  int& from_vertex, int& to_vertex,
							  vector<tuple<int, int, int>>& temporary_arc_copies) {
	for (tuple<int, int, int> arc : temporary_arc_copies) {
		if (from_vertex == get<0>(arc) and to_vertex == get<1>(arc))
			return { get<2>(arc), 1 };
	}
	for (int i = 0; i < source_graph[from_vertex].size(); ++i) {
		if (to_vertex == source_graph[from_vertex][i].first) {
			int cost = source_graph[from_vertex][i].second;
			temporary_arc_copies.push_back({ from_vertex, to_vertex, cost });
			source_graph[from_vertex].erase(source_graph[from_vertex].begin() + i);
			return { cost, 0 };
		}
	}
}

double splitGraph(vector<vector<pair<int, int>>>& source_graph, 
				vector<vector<pair<int, int>>>& splitted_graph,
				vector<tuple<int, int, int>>& banned_turns) {
	clock_t start_time = clock();
	if (banned_turns.size() != 0) {
		vector<tuple<int, int, int>> temporary_arc_copies;
		int source_graph_size = source_graph.size();
		splitted_graph.resize(2 * source_graph_size - 1);
		for (int i = 0; i < banned_turns.size(); ++i) {
			int from_vertex = get<0>(banned_turns[i]);
			int by_vertex = get<1>(banned_turns[i]);
			int to_vertex = get<2>(banned_turns[i]);
			pair<int, int> cost_and_bool;
			cost_and_bool = findArcsWeight(source_graph, from_vertex, by_vertex, temporary_arc_copies);
			if (!cost_and_bool.second)
				splitted_graph[from_vertex].push_back({ by_vertex + source_graph_size - 1, cost_and_bool.first });
			cost_and_bool = findArcsWeight(source_graph, by_vertex, to_vertex, temporary_arc_copies);
			if (!cost_and_bool.second)
				splitted_graph[by_vertex].push_back({ to_vertex + source_graph_size - 1, cost_and_bool.first });
		}
		for (int i = 1; i < source_graph.size(); ++i) {
			for (int j = 0; j < source_graph[i].size(); ++j) {
				pair<int, int> arc = source_graph[i][j];
				splitted_graph[i].push_back({ arc.first, arc.second });
				splitted_graph[i + source_graph_size - 1].push_back({ arc.first, arc.second });
			}
		}
	}
	else {
		splitted_graph = source_graph;
	}
	clock_t end_time = clock();
	return (end_time - start_time);
}

int generateRandomNumber(int start_range, int end_range) {
	random_device                  rand_dev;
	mt19937                        generator(rand_dev());
	uniform_int_distribution<int>  distr(start_range, end_range);
	return distr(generator);
}

bool findInBannedTurnsVector(vector<tuple<int, int, int>>& banned_turns,
							 tuple<int, int, int> banned_turn) {
	for (tuple<int, int, int> temp : banned_turns) {
		if ((get<0>(temp) == get<0>(banned_turn))
			and (get<1>(temp) == get<1>(banned_turn))
			and (get<2>(temp) == get<2>(banned_turn))) {
				return true;
		}
	}
	return false;
}

double generateTurnRestrictions(vector<vector<pair<int, int>>>& graph,
							  vector<tuple<int, int, int>>& banned_turns,
							  int& vertex_count, int& arc_count) {
	clock_t start_time = clock();
	char source_file[] = "random_turn_r.txt";
	ofstream fin;
	fin.open(source_file, ios_base::app);
	int turn_restrictions_count = arc_count / 50;
	for (int i = 1; i <= turn_restrictions_count; ++i) {
		int from_vertex, by_vertex, to_vertex;
		from_vertex = generateRandomNumber(1, vertex_count);
		if (graph[from_vertex].size()) {
			int by_vertex_index = generateRandomNumber(0, graph[from_vertex].size()-1);
			by_vertex = graph[from_vertex][by_vertex_index].first;
			int to_vertex_index = generateRandomNumber(0, graph[by_vertex].size()-1);
			to_vertex = graph[by_vertex][to_vertex_index].first;
			tuple<int, int, int> banned_turn = make_tuple(from_vertex, by_vertex, to_vertex);
			if (!findInBannedTurnsVector(banned_turns, banned_turn)) {
				banned_turns.push_back(banned_turn);
			}
		}
	}
	fin << banned_turns.size() << endl;
	for (tuple<int, int, int> banned_turn : banned_turns) {
		fin << get<0>(banned_turn) << " " << get<1>(banned_turn) << " " << get<2>(banned_turn) << endl;
	}
	fin.close();
	clock_t end_time = clock();
	return (end_time - start_time);
}

void printPath(vector<int>& ancestors, 
			   int& source_vertex, 
			   int& target_vertex, 
			   vector<int>& distances,
			   int& vertex_count,
			   bool& graph_is_splitted,
			   vector<pair<long double, long double>>& vertex_coordinates) {
	vector<int> path;
	int v = target_vertex;
	if (ancestors[v] == -1) {
		if (ancestors[v + vertex_count] == -1) {
			printf("There is no path to target_vertex = %d from the source_vertex = %d!\n", target_vertex, source_vertex);
			return;
		}
	}
	if (distances.size() > vertex_count + 1) {
		if (graph_is_splitted and (distances[v] >= distances[v + vertex_count])) {
			v = v + vertex_count;
		}
	}
	
	cout << "Distance to target_vertex: " << distances[v] << endl;
	for (; v != source_vertex; v = ancestors[v]) {
		if (v > vertex_count) {
			path.push_back(v - vertex_count);
		}
		else {
			path.push_back(v);
		}
	}
	path.push_back(source_vertex);
	reverse(path.begin(), path.end());
	cout << "Path from " << source_vertex << " to " << target_vertex << ":\n";
	//ofstream fout;
	//char pathh[] = "path2.csv";
	//fout.open(pathh);
	//fout << setprecision(9) << endl;
	for (size_t i = 0; i < path.size() - 1; ++i) {
		cout << path[i] << '-';
		//fout << vertex_coordinates[i + 1].first << "," << vertex_coordinates[i + 1].second << endl;
		if (i != 0 and i % 14 == 0)
			cout << endl;
	}
	cout << path[path.size() - 1] << endl;
	//fout << vertex_coordinates[path.size()].first << "," << vertex_coordinates[path.size() - 1].second << endl;
	//fout.close();
}

void printPath2(vector<pair<int, int>>& ancestors,
	int& source_vertex,
	int& target_vertex,
	vector<int>& distances,
	int& vertex_count,
	bool& graph_is_splitted,
	vector<pair<long double, long double>>& vertex_coordinates) {
	vector<int> path;
	int v = target_vertex;
	if (distances[v] == INF) {
		printf("There is no path to target_vertex = %d from the source_vertex = %d!\n", target_vertex, source_vertex);
		return;
	}
	cout << "Distance to target_vertex: " << distances[v] << endl;
	path.push_back(v);
	for (; v != source_vertex;) {
		pair<int, int> temp = ancestors[v];
		path.push_back(temp.second);
		if (temp.first != -5) {
			path.push_back(temp.first);
			v = ancestors[v].first;
		}
		else {
			v = ancestors[v].second;
		}
		
	}
	//path.push_back(source_vertex);
	reverse(path.begin(), path.end());
	cout << "Path from " << source_vertex << " to " << target_vertex << ":\n";
	//ofstream fout;
	//char pathh[] = "path2.csv";
	//fout.open(pathh);
	//fout << setprecision(9) << endl;
	for (size_t i = 0; i < path.size() - 1; ++i) {
		cout << path[i] << '-';
		//fout << vertex_coordinates[i + 1].first << "," << vertex_coordinates[i + 1].second << endl;
		if (i != 0 and i % 14 == 0)
			cout << endl;
	}
	cout << path[path.size() - 1] << endl;
	//fout << vertex_coordinates[path.size()].first << "," << vertex_coordinates[path.size() - 1].second << endl;
	//fout.close();
}

double Dijkstra(vector<vector<pair<int, int>>>& adjacency_list,
				vector<int>& distances, 
				vector<int>& ancestors, 
				int& source_vertex, 
				int& target_vertex,
				int& vertex_count) {
	clock_t start_time = clock();
	clock_t end_time;
	priority_queue<pair<int, int>> Queue;
	distances.assign(adjacency_list.size(), INF);
	ancestors.assign(adjacency_list.size(), -1);
	distances[source_vertex] = 0;
	Queue.push({ distances[source_vertex], source_vertex });
	while (!Queue.empty()) {
		pair<int, int> u = Queue.top(); Queue.pop();
		if ((u.second == target_vertex) or (u.second == target_vertex + vertex_count)) {
			end_time = clock();
			return (end_time - start_time);
		}
		for (pair<int, int> arc : adjacency_list[u.second]) {
			int v = arc.first;
			int alt = (distances[u.second] + arc.second);
			if (distances[v] > alt) {
				distances[v] = alt;
				Queue.push({ -alt, v });
				ancestors[v] = u.second;
			}
		}
	}
	end_time = clock();
	return (end_time - start_time);
}

void getAllowedTargets(vector<vector<pair<int, int>>>& adjacency_list,
				vector<tuple<int, int, int>>& banned_turns,
				int& s, int& t,
				vector<pair<int, int>>& targets) {
	for (int i = 0; i < adjacency_list[t].size(); ++i) {
		tuple<int, int, int> turn = make_tuple(s, t, adjacency_list[t][i].first);
		if (!findInBannedTurnsVector(banned_turns, turn)) {
			targets.push_back({adjacency_list[t][i].first, adjacency_list[t][i].second });
		}
	}
}

double DijkstraModified(vector<vector<pair<int, int>>>& adjacency_list,
				vector<int>& distances,
				vector<pair<int, int>>& ancestors,
				int& source_vertex,
				int& target_vertex,
				int& vertex_count,
				vector<tuple<int, int, int>>& banned_turns) {
	clock_t start_time = clock();
	clock_t end_time;
	priority_queue<pair<int, int>> Queue;
	vector<pair<int, int>> targets;
	distances.assign(adjacency_list.size(), INF);
	ancestors.assign(adjacency_list.size(), {});
	distances[source_vertex] = 0;
	Queue.push({ distances[source_vertex], source_vertex });
	while (!Queue.empty()) {
		pair<int, int> u = Queue.top(); Queue.pop();
		if ((u.second == target_vertex)) {
			end_time = clock();
			return (end_time - start_time);
		}

		for (pair<int, int> arc : adjacency_list[u.second]) {
			int v = arc.first;
			targets.clear();		
			//distances[v] = distances[u.second] + arc.second;
			//ancestors[v] = { ancestors[u.second].second, u.second };
			/*if (v == target_vertex) {
				distances[v] = distances[u.second] + arc.second;
				ancestors[v] = { -5, u.second };
				end_time = clock();
				return (end_time - start_time);
			}*/
			//Queue.push({ -(distances[u.second] + arc.second), v });
			if (v == target_vertex) {
				if (ancestors[u.second].second != 0) {
					if (!findInBannedTurnsVector(banned_turns, { ancestors[u.second].second, u.second, v })) {
						if (distances[v] > distances[u.second] + arc.second) {
							distances[v] = distances[u.second] + arc.second;
							ancestors[v] = { -5, u.second };
							end_time = clock();
							return (end_time - start_time);
						}
						
					}
					else {
						continue;
					}
				}
				/*else {
					if (distances[v] > distances[u.second] + arc.second) {
						distances[v] = distances[u.second] + arc.second;
						ancestors[v] = { -5, u.second };
						end_time = clock();
						return (end_time - start_time);
					}
				}*/
			}
			getAllowedTargets(adjacency_list, banned_turns, u.second, v, targets);
			for (int i = 0; i < targets.size(); ++i) {
				int alt = distances[u.second] + arc.second + targets[i].second;
				if (distances[targets[i].first] > alt ) {
					distances[targets[i].first] = alt;
					Queue.push({ -alt, targets[i].first });
					ancestors[targets[i].first] = make_pair(u.second, v);
					//ancestors[targets[i].first] = v;
				}
				/*if ((targets[i].first == target_vertex) or (targets[i].first == target_vertex + vertex_count)) {
					end_time = clock();
					return (end_time - start_time);
				}*/
			}
		}
	}
	end_time = clock();
	return (end_time - start_time);
}

double BellmanFord(vector<vector<pair<int, int>>>& adjacency_list,
				   vector<int>& distances,
				   vector<int>& ancestors,
				   int& source_vertex,
				   int& target_vertex) {
	clock_t start_time = clock();
	distances.assign(adjacency_list.size(), INF);
	ancestors.assign(adjacency_list.size(), -1);
	distances[source_vertex] = 0;
	for (;;) {
		bool flag = false;
		for (int i = 1; i < adjacency_list.size(); ++i) {
			for (pair<int, int> arc : adjacency_list[i]) {
				int from = i;
				int to = arc.first;
				int cost = arc.second;
				if (distances[to] > distances[from] + cost) {
					distances[to] = distances[from] + cost;
					ancestors[to] = from;
					flag = true;
				}
			}
		}
		if (not flag)
			break;
	}
	clock_t end_time = clock();
	return (double)(end_time - start_time);
}

long double toRadians(const long double& degree) {
	return ((long double)(M_PI/180) * degree);
}

int heuristic(vector<pair<long double, long double>>& vertex_coordinates, int vertex_1, int vertex_2) {
	if (vertex_1 >= vertex_coordinates.size()) {
		vertex_1 -= vertex_coordinates.size();
	}
	if (vertex_2 >= vertex_coordinates.size()) {
		vertex_2 -= vertex_coordinates.size();
	}
	long double latitude_1; 
	long double longitude_1; 
	long double latitude_2; 
	long double longitude_2;

	latitude_1 = toRadians(vertex_coordinates[vertex_1].first);
	longitude_1 = toRadians(vertex_coordinates[vertex_1].second);
	latitude_2 = toRadians(vertex_coordinates[vertex_2].first);
	longitude_2 = toRadians(vertex_coordinates[vertex_2].second);

	long double d_longitude = longitude_2 - longitude_1;
	long double d_latitude = latitude_2 - latitude_1;

	long double ans = pow(sin(d_latitude / 2), 2) + 
		cos(latitude_1) * cos(latitude_2) * pow(sin(d_longitude / 2), 2);

	ans = 2 * asin(sqrt(ans));
	long double R = 6371;
	ans = 1000 * ans * R;
	return (int)ans;
}

double Astar(vector<vector<pair<int, int>>>& adjacency_list,
			 vector<int>& distances,
			 vector<int>& ancestors,
			 int& source_vertex,
			 int& target_vertex,
			 int& vertex_count,
			 vector<pair<long double, long double>>& vertex_coordinates) {
	clock_t start_time = clock();
	clock_t end_time;
	priority_queue<pair<int, int>> Queue;
	distances.assign(adjacency_list.size(), INF);
	ancestors.assign(adjacency_list.size(), -1);
	distances[source_vertex] = 0;
	Queue.push({ distances[source_vertex], source_vertex });
	while (!Queue.empty()) {
		pair<int, int> u = Queue.top(); Queue.pop();
		if ((u.second == target_vertex) or (u.second == target_vertex + vertex_count)) {
			end_time = clock();
			return (end_time - start_time);
		}
		for (pair<int, int> arc : adjacency_list[u.second]) {
			int v = arc.first;
			int alt = (distances[u.second] + arc.second);
			if (distances[v] > alt) {
				distances[v] = alt;
				Queue.push({ -(alt + heuristic(vertex_coordinates, v, target_vertex)), v });
				ancestors[v] = u.second;
			}
		}
	}
	end_time = clock();
	return (end_time - start_time);
}

double AstarModified(vector<vector<pair<int, int>>>& adjacency_list,
					vector<int>& distances,
					vector<pair<int, int>>& ancestors,
					int& source_vertex,
					int& target_vertex,
					int& vertex_count,
					vector<pair<long double, long double>>& vertex_coordinates,
					vector<tuple<int, int, int>>& banned_turns) {
	clock_t start_time = clock();
	clock_t end_time;
	priority_queue<pair<int, int>> Queue;
	vector<pair<int, int>> targets;
	distances.assign(adjacency_list.size(), INF);
	ancestors.assign(adjacency_list.size(), {});
	distances[source_vertex] = 0;
	Queue.push({ distances[source_vertex], source_vertex });
	while (!Queue.empty()) {
		pair<int, int> u = Queue.top(); Queue.pop();
		if ((u.second == target_vertex) or (u.second == target_vertex + vertex_count)) {
			end_time = clock();
			return (end_time - start_time);
		}
		for (pair<int, int> arc : adjacency_list[u.second]) {
			int v = arc.first;
			targets.clear();
			/*if (v == target_vertex) {
				distances[v] = distances[u.second] + arc.second;
				ancestors[v] = { -5, u.second };
				end_time = clock();
				return (end_time - start_time);
			}*/
			if (v == target_vertex) {
				if (ancestors[u.second].second != 0) {
					if (!findInBannedTurnsVector(banned_turns, { ancestors[u.second].second, u.second, v })) {
						if (distances[v] > distances[u.second] + arc.second) {
							distances[v] = distances[u.second] + arc.second;
							ancestors[v] = { -5, u.second };
							end_time = clock();
							return (end_time - start_time);
						}

					}
					else {
						continue;
					}
				}
				/*else {
					if (distances[v] > distances[u.second] + arc.second) {
						distances[v] = distances[u.second] + arc.second;
						ancestors[v] = { -5, u.second };
						end_time = clock();
						return (end_time - start_time);
					}
				}*/
			}
			getAllowedTargets(adjacency_list, banned_turns, u.second, v, targets);
			for (int i = 0; i < targets.size(); ++i) {
				int alt = distances[u.second] + arc.second + targets[i].second;
				if (distances[targets[i].first] > alt) {
					distances[targets[i].first] = alt;
					Queue.push({ -(alt + heuristic(vertex_coordinates, targets[i].first, target_vertex)), targets[i].first });
					ancestors[targets[i].first] = make_pair(u.second, v);
				}
				/*if ((targets[i].first == target_vertex) or (targets[i].first == target_vertex + vertex_count)) {
					end_time = clock();
					return (end_time - start_time);
				}*/
			}
		}
	}
	end_time = clock();
	return (end_time - start_time);
}

int main() {
	vector<vector<pair<int, int>>> graph;
	vector<tuple<int, int, int>> banned_turns;
	vector<pair<long double, long double>> vertex_coordinates;
	bool data_with_coordinates = true;	// <-
	bool data_without_coordinates = false;	// <-
	bool graph_is_splitted = true;		// <-
	bool graph_is_not_splitted = false;
	int vertex_count, arc_count;
	vector<int> distances;
	vector<int> ancestors1;
	vector<pair<int, int>> ancestors2;

	int source_vertex = 1; int target_vertex = 4;
	double runtime = readInputData(graph, banned_turns, vertex_count, arc_count, vertex_coordinates, data_with_coordinates); // <-
	cout << "Runtime(Reading and initializing graph): " << runtime / CLOCKS_PER_SEC << "s\n\n";
	
	/****************** Without turn restrictions ******************/
	cout << "Heuristic between source and target: " << heuristic(vertex_coordinates, source_vertex, target_vertex) << endl; // <-
	runtime = Dijkstra(graph, distances, ancestors1, source_vertex, target_vertex, vertex_count);
	cout << "\n*** Graph without restrictions ***" << endl;
	cout << "Runtime(Dijkstra's algorithm): " << runtime / CLOCKS_PER_SEC << "s \n";
	printPath(ancestors1, source_vertex, target_vertex, distances, vertex_count, graph_is_not_splitted, vertex_coordinates);
	
	runtime = Astar(graph, distances, ancestors1, source_vertex, target_vertex, vertex_count, vertex_coordinates);
	cout << "\n*** Graph without restrictions ***" << endl;
	cout << "Runtime(A* algorithm): " << runtime / CLOCKS_PER_SEC << "s \n";
	printPath(ancestors1, source_vertex, target_vertex, distances, vertex_count, graph_is_not_splitted, vertex_coordinates);

	runtime = BellmanFord(graph, distances, ancestors1, source_vertex, target_vertex);
	cout << "\n*** Graph without restrictions ***" << endl;
	cout << "Runtime(Bellman-Ford algorithm): " << runtime / CLOCKS_PER_SEC << "s \n";
	printPath(ancestors1, source_vertex, target_vertex, distances, vertex_count, graph_is_not_splitted, vertex_coordinates);

	/********************* With turn restrictions ******************/
	vector<vector<pair<int, int>>> splitted_graph, copy_of_graph;
	copy_of_graph = graph;
	runtime = generateTurnRestrictions(graph, banned_turns, vertex_count, arc_count);
	cout << "\nRuntime(Generating turn restriction): " << runtime / CLOCKS_PER_SEC << "s \n";
	cout << "Count of generated turns: " << arc_count / 70 << endl;
	runtime = DijkstraModified(graph, distances, ancestors2, source_vertex, target_vertex, vertex_count, banned_turns);
	cout << "\n*** Graph with turn restrictions ***" << endl;
	cout << "Runtime(Modified Dijkstra's algorithm): " << runtime / CLOCKS_PER_SEC << "s \n";
	printPath2(ancestors2, source_vertex, target_vertex, distances, vertex_count, graph_is_not_splitted, vertex_coordinates);	

	runtime = AstarModified(graph, distances, ancestors2, source_vertex, target_vertex, vertex_count, vertex_coordinates, banned_turns);
	cout << "\n*** Graph with turn restrictions ***" << endl;
	cout << "Runtime(Modified A* algorithm): " << runtime / CLOCKS_PER_SEC << "s \n";
	printPath2(ancestors2, source_vertex, target_vertex, distances, vertex_count, graph_is_not_splitted, vertex_coordinates);

	/******************* With forbidden sequences ******************/
	runtime = splitGraph(copy_of_graph, splitted_graph, banned_turns);
	cout << "\nRuntime(Splitting graph): " << runtime / CLOCKS_PER_SEC << "s \n";

	runtime = Dijkstra(splitted_graph, distances, ancestors1, source_vertex, target_vertex, vertex_count);
	cout << "\n*** Graph with restrictions (forbidden subsequences) ***" << endl;
	cout << "Runtime(Dijkstra's algorithm for graph with forbidden sequences - splitted): " << runtime / CLOCKS_PER_SEC << "s \n";
	printPath(ancestors1, source_vertex, target_vertex, distances, vertex_count, graph_is_splitted, vertex_coordinates);

	runtime = BellmanFord(splitted_graph, distances, ancestors1, source_vertex, target_vertex);
	cout << "\n*** Graph with restrictions (forbidden subsequences) ***" << endl;
	cout << "Runtime(Bellman-Ford algorithm for graph with forbidden sequences - splitted): " << runtime / CLOCKS_PER_SEC << "s \n";
	printPath(ancestors1, source_vertex, target_vertex, distances, vertex_count, graph_is_splitted, vertex_coordinates);
	
	runtime = Astar(splitted_graph, distances, ancestors1, source_vertex, target_vertex, vertex_count, vertex_coordinates);
	cout << "\n*** Graph with restrictions (forbidden subsequences) ***" << endl;
	cout << "Runtime(A* algorithm for graph with forbidden sequences - splitted): " << runtime / CLOCKS_PER_SEC << "s \n";
	printPath(ancestors1, source_vertex, target_vertex, distances, vertex_count, graph_is_splitted, vertex_coordinates);

	cout << endl;
	system("pause");
	return 0;
}