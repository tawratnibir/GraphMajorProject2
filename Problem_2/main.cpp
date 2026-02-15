#include <bits/stdc++.h>
using namespace std;

#define M_PI 3.14159265358979323846

// Cost parameters
const double CAR_FARE = 20.0;    // Taka per km
const double METRO_FARE = 5.0;   // Taka per km

struct Point {
    double lat;
    double lon;

    Point() {}
    Point(double pLat, double pLan) {
        lat = pLat;
        lon = pLan;
    }
};

struct EdgeInfo {
    string type;         // "car" or "metro"
    double distance;     // km
    double cost;         // Taka
    string startStation; // for metro only
    string endStation;   // for metro only
    
    EdgeInfo() {}
    EdgeInfo(string t, double d, double c, string s1 = "", string s2 = "") {
        type = t;
        distance = d;
        cost = c;
        startStation = s1;
        endStation = s2;
    }
};

unordered_map<string, int> idMap;
int nextId = 0;
vector<Point> nodes;
map<pair<int,int>, EdgeInfo> edgeInfo;  // Store edge metadata
string key(Point p) {
    char buffer[100];
    sprintf(buffer, "%.6f,%.6f", p.lat, p.lon);
    return string(buffer);
}
int getNodeId(Point p) {
    string k = key(p);
    if(idMap.find(k) == idMap.end()) {
        idMap[k] = nextId;
        nodes.push_back(p);
        nextId++;
    }
    return idMap[k];
}

// Dijkstra's algorithm for shortest path
// Parameters:
//   adj - adjacency list with double weights
//   source - source node ID
//   nodeCount - total number of nodes
// Returns: pair of (distance map, parent map)
pair<map<int, double>, map<int, int>> dijkstra(
    map<int, vector<pair<int, double>>>& adj, 
    int source, 
    int nodeCount
) {
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    
    map<int, double> dist;
    map<int, int> parent;
    
    // Initialize distances to infinity and parent to -1
    for(int i = 0; i < nodeCount; i++) {
        dist[i] = 1e18;
        parent[i] = -1;
    }

    dist[source] = 0.0;
    pq.push({0.0, source});

    while (!pq.empty()) {
        double d = pq.top().first;
        int node = pq.top().second;
        pq.pop();

        // Skip if we've found a better path already
        if(d > dist[node]) continue;

        for(auto& edge : adj[node]) {
            int neighbor = edge.first;
            double weight = edge.second;

            if(dist[node] + weight < dist[neighbor]) {
                dist[neighbor] = dist[node] + weight;
                parent[neighbor] = node;
                pq.push({dist[neighbor], neighbor});
            }
        }
    }
    
    return {dist, parent};
}

double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
    const double R = 6371.0; // Earth radius in km
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    
    double a = sin(dLat/2) * sin(dLat/2) + 
               cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) * 
               sin(dLon/2) * sin(dLon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    
    return R * c;
}

// Parse one row of CSV and extract coordinate pairs
// Returns vector of Points from the row
vector<Point> parseCSVRow(string line) {
    vector<Point> coords;
    stringstream ss(line);
    string token;
    vector<string> tokens;
    
    // Split by comma and space
    while(getline(ss, token, ',')) {
        // Remove leading/trailing spaces
        token.erase(0, token.find_first_not_of(" \t\r\n"));
        token.erase(token.find_last_not_of(" \t\r\n") + 1);
        if(!token.empty()) {
            tokens.push_back(token);
        }
    }
    
    // Skip first token (DhakaStreet label)
    // Parse coordinate pairs: lon1, lat1, lon2, lat2, ... (CSV stores lon,lat)
    // Skip last two tokens (altitude and distance)
    for(size_t i = 1; i < tokens.size() - 2; i += 2) {
        if(i+1 < tokens.size()) {
            try {
                double lon = stod(tokens[i]);
                double lat = stod(tokens[i+1]);
                coords.push_back(Point(lat, lon));  // Store as lat, lon internally
            } catch(...) {
                // Skip invalid coordinates
                continue;
            }
        }
    }
    
    return coords;
}

// Parse entire CSV file
// Parameters: filename - path to Roadmap-Dhaka.csv
// Returns: vector of road polylines, each polyline is a vector of Points
vector<vector<Point>> parseCSV(string filename) {
    vector<vector<Point>> roads;
    ifstream file(filename);
    
    if(!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return roads;
    }
    
    string line;
    while(getline(file, line)) {
        if(line.empty()) continue;
        
        vector<Point> coords = parseCSVRow(line);
        if(coords.size() >= 2) {
            roads.push_back(coords);
        }
    }
    
    file.close();
    cout << "Parsed " << roads.size() << " road segments from CSV" << endl;
    return roads;
}

// Parse metro CSV row and extract coordinates and station names
// Returns: tuple of (coordinates, startStation, endStation)
tuple<vector<Point>, string, string> parseMetroRow(string line) {
    vector<Point> coords;
    string startStation = "";
    string endStation = "";
    
    stringstream ss(line);
    string token;
    vector<string> tokens;
    
    // Split by comma
    while(getline(ss, token, ',')) {
        token.erase(0, token.find_first_not_of(" \t\r\n"));
        token.erase(token.find_last_not_of(" \t\r\n") + 1);
        if(!token.empty()) {
            tokens.push_back(token);
        }
    }
    
    // Skip first token (DhakaMetroRail label)
    // Parse coordinate pairs: lon1, lat1, lon2, lat2, ...
    // Last two tokens are station names
    size_t i = 1;
    while(i < tokens.size() - 2) {
        // Check if this is a number (coordinate) or text (station name)
        try {
            double lon = stod(tokens[i]);
            double lat = stod(tokens[i+1]);
            coords.push_back(Point(lat, lon));
            i += 2;
        } catch(...) {
            // Reached station names
            break;
        }
    }
    
    // Get station names (last two tokens)
    if(tokens.size() >= 2) {
        startStation = tokens[tokens.size() - 2];
        endStation = tokens[tokens.size() - 1];
    }
    
    return make_tuple(coords, startStation, endStation);
}

// Parse metro CSV file
vector<tuple<vector<Point>, string, string>> parseMetroCSV(string filename) {
    vector<tuple<vector<Point>, string, string>> metroLines;
    ifstream file(filename);
    
    if(!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return metroLines;
    }
    
    string line;
    while(getline(file, line)) {
        if(line.empty()) continue;
        
        auto result = parseMetroRow(line);
        vector<Point> coords = get<0>(result);
        
        if(coords.size() >= 2) {
            metroLines.push_back(result);
        }
    }
    
    file.close();
    cout << "Parsed " << metroLines.size() << " metro segments from CSV" << endl;
    return metroLines;
}

// Build weighted graph from parsed road data
// Parameters: 
//   roads - vector of road polylines from parseCSV
//   adj - adjacency list to populate (pass by reference)
// Returns: void (modifies adj in place)
void buildGraph(vector<vector<Point>>& roads, map<int, vector<pair<int, double>>>& adj) {
    int edgeCount = 0;
    
    for(auto& road : roads) {
        // For each consecutive pair of points in the road
        for(size_t i = 0; i < road.size() - 1; i++) {
            Point p1 = road[i];
            Point p2 = road[i + 1];
            
            // Get or create node IDs
            int u = getNodeId(p1);
            int v = getNodeId(p2);
            
            // Calculate edge distance
            double distance = haversineDistance(p1.lat, p1.lon, p2.lat, p2.lon);
            double cost = CAR_FARE * distance;
            
            // Add bidirectional edge
            adj[u].push_back({v, cost});
            adj[v].push_back({u, cost});
            
            // Store edge info
            edgeInfo[{u, v}] = EdgeInfo("car", distance, cost);
            edgeInfo[{v, u}] = EdgeInfo("car", distance, cost);
            
            edgeCount++;
        }
    }
    
    cout << "Built road graph with " << nodes.size() << " nodes and " 
         << edgeCount << " edges" << endl;
}

// Add metro lines to the graph
void addMetroToGraph(vector<tuple<vector<Point>, string, string>>& metroLines, 
                     map<int, vector<pair<int, double>>>& adj) {
    int edgeCount = 0;
    
    for(auto& metroLine : metroLines) {
        vector<Point> coords = get<0>(metroLine);
        string startStation = get<1>(metroLine);
        string endStation = get<2>(metroLine);
        
        // Get start and end nodes
        Point p1 = coords[0];
        Point p2 = coords[coords.size() - 1];
        
        int u = getNodeId(p1);
        int v = getNodeId(p2);
        
        // Calculate total distance through all waypoints
        double distance = 0.0;
        for(size_t i = 0; i < coords.size() - 1; i++) {
            distance += haversineDistance(coords[i].lat, coords[i].lon, 
                                         coords[i+1].lat, coords[i+1].lon);
        }
        
        double cost = METRO_FARE * distance;
        
        // Add bidirectional edge
        adj[u].push_back({v, cost});
        adj[v].push_back({u, cost});
        
        // Store edge info with station names
        edgeInfo[{u, v}] = EdgeInfo("metro", distance, cost, startStation, endStation);
        edgeInfo[{v, u}] = EdgeInfo("metro", distance, cost, endStation, startStation);
        
        edgeCount++;
    }
    
    cout << "Added " << edgeCount << " metro edges to graph" << endl;
}

// Find nearest node to a given GPS coordinate
// Parameters: lat, lon - query point coordinates
// Returns: pair<int, double> - (nearest node ID, distance to it in km)
pair<int, double> findNearestNode(double lat, double lon) {
    int nearestNode = -1;
    double minDist = 1e18;
    
    for(int i = 0; i < nodes.size(); i++) {
        double dist = haversineDistance(lat, lon, nodes[i].lat, nodes[i].lon);
        if(dist < minDist) {
            minDist = dist;
            nearestNode = i;
        }
    }
    
    return {nearestNode, minDist};
}

// Reconstruct path from parent pointers
// Parameters:
//   parent - parent map from Dijkstra
//   source - source node ID
//   dest - destination node ID
// Returns: vector<int> - path from source to dest as node IDs
vector<int> reconstructPath(map<int, int>& parent, int source, int dest) {
    vector<int> path;
    
    if(parent[dest] == -1 && dest != source) {
        // No path exists
        return path;
    }
    
    int current = dest;
    while(current != -1) {
        path.push_back(current);
        current = parent[current];
    }
    
    reverse(path.begin(), path.end());
    return path;
}

// Generate text directions from path
// Parameters:
//   path - vector of node IDs
//   srcLat, srcLon - original source coordinates
//   destLat, destLon - original destination coordinates
//   srcSnapped, destSnapped - whether source/dest are on graph nodes
//   totalCost - total cost in Taka
//   totalDist - total distance in km
//   filename - output text file name
void generateTextDirections(vector<int>& path, double srcLat, double srcLon, 
                          double destLat, double destLon, bool srcSnapped, bool destSnapped,
                          double totalCost, double totalDist, string filename) {
    ofstream file(filename);
    
    if(!file.is_open()) {
        cerr << "Error: Cannot create file " << filename << endl;
        return;
    }
    
    file << fixed << setprecision(6);
    
    file << "Problem no : 2" << endl;
    file << "Minimum cost is " << setprecision(16) << totalCost << setprecision(6) << "৳ " << endl;
    file << "Source:  (" << srcLon << ", " << srcLat << ")" << endl;
    file << "Destination:  (" << destLon << ", " << destLat << ")" << endl;
    
    // Handle walking from source to first graph node
    if(!srcSnapped && path.size() > 0) {
        Point firstNode = nodes[path[0]];
        file << "Walk from Source (" << srcLon << ", " << srcLat 
             << ") to (" << firstNode.lon << ", " << firstNode.lat << ")" << endl;
    }
    
    // Generate turn-by-turn directions with transport type
    for(int i = 0; i < path.size() - 1; i++) {
        int u = path[i];
        int v = path[i + 1];
        Point from = nodes[u];
        Point to = nodes[v];
        
        // Get edge info
        EdgeInfo info = edgeInfo[{u, v}];
        
        // Output cost first
        file << "Cost: " << setprecision(16) << info.cost << setprecision(6) << "৳ . ";
        
        if(info.type == "metro") {
            // Metro segment
            if(i == 0 && srcSnapped) {
                file << " Ride Metro from Source " << info.startStation 
                     << "(" << from.lon << ", " << from.lat << ")";
            } else {
                file << " Ride Metro from " << info.startStation 
                     << "(" << from.lon << ", " << from.lat << ")";
            }
            
            if(i == path.size() - 2 && destSnapped) {
                file << " to Destination " << info.endStation 
                     << " (" << to.lon << ", " << to.lat << ")" << endl;
            } else {
                file << " to " << info.endStation 
                     << " (" << to.lon << ", " << to.lat << ")" << endl;
            }
        } else {
            // Car segment
            if(i == 0 && srcSnapped) {
                file << "Ride Car from Source (" << from.lon << ", " << from.lat << ")";
            } else {
                file << "Ride Car from  (" << from.lon << ", " << from.lat << ")";
            }
            
            if(i == path.size() - 2 && destSnapped) {
                file << " to Destination (" << to.lon << ", " << to.lat << ")" << endl;
            } else {
                file << " to (" << to.lon << ", " << to.lat << ")" << endl;
            }
        }
    }
    
    // Handle walking from last graph node to destination
    if(!destSnapped && path.size() > 0) {
        Point lastNode = nodes[path[path.size() - 1]];
        file << "Walk from  (" << lastNode.lon << ", " << lastNode.lat 
             << ") to Destination (" << destLon << ", " << destLat << ")" << endl;
    }
    
    file.close();
    cout << "Text directions saved to " << filename << endl;
}

// Generate KML file for map visualization
// Parameters:
//   path - vector of node IDs
//   filename - output KML file name
void generateKML(vector<int>& path, string filename) {
    ofstream file(filename);
    
    if(!file.is_open()) {
        cerr << "Error: Cannot create file " << filename << endl;
        return;
    }
    
    // KML header
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    file << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">" << endl;
    file << "  <Document>" << endl;
    file << "    <name>Route</name>" << endl;
    
    // Style for car route (red)
    file << "    <Style id=\"carLine\">" << endl;
    file << "      <LineStyle>" << endl;
    file << "        <color>ff0000ff</color>" << endl;
    file << "        <width>3</width>" << endl;
    file << "      </LineStyle>" << endl;
    file << "    </Style>" << endl;
    
    // Style for metro route (green)
    file << "    <Style id=\"metroLine\">" << endl;
    file << "      <LineStyle>" << endl;
    file << "        <color>ff00ff00</color>" << endl;
    file << "        <width>4</width>" << endl;
    file << "      </LineStyle>" << endl;
    file << "    </Style>" << endl;
    
    // Style for source marker (black)
    file << "    <Style id=\"blackMarker\">" << endl;
    file << "      <IconStyle>" << endl;
    file << "        <color>ff000000</color>" << endl;
    file << "        <scale>1.2</scale>" << endl;
    file << "      </IconStyle>" << endl;
    file << "    </Style>" << endl;
    
    // Style for destination marker (blue)
    file << "    <Style id=\"blueMarker\">" << endl;
    file << "      <IconStyle>" << endl;
    file << "        <color>ffff0000</color>" << endl;
    file << "        <scale>1.2</scale>" << endl;
    file << "      </IconStyle>" << endl;
    file << "    </Style>" << endl;
    
    file << fixed << setprecision(6);
    
    // Source marker
    Point source = nodes[path[0]];
    file << "    <Placemark>" << endl;
    file << "      <name>Source</name>" << endl;
    file << "      <styleUrl>#blackMarker</styleUrl>" << endl;
    file << "      <Point>" << endl;
    file << "        <coordinates>" << source.lon << "," << source.lat << ",0</coordinates>" << endl;
    file << "      </Point>" << endl;
    file << "    </Placemark>" << endl;
    
    // Destination marker
    Point destination = nodes[path[path.size() - 1]];
    file << "    <Placemark>" << endl;
    file << "      <name>Destination</name>" << endl;
    file << "      <styleUrl>#blueMarker</styleUrl>" << endl;
    file << "      <Point>" << endl;
    file << "        <coordinates>" << destination.lon << "," << destination.lat << ",0</coordinates>" << endl;
    file << "      </Point>" << endl;
    file << "    </Placemark>" << endl;
    
    // Route segments (different colors for car and metro)
    for(int i = 0; i < path.size() - 1; i++) {
        int u = path[i];
        int v = path[i + 1];
        Point from = nodes[u];
        Point to = nodes[v];
        
        EdgeInfo info = edgeInfo[{u, v}];
        string styleUrl = (info.type == "metro") ? "#metroLine" : "#carLine";
        string name = (info.type == "metro") ? 
                      "Metro: " + info.startStation + " to " + info.endStation :
                      "Car";
        
        file << "    <Placemark>" << endl;
        file << "      <name>" << name << "</name>" << endl;
        file << "      <styleUrl>" << styleUrl << "</styleUrl>" << endl;
        file << "      <LineString>" << endl;
        file << "        <coordinates>" << endl;
        file << "          " << from.lon << "," << from.lat << ",0" << endl;
        file << "          " << to.lon << "," << to.lat << ",0" << endl;
        file << "        </coordinates>" << endl;
        file << "      </LineString>" << endl;
        file << "    </Placemark>" << endl;
    }
    
    file << "  </Document>" << endl;
    file << "</kml>" << endl;
    
    file.close();
    cout << "KML file saved to " << filename << endl;
}

int main() {
    cout << fixed << setprecision(7);
    
    // Parse CSV and build graph
    cout << "=== Building Graph for Problem 2 ===" << endl;
    string csvFile = "../Dataset/Roadmap-Dhaka.csv";
    string metroFile = "../Dataset/Routemap-DhakaMetroRail.csv";
    
    vector<vector<Point>> roads = parseCSV(csvFile);
    map<int, vector<pair<int, double>>> adj;
    buildGraph(roads, adj);
    
    // Add metro lines to graph
    vector<tuple<vector<Point>, string, string>> metroLines = parseMetroCSV(metroFile);
    addMetroToGraph(metroLines, adj);
    
    cout << "\nTotal nodes in combined graph: " << nodes.size() << endl;
    
    // Test with manual input
    cout << "\nEnter source coordinates (lat lon): ";
    double srcLat, srcLon;
    cin >> srcLat >> srcLon;
    
    cout << "Enter destination coordinates (lat lon): ";
    double destLat, destLon;
    cin >> destLat >> destLon;
    
    // Find nearest nodes for source and destination
    cout << "\nFinding nearest nodes..." << endl;
    pair<int, double> srcResult = findNearestNode(srcLat, srcLon);
    pair<int, double> destResult = findNearestNode(destLat, destLon);
    
    int srcNode = srcResult.first;
    double srcDist = srcResult.second;
    int destNode = destResult.first;
    double destDist = destResult.second;
    
    // Check if source/destination are exactly on graph nodes (threshold: 0.001 km)
    bool srcSnapped = (srcDist < 0.001);
    bool destSnapped = (destDist < 0.001);
    
    cout << "Source snapped to node " << srcNode 
         << " at (" << nodes[srcNode].lat << ", " << nodes[srcNode].lon << ")"
         << " [distance: " << srcDist << " km]" << endl;
    cout << "Destination snapped to node " << destNode 
         << " at (" << nodes[destNode].lat << ", " << nodes[destNode].lon << ")"
         << " [distance: " << destDist << " km]" << endl;
    
    // Run Dijkstra (using cost as weight)
    cout << "\nRunning Dijkstra's algorithm (minimizing cost)..." << endl;
    pair<map<int, double>, map<int, int>> result = dijkstra(adj, srcNode, nodes.size());
    map<int, double> cost = result.first;
    map<int, int> parent = result.second;
    
    // Check if destination is reachable
    if(cost[destNode] >= 1e17) {
        cout << "No path found from source to destination!" << endl;
        return 1;
    }
    
    // Reconstruct path
    vector<int> path = reconstructPath(parent, srcNode, destNode);
    
    // Calculate total distance
    double totalDist = 0.0;
    for(int i = 0; i < path.size() - 1; i++) {
        int u = path[i];
        int v = path[i + 1];
        EdgeInfo info = edgeInfo[{u, v}];
        totalDist += info.distance;
    }
    
    cout << "\n=== ROUTE FOUND ===" << endl;
    cout << "Total cost: " << fixed << setprecision(2) 
         << cost[destNode] << " Taka" << endl;
    cout << "Total distance: " << setprecision(3) 
         << totalDist << " km" << endl;
    cout << "Number of nodes in path: " << path.size() << endl;
    
    cout << "\nPath (first 10 segments):" << endl;
    for(int i = 0; i < min(10, (int)path.size() - 1); i++) {
        int u = path[i];
        int v = path[i + 1];
        EdgeInfo info = edgeInfo[{u, v}];
        cout << "  " << info.type << ": ";
        if(info.type == "metro") {
            cout << info.startStation << " -> " << info.endStation;
        } else {
            cout << "(" << nodes[u].lat << ", " << nodes[u].lon << ") -> "
                 << "(" << nodes[v].lat << ", " << nodes[v].lon << ")";
        }
        cout << endl;
    }
    if(path.size() > 11) {
        cout << "  ... (" << (path.size() - 11) << " more segments)" << endl;
    }
    
    // Generate outputs
    cout << "\nGenerating output files..." << endl;
    generateTextDirections(path, srcLat, srcLon, destLat, destLon, 
                          srcSnapped, destSnapped, cost[destNode], totalDist, "route_output.txt");
    generateKML(path, "route.kml");
    
    cout << "\nDone! You can now:" << endl;
    cout << "  - View route_output.txt for turn-by-turn directions" << endl;
    cout << "  - Open route.kml in Google Earth or Google MyMaps" << endl;
    cout << "  - Red lines = Car segments, Green lines = Metro segments" << endl;
    
    return 0;
}