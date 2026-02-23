#include <bits/stdc++.h>
using namespace std;

#define M_PI 3.14159265358979323846

// Cost parameters for Problem 3
const double CAR_FARE = 20.0;       // Taka per km
const double METRO_FARE = 5.0;      // Taka per km  
const double UTTARA_FARE = 10.0;    // Taka per km
const double BIKOLPO_FARE = 7.0;    // Taka per km

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
    string type;         // "car", "metro", "uttara_bus", or "bikolpo_bus"
    double distance;     // km
    double cost;         // Taka
    string startStation; // for public transport
    string endStation;   // for public transport
    
    EdgeInfo() {}
    EdgeInfo(string t, double d, double c, string s1 = "", string s2 = "") {
        type = t;
        distance = d;
        cost = c;
        startStation = s1;
        endStation = s2;
    }
};

// Time structure for Problem 4
struct TimeInfo {
    int hour;
    int minute;
    
    TimeInfo() : hour(0), minute(0) {}
    TimeInfo(int h, int m) : hour(h), minute(m) {}
};

// Add minutes to time (for Problem 4)
TimeInfo addMinutes(TimeInfo t, double minutes) {
    int totalMinutes = t.hour * 60 + t.minute + (int)minutes;
    int newHour = (totalMinutes / 60) % 24;
    int newMinute = totalMinutes % 60;
    return TimeInfo(newHour, newMinute);
}

// Get waiting time until next 15-minute interval (for Problem 4)
// Returns minutes to wait
int getWaitingTime(TimeInfo t) {
    if(t.minute <= 15) return 15 - t.minute;
    if(t.minute <= 30) return 30 - t.minute;
    if(t.minute <= 45) return 45 - t.minute;
    return 60 - t.minute;
}

// Check if public transport is operating (5:30 AM - 11:59 PM) for Problem 4
bool isPublicTransportOperating(TimeInfo t) {
    if(t.hour == 5) {
        return t.minute >= 30;
    }
    return t.hour >= 6 && t.hour <= 23;
}

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

// Dijkstra's algorithm with time constraints for Problem 5
// Parameters:
//   adj - adjacency list with double weights (time in hours)
//   source - source node ID
//   nodeCount - total number of nodes
//   startTime - journey start time
// Returns: tuple of (time map in hours, parent map, arrival time map)
tuple<map<int, double>, map<int, int>, map<int, TimeInfo>> dijkstraWithTimeP5(
    map<int, vector<pair<int, double>>>& adj, 
    int source, 
    int nodeCount,
    TimeInfo startTime
) {
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    
    map<int, double> totalTime;  // Total time in hours
    map<int, int> parent;
    map<int, TimeInfo> arrivalTime;
    
    // Initialize times to infinity and parent to -1
    for(int i = 0; i < nodeCount; i++) {
        totalTime[i] = 1e18;
        parent[i] = -1;
    }

    totalTime[source] = 0.0;
    arrivalTime[source] = startTime;
    pq.push({0.0, source});

    while (!pq.empty()) {
        double t = pq.top().first;
        int node = pq.top().second;
        pq.pop();

        // Skip if we've found a better path already
        if(t > totalTime[node]) continue;
        
        TimeInfo timeAtNode = arrivalTime[node];

        for(auto& edge : adj[node]) {
            int neighbor = edge.first;
            double edgeTravelTimeHours = edge.second;  // Travel time from edge weight
            
            // Get edge metadata
            EdgeInfo info = edgeInfo[{node, neighbor}];
            
            // Calculate time components for this edge
            double waitTimeHours = 0.0;
            TimeInfo timeAtNeighbor = timeAtNode;
            
            // If public transport, add waiting time and check operating hours
            if(info.type != "car") {
                int waitMinutes = getWaitingTime(timeAtNeighbor);
                timeAtNeighbor = addMinutes(timeAtNeighbor, waitMinutes);
                waitTimeHours = waitMinutes / 60.0;
                
                // Check if public transport is operating at this time
                if(!isPublicTransportOperating(timeAtNeighbor)) {
                    continue; // Skip this edge - public transport not operating
                }
            }
            
            // Add travel time
            int travelMinutes = (int)(edgeTravelTimeHours * 60);
            timeAtNeighbor = addMinutes(timeAtNeighbor, travelMinutes);
            
            // Total time for this edge = wait time + travel time
            double edgeTotalTime = waitTimeHours + edgeTravelTimeHours;

            // Check if this is a better path (by total time)
            double newTotalTime = totalTime[node] + edgeTotalTime;
            if(newTotalTime < totalTime[neighbor]) {
                totalTime[neighbor] = newTotalTime;
                parent[neighbor] = node;
                arrivalTime[neighbor] = timeAtNeighbor;
                pq.push({totalTime[neighbor], neighbor});
            }
        }
    }
    
    return {totalTime, parent, arrivalTime};
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

// Parse bus CSV files (same format as metro)
// Returns: vector of (coordinates, startStation, endStation)
vector<tuple<vector<Point>, string, string>> parseBusCSV(string filename) {
    vector<tuple<vector<Point>, string, string>> busLines;
    ifstream file(filename);
    
    if(!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return busLines;
    }
    
    string line;
    while(getline(file, line)) {
        if(line.empty()) continue;
        
        // Reuse metro parsing function (same format)
        auto result = parseMetroRow(line);
        vector<Point> coords = get<0>(result);
        
        if(coords.size() >= 2) {
            busLines.push_back(result);
        }
    }
    
    file.close();
    cout << "Parsed " << busLines.size() << " bus segments from CSV" << endl;
    return busLines;
}

// ============== TIME-OPTIMIZED GRAPH BUILDING (Problem 5) ==============
// Build time-optimized graph (minimize travel time)
// Weight = distance/10 hours (assumes ~10 km/h base speed)
void buildGraphTime(vector<vector<Point>>& roads, map<int, vector<pair<int, double>>>& adj) {
    int edgeCount = 0;
    
    for(auto& road : roads) {
        for(size_t i = 0; i < road.size() - 1; i++) {
            Point p1 = road[i];
            Point p2 = road[i + 1];
            
            int u = getNodeId(p1);
            int v = getNodeId(p2);
            
            double distance = haversineDistance(p1.lat, p1.lon, p2.lat, p2.lon);
            double timeWeight = distance / 10.0;  // hours (10 km/h speed)
            double cost = CAR_FARE * distance;    // still track cost for output
            
            adj[u].push_back({v, timeWeight});
            adj[v].push_back({u, timeWeight});
            
            edgeInfo[{u, v}] = EdgeInfo("car", distance, cost);
            edgeInfo[{v, u}] = EdgeInfo("car", distance, cost);
            
            edgeCount++;
        }
    }
    
    cout << "Built time-optimized road graph with " << nodes.size() << " nodes and " 
         << edgeCount << " edges" << endl;
}

void addMetroToGraphTime(vector<tuple<vector<Point>, string, string>>& metroLines, 
                         map<int, vector<pair<int, double>>>& adj) {
    int edgeCount = 0;
    
    for(auto& metroLine : metroLines) {
        vector<Point> coords = get<0>(metroLine);
        string startStation = get<1>(metroLine);
        string endStation = get<2>(metroLine);
        
        Point p1 = coords[0];
        Point p2 = coords[coords.size() - 1];
        
        int u = getNodeId(p1);
        int v = getNodeId(p2);
        
        double distance = 0.0;
        for(size_t i = 0; i < coords.size() - 1; i++) {
            distance += haversineDistance(coords[i].lat, coords[i].lon, 
                                         coords[i+1].lat, coords[i+1].lon);
        }
        
        double timeWeight = distance / 10.0;  // hours
        double cost = METRO_FARE * distance;
        
        adj[u].push_back({v, timeWeight});
        adj[v].push_back({u, timeWeight});
        
        edgeInfo[{u, v}] = EdgeInfo("metro", distance, cost, startStation, endStation);
        edgeInfo[{v, u}] = EdgeInfo("metro", distance, cost, endStation, startStation);
        
        edgeCount++;
    }
    
    cout << "Added " << edgeCount << " metro edges (time-optimized)" << endl;
}

void addUttaraBusToGraphTime(vector<tuple<vector<Point>, string, string>>& busLines, 
                             map<int, vector<pair<int, double>>>& adj) {
    int edgeCount = 0;
    
    for(auto& busLine : busLines) {
        vector<Point> coords = get<0>(busLine);
        string startStation = get<1>(busLine);
        string endStation = get<2>(busLine);
        
        Point p1 = coords[0];
        Point p2 = coords[coords.size() - 1];
        
        int u = getNodeId(p1);
        int v = getNodeId(p2);
        
        double distance = 0.0;
        for(size_t i = 0; i < coords.size() - 1; i++) {
            distance += haversineDistance(coords[i].lat, coords[i].lon, 
                                         coords[i+1].lat, coords[i+1].lon);
        }
        
        double timeWeight = distance / 10.0;  // hours
        double cost = UTTARA_FARE * distance;
        
        adj[u].push_back({v, timeWeight});
        adj[v].push_back({u, timeWeight});
        
        edgeInfo[{u, v}] = EdgeInfo("uttara_bus", distance, cost, startStation, endStation);
        edgeInfo[{v, u}] = EdgeInfo("uttara_bus", distance, cost, endStation, startStation);
        
        edgeCount++;
    }
    
    cout << "Added " << edgeCount << " Uttara bus edges (time-optimized)" << endl;
}

void addBikolpoBusToGraphTime(vector<tuple<vector<Point>, string, string>>& busLines, 
                              map<int, vector<pair<int, double>>>& adj) {
    int edgeCount = 0;
    
    for(auto& busLine : busLines) {
        vector<Point> coords = get<0>(busLine);
        string startStation = get<1>(busLine);
        string endStation = get<2>(busLine);
        
        Point p1 = coords[0];
        Point p2 = coords[coords.size() - 1];
        
        int u = getNodeId(p1);
        int v = getNodeId(p2);
        
        double distance = 0.0;
        for(size_t i = 0; i < coords.size() - 1; i++) {
            distance += haversineDistance(coords[i].lat, coords[i].lon, 
                                         coords[i+1].lat, coords[i+1].lon);
        }
        
        double timeWeight = distance / 10.0;  // hours
        double cost = BIKOLPO_FARE * distance;
        
        adj[u].push_back({v, timeWeight});
        adj[v].push_back({u, timeWeight});
        
        edgeInfo[{u, v}] = EdgeInfo("bikolpo_bus", distance, cost, startStation, endStation);
        edgeInfo[{v, u}] = EdgeInfo("bikolpo_bus", distance, cost, endStation, startStation);
        
        edgeCount++;
    }
    
    cout << "Added " << edgeCount << " Bikolpo bus edges (time-optimized)" << endl;
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

// ============== TEXT DIRECTIONS OUTPUT (Problem 5) ==============

// Generate text directions with time tracking (for Problem 5 - minimizing time)
// Parameters similar to Problem 4 but shows TIME as optimization metric
void generateTextDirectionsMinTime(vector<int>& path, double srcLat, double srcLon, 
                          double destLat, double destLon, bool srcSnapped, bool destSnapped,
                          double totalTime, double totalDist, TimeInfo startTime, string filename) {
    // Calculate end time first by simulating the journey
    TimeInfo endTime = startTime;
    
    // Calculate initial walking time if needed
    if(!srcSnapped) {
        Point firstNode = nodes[path[0]];
        double walkDist = haversineDistance(srcLat, srcLon, firstNode.lat, firstNode.lon);
        int walkTime = (int)((walkDist / 2.0) * 60);
        endTime = addMinutes(endTime, walkTime);
    }
    
    // Calculate time through all path segments
    for(int i = 0; i < path.size() - 1; i++) {
        int u = path[i];
        int v = path[i + 1];
        EdgeInfo info = edgeInfo[{u, v}];
        
        // Add waiting time for public transport
        if(info.type != "car") {
            int waitTime = getWaitingTime(endTime);
            endTime = addMinutes(endTime, waitTime);
        }
        
        // Add travel time (10 km/h uniform speed)
        int travelTime = (int)((info.distance / 10.0) * 60);
        endTime = addMinutes(endTime, travelTime);
    }
    
    // Calculate final walking time if needed
    if(!destSnapped) {
        Point lastNode = nodes[path[path.size() - 1]];
        double walkDist = haversineDistance(lastNode.lat, lastNode.lon, destLat, destLon);
        int walkTime = (int)((walkDist / 2.0) * 60);
        endTime = addMinutes(endTime, walkTime);
    }
    
    // Now write the file
    ofstream file(filename);
    
    if(!file.is_open()) {
        cerr << "Error: Cannot create file " << filename << endl;
        return;
    }
    
    file << fixed << setprecision(6);
    
    file << "Problem no : 5" << endl;
    file << "Source: (" << srcLat << ", " << srcLon << ")" << endl;
    file << "Destination: (" << destLat << ", " << destLon << ")" << endl;
    file << "Start Time: " << startTime.hour << ":" << setfill('0') << setw(2) << startTime.minute << endl;
    file << "End Time: " << endTime.hour << ":" << setfill('0') << setw(2) << endTime.minute << endl;
    file << "Total Time: " << setprecision(2) << (totalTime * 60.0) << " minutes" << endl;
    file << "Total Distance: " << setprecision(3) << totalDist << " km" << endl;
    file << endl;
    
    file << setprecision(6);
    
    TimeInfo currentTime = startTime;
    
    // Generate turn-by-turn directions with time tracking
    for(int i = 0; i < path.size() - 1; i++) {
        int u = path[i];
        int v = path[i + 1];
        Point from = nodes[u];
        Point to = nodes[v];
        
        EdgeInfo info = edgeInfo[{u, v}];
        
        // Handle walking from source to first graph node
        if(i == 0 && !srcSnapped) {
            double walkDist = haversineDistance(srcLat, srcLon, from.lat, from.lon);
            int walkTime = (int)((walkDist / 2.0) * 60);
            currentTime = addMinutes(currentTime, walkTime);
            file << "Time: " << currentTime.hour << ":" << setfill('0') << setw(2) << currentTime.minute 
                 << ". Cost: " << setprecision(2) << (walkDist / 2.0 * 60) << " minutes. Walk from Source (" 
                 << srcLat << ", " << srcLon << ") to (" << from.lat << ", " << from.lon << ")" << endl;
        }
        
        // Handle public transport waiting time
        if(info.type != "car") {
            int waitTime = getWaitingTime(currentTime);
            currentTime = addMinutes(currentTime, waitTime);
        }
        
        // Calculate travel time for this segment (distance / 10 hours = distance * 6 minutes)
        double segmentTimeMinutes = info.distance * 6.0;
        
        // Output the segment action
        file << "Time: " << currentTime.hour << ":" << setfill('0') << setw(2) << currentTime.minute 
             << ". Cost: " << setprecision(2) << segmentTimeMinutes << " minutes. ";
        
        if(info.type == "metro") {
            if(i == 0 && srcSnapped) {
                file << "Ride Metro from Source " << info.startStation;
            } else {
                file << "Ride Metro from " << info.startStation;
            }
            file << " (" << from.lat << ", " << from.lon << ")";
            
            if(i == path.size() - 2 && destSnapped) {
                file << " to Destination " << info.endStation;
            } else {
                file << " to " << info.endStation;
            }
            file << " (" << to.lat << ", " << to.lon << ")" << endl;
        } else if(info.type == "uttara_bus") {
            if(i == 0 && srcSnapped) {
                file << "Ride Uttara Bus from Source " << info.startStation;
            } else {
                file << "Ride Uttara Bus from " << info.startStation;
            }
            file << " (" << from.lat << ", " << from.lon << ")";
            
            if(i == path.size() - 2 && destSnapped) {
                file << " to Destination " << info.endStation;
            } else {
                file << " to " << info.endStation;
            }
            file << " (" << to.lat << ", " << to.lon << ")" << endl;
        } else if(info.type == "bikolpo_bus") {
            if(i == 0 && srcSnapped) {
                file << "Ride Bikolpo Bus from Source " << info.startStation;
            } else {
                file << "Ride Bikolpo Bus from " << info.startStation;
            }
            file << " (" << from.lat << ", " << from.lon << ")";
            
            if(i == path.size() - 2 && destSnapped) {
                file << " to Destination " << info.endStation;
            } else {
                file << " to " << info.endStation;
            }
            file << " (" << to.lat << ", " << to.lon << ")" << endl;
        } else {
            // Car segment
            if(i == 0 && srcSnapped) {
                file << "Ride Car from Source (" << from.lat << ", " << from.lon << ")";
            } else {
                file << "Ride Car from (" << from.lat << ", " << from.lon << ")";
            }
            
            if(i == path.size() - 2 && destSnapped) {
                file << " to Destination (" << to.lat << ", " << to.lon << ")" << endl;
            } else {
                file << " to (" << to.lat << ", " << to.lon << ")" << endl;
            }
        }
        
        // Add travel time
        int travelTime = (int)segmentTimeMinutes;
        currentTime = addMinutes(currentTime, travelTime);
        
        // Handle walking from last graph node to destination
        if(i == path.size() - 2 && !destSnapped) {
            double walkDist = haversineDistance(to.lat, to.lon, destLat, destLon);
            int walkTime = (int)((walkDist / 2.0) * 60);
            currentTime = addMinutes(currentTime, walkTime);
            file << "Time: " << currentTime.hour << ":" << setfill('0') << setw(2) << currentTime.minute 
                 << ". Cost: " << setprecision(2) << (walkDist / 2.0 * 60) << " minutes. Walk from (" 
                 << to.lat << ", " << to.lon << ") to Destination (" << destLat << ", " << destLon << ")" << endl;
        }
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
    
    // Style for Uttara bus route (blue)
    file << "    <Style id=\"uttaraBusLine\">" << endl;
    file << "      <LineStyle>" << endl;
    file << "        <color>ffff0000</color>" << endl;
    file << "        <width>3</width>" << endl;
    file << "      </LineStyle>" << endl;
    file << "    </Style>" << endl;
    
    // Style for Bikolpo bus route (yellow)
    file << "    <Style id=\"bikolpoBusLine\">" << endl;
    file << "      <LineStyle>" << endl;
    file << "        <color>ff00ffff</color>" << endl;
    file << "        <width>3</width>" << endl;
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
        
        // Determine style and name based on transport type
        string styleUrl, name;
        if(info.type == "metro") {
            styleUrl = "#metroLine";
            name = "Metro: " + info.startStation + " to " + info.endStation;
        } else if(info.type == "uttara_bus") {
            styleUrl = "#uttaraBusLine";
            name = "Uttara Bus: " + info.startStation + " to " + info.endStation;
        } else if(info.type == "bikolpo_bus") {
            styleUrl = "#bikolpoBusLine";
            name = "Bikolpo Bus: " + info.startStation + " to " + info.endStation;
        } else {
            styleUrl = "#carLine";
            name = "Car";
        }
        
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
    cout << "=== Building Graph for Problem 5 ===" << endl;
    string csvFile = "../Dataset/Roadmap-Dhaka.csv";
    string metroFile = "../Dataset/Routemap-DhakaMetroRail.csv";
    string uttaraBusFile = "../Dataset/Routemap-UttaraBus.csv";
    string bikolpoBusFile = "../Dataset/Routemap-BikolpoBus.csv";
    
    vector<vector<Point>> roads = parseCSV(csvFile);
    map<int, vector<pair<int, double>>> adj;
    
    // Build graph - Problem 5: Time-optimized graph (uniform 10 km/h speed)
    buildGraphTime(roads, adj);
    
    vector<tuple<vector<Point>, string, string>> metroLines = parseMetroCSV(metroFile);
    addMetroToGraphTime(metroLines, adj);
    
    vector<tuple<vector<Point>, string, string>> uttaraBusLines = parseBusCSV(uttaraBusFile);
    addUttaraBusToGraphTime(uttaraBusLines, adj);
    
    vector<tuple<vector<Point>, string, string>> bikolpoBusLines = parseBusCSV(bikolpoBusFile);
    addBikolpoBusToGraphTime(bikolpoBusLines, adj);
    
    cout << "\nTotal nodes in combined graph: " << nodes.size() << endl;
    
    // Test with manual input
    cout << "\nEnter source coordinates (lat lon): ";
    double srcLat, srcLon;
    cin >> srcLat >> srcLon;
    
    cout << "Enter destination coordinates (lat lon): ";
    double destLat, destLon;
    cin >> destLat >> destLon;
    
    cout << "Enter starting time (hour minute, 24-hour format): ";
    int startHour, startMinute;
    cin >> startHour >> startMinute;
    TimeInfo startTime(startHour, startMinute);
    
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
    
    // Run time-aware Dijkstra (minimizing time with constraints)
    cout << "\nRunning time-aware Dijkstra's algorithm (minimizing travel time)..." << endl;
    map<int, double> totalTime;
    map<int, int> parent;
    map<int, TimeInfo> arrivalTimeMap;
    tie(totalTime, parent, arrivalTimeMap) = dijkstraWithTimeP5(adj, srcNode, nodes.size(), startTime);
    
    // Check if destination is reachable
    if(totalTime[destNode] >= 1e17) {
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
    cout << "Total time: " << fixed << setprecision(2) 
         << (totalTime[destNode] * 60.0) << " minutes" << endl;
    cout << "Total distance: " << setprecision(3) 
         << totalDist << " km" << endl;
    cout << "Number of nodes in path: " << path.size() << endl;
    
    cout << "\nPath (first 10 segments):" << endl;
    for(int i = 0; i < min(10, (int)path.size() - 1); i++) {
        int u = path[i];
        int v = path[i + 1];
        EdgeInfo info = edgeInfo[{u, v}];
        cout << "  " << info.type << ": ";
        if(info.type != "car") {
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
    
    // Generate output files
    cout << "\nGenerating output files..." << endl;
    generateTextDirectionsMinTime(path, srcLat, srcLon, destLat, destLon, 
                          srcSnapped, destSnapped, totalTime[destNode], totalDist, 
                          startTime, "route_output.txt");
    generateKML(path, "route.kml");
    
    cout << "\nDone! You can now:" << endl;
    cout << "  - View route_output.txt for time-minimized directions" << endl;
    cout << "  - Open route.kml in Google Earth or Google MyMaps" << endl;
    cout << "  - Red=Car, Green=Metro, Blue=Uttara Bus, Yellow=Bikolpo Bus" << endl;
    
    return 0;
}