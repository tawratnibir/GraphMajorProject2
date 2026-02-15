## Route Planning System - Problem 6: Minimize Time with Realistic Speeds

### Current Implementation:

**Problem 6**: Minimize Time with Realistic Speeds
- Optimizes for fastest route using transport-specific speeds
- Uses time-based Dijkstra with realistic speed parameters:
  * **Car**: 20 km/h
  * **Metro**: 15 km/h
  * **Uttara Bus**: 12 km/h
  * **Bikolpo Bus**: 10 km/h
- Public transport operates 5:30 AM - 11:59 PM
- 15-minute waiting intervals for metro/buses
- Walking speed: 2 km/h
- Input: Source, destination, start time
- Output: Total time in minutes + time at each step

### Comparison with Problem 5:

| Feature | Problem 5 | Problem 6 |
|---------|-----------|-----------|
| Optimization | Minimize Time | Minimize Time |
| Speed Model | Uniform (10 km/h) | Realistic per transport |
| Car Speed | 10 km/h | 20 km/h |
| Metro Speed | 10 km/h | 15 km/h |
| Bus Speeds | 10 km/h | 10-12 km/h |

Problem 6 produces more realistic routing as faster transport modes (car, metro) are properly weighted.

### How to Run:

```powershell
# Compile
g++ -std=c++17 -o dij.exe dij.cpp

# Run with test input
Get-Content test_input_p6.txt | .\dij.exe
```

### Input Format:

```
<source_lat> <source_lon>
<dest_lat> <dest_lon>
<hour> <minute>
```

Example:
```
23.8103 90.4125
23.7504 90.3884
9 30
```

### Output Files:
- `route_output.txt` - Turn-by-turn directions with realistic time tracking
- `route.kml` - Color-coded map visualization
  - Red = Car (20 km/h)
  - Green = Metro (15 km/h)
  - Blue = Uttara Bus (12 km/h)
  - Yellow = Bikolpo Bus (10 km/h)

### Test Coordinates:
- Test route: From (23.8103, 90.4125) to (23.7504, 90.3884)
- Start time: 9:30 AM

### Note:
To change which problem is solved, edit `problemNum` variable in `main()` function:
- `problemNum = 3` - Minimize cost (no time constraints)
- `problemNum = 4` - Minimize cost (with time constraints)
- `problemNum = 5` - Minimize time (uniform 10 km/h speed)
- `problemNum = 6` - Minimize time (realistic speeds) - **current default**

### All 6 Problems Summary:

1. **Problem 1**: Shortest distance (car only)
2. **Problem 2**: Minimum cost (car + metro)
3. **Problem 3**: Minimum cost (all 4 transport modes)
4. **Problem 4**: Minimum cost with time constraints (time-aware)
5. **Problem 5**: Minimum time (uniform speed assumption)
6. **Problem 6**: Minimum time (realistic speeds per transport) ✓
