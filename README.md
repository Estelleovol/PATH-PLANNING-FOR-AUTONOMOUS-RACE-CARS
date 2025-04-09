# PATH-PLANNING-FOR-AUTONOMOUS-RACE-CARS
# Path Planning for Autonomous Race Cars

**Author:** Estelle Liu

This repository implements a path-planning pipeline designed for autonomous race cars, optimizing racing trajectories based on distance, curvature, and execution time.

## Workflow Overview

1. **Load Track Data**
   - Import track information from a CSV file.

2. **Preprocessing**
   - Interpolate track points and compute edge boundaries.

3. **Build Quadratic Programming (QP) Matrices**
   - Generate matrices required for optimization.

4. **Optimization (QP Solver)**
   - Solve for the shortest path while minimizing curvature using Quadratic Programming.

5. **Grey Wolf Optimizer**
   - Further optimize trajectory based on distance and time constraints.

6. **Extract Optimal Path**
   - Determine optimal racing line coordinates (alpha → x,y).

7. **Evaluation**
   - Assess trajectory performance (distance, curvature, etc.).

8. **Visualization**
   - Generate insightful plots and heatmaps for path visualization.


## Visualization
Sample visualizations (plots and heatmaps) illustrating optimized paths are provided.


---

*Estelle Liu © 2024*

