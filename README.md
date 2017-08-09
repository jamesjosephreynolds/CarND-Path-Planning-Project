# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

For the controller/trajectory portion of this project, I relied on the guidance from the project walkthrough.  Prior to the walkthrough I tried many different iterations using JMT, and was unable to create a smooth transition between old path and new path.  As in the walkthrough, JMT is not used, in favor of a trajectory spline.

## Behavior

I created two custom classes to determine the behavior of the ego vehicle: `Lane` and `Behavior`.  `Lane` describes what is going on in each of the lanes: feasible speed, following distance, and obstruction.  `Behavior` takes `Lane lane` as an input, and applies a cost function to determine the target lane.

The basic logic for the `Lane` class is, for each lane, iterate through all of the objects in the sensor fusion, determine which one is the closest to the ego vehicles, and set the lane parameters according to that object.

```C++
// reset lane status
this->status[i] = "open";
this->speed[i] = 99*MPH2MPS;
    
double nearest = 200;
    
// check each lane for nearby vehicles
int num_objs = sf.size();

for (int j = 0; j < num_objs; ++j){
      
  if (sf[j][SF.d] >= this->lines[i][0]) {
    
    // current object is to the right of the left lane line of ith lane
        
    if (sf[j][SF.d] <= this->lines[i][1]) {
          
      // current object is to the left of the right lane line of ith lane
          
      // predict object location in the future
      double vx = sf[j][SF.vx];
      double vy = sf[j][SF.vy];
      double v = sqrt(vx*vx + vy*vy);
      double s = sf[j][SF.s];
      s += v*TS*N;
          
      if (fabs(s - s_ego) < 10) {
		  
	     // the ith lane is blocked, no possiblity to switch here
        this->status[i] = "blocked";
        nearest = -999999;
        this->speed[i] = 0;
        this->distance[i] = 0;
            
      } else if ((fabs(s - s_ego) < nearest) && ((s - s_ego) > 0)) {
		  
	     // the ith lane has a vehicle somewhere on the horizon, but not immediately in the way
        this->speed[i] = v;
        this->distance[i] = s - s_ego;
        nearest = s - s_ego;
            
      }
    }
  }
}
```

The behavior is biased toward staying in the current lane.

```C++
/* 
cost[0] = cost to maintain lane
cost[1] = cost to increase lane (shift right)
cost[2] = cost to decrease lane (shift left)
*/
vector<double> cost = {0,0.1,0.1}; // biased against lane change
```

This bias is to prevent the car from constantly switching back and forth due to small numerical rounding errors.  As you can see in the video, there are still some unnecessary lane changes that could be removed by increasing the bias.

There are three major components to the cost for each action: velocity, follow distance, and obstruction.  If there is an obstacle in an adjacent lane, or if there is no adjacent lane (e.g. right lane shift from the right-most lane), the cost is set to an impossibly high value: `999`.

```C++
int num_lanes_right = (num_lanes - 1) - curr_lane;
int num_lanes_left = curr_lane;

  
if (num_lanes_right == 0) {
    
  // right shift not possible (no lane exists)
  cost[1] = 999;
    
} else if (lane.status[curr_lane+1] == "blocked") {
    
  // right shift not possible (obstacle)
  cost[1] = 999;
    
}
```

The latter of the three cost components is omitted for "keep current lane".

```C++
if (lane.distance[curr_lane] <= 100) {

  cost[0] += fmax(0.0, (100 - lane.distance[curr_lane])/200);
  
}

if (lane.speed[curr_lane] <= V_MAX) {

  cost[0] += fmax(0.0, (V_MAX - lane.speed[curr_lane]))/(2*V_MAX);
  
}
```

## Results

[![Whoops, there should be a picture here!](https://img.youtube.com/vi/FIU6CmpPsY0/0.jpg)](https://youtu.be/FIU6CmpPsY0)
