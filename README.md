# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

For the controller/trajectory portion of this project, I relied on the guidance from the project walkthrough.  Prior to the walkthrough I tried many different iterations using JMT, and was unable to create a smooth transition between old path and new path.  As in the walkthrough, JMT is not used, in favor of a trajectory spline.

## Behavior

I created two custom classes to determine the behavior of the ego vehicle: `Lane` and `Behavior`.  `Lane` describes what is going on in each of the lanes: feasible speed, following distance, and obstruction.  `Behavior` takes `Lane lane` as an input, and applies a cost function to determine the target lane.

The behavior is biased toward staying in the current lane:

```C++
/* 
   cost[0] = cost to maintain lane
   cost[1] = cost to increase lane (shift right)
   cost[2] = cost to decrease lane (shift left)
   */
  vector<double> cost = {0,0.1,0.1}; // biased against lane change
```
   
## Results

[![Whoops, there should be a picture here!](https://img.youtube.com/vi/FIU6CmpPsY0/0.jpg)](https://youtu.be/FIU6CmpPsY0)
