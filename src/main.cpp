#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

/* Custom variables */
#define MPH2MPS 0.44704                 // miles per hour to meters per second
const double TS = 0.02;                // software timestep seconds
const int N_TRAJ_PTS = int(1/TS);      // corresponds to 1 second time horizon
const double V_MAX_MPH = 49;            // miles per hour
double V_MAX= V_MAX_MPH * MPH2MPS; // meters per second

/* Sensor fusion indices */
const struct sensor_fusion_indices {
  int id = 0;
  int x = 1;
  int y = 2;
  int vx = 3;
  int vy = 4;
  int s = 5;
  int d = 6;
  int size = 7;
} SF;

/* Custom classes */

/*
 A class to keep track of the status of the road:
 where are the objects, how fast are the moving, etc
 */
class Lane
{
public:
  vector<string> status;        // text description of the lane status
  vector<vector<double>> lines; // an array of pairs, containing lane borders in 'd' coordinates
  vector<double> speed;         // maximum speed of each lane
  void init(int dimension, double lane_width);
  void update(vector<vector<double>> sf, double s_ego, double d_ego, int N);
  Lane();
};

Lane::Lane(void) {};

void Lane::init(int dimension, double lane_width){
  
  vector<string> ini;
  
  /*
   Initial state: assume no blockage, no speed limitation
   */
  for (int i = 0; i < dimension; ++i){
    
    ini.push_back("open");
    double left = double(i)*lane_width;
    double right = double(i+1)*lane_width;
    vector<double> lr_boundary;             // lane boundaries
    lr_boundary.push_back(left);            // left boundary for ith lane
    lr_boundary.push_back(right);           // right boundard for ith lane
    this->lines.push_back(lr_boundary);     // boundaries for ith lane
    this->speed.push_back(99*MPH2MPS);      // maximum speed for ith lane
    
  }
  
  // set init values
  this->status = ini;

}

void Lane::update(vector<vector<double>> sf, double s_ego, double d_ego, int N){
  
  int dim = this->status.size();
  for (int i = 0; i < dim; ++i) {
    
    // reset lane status
    this->status[i] = "open";
    this->speed[i] = 99*MPH2MPS;
    
    double nearest = 999999;
    
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
          
          if (fabs(s - s_ego) < 20) { // 3 car lengths between
		  
	          // the ith lane is blocked, no possiblity to switch here
            this->status[i] = "blocked";
            nearest = -999999;
            this->speed[i] = 0;
            
          } else if (((sf[j][SF.s] - s_ego) < nearest) && ((sf[j][SF.s] - s_ego) > 0)) {
		  
	          // the ith lane has a vehicle somewhere on the horizon, but not immediately in the way
            this->speed[i] = sqrt(sf[j][SF.vx]*sf[j][SF.vx] + sf[j][SF.vy]*sf[j][SF.vy]);
            
          }
        }
      }
    }
  }
}

/*
 A class to determine the optimal behavior of the vehicle
 */
class Behavior
{
public:
  string next; // text description of the next action to take
  void init(double car_d);
  void update(Lane lane, double car_d);
  double target_lane;
  Behavior();
};

Behavior::Behavior(void) {};

void Behavior::init(double car_d) {
  this->next = "KL"; // keep current lane
  this->target_lane = car_d;
}

void Behavior::update(Lane lane, double car_d) {

  /* 
   cost[0] = cost to maintain lane
   cost[1] = cost to increase lane (shift right)
   cost[2] = cost to decrease lane (shift left)
   */
  vector<double> cost = {0,0.03,0.03}; // biased against lane change

  /*
   Identify current lane
   */
  int num_lanes = lane.speed.size();
  int curr_lane;
  
  for (int i = 0; i < num_lanes; ++i) {
    
    if ((car_d >= lane.lines[i][0]) && (car_d <= lane.lines[i][1])) {
      
      // ego car is currently in the ith lane
      curr_lane = i;
      
    } else if (car_d < 0.0) {
      
      // prevent segmentation errors
      curr_lane = 0;
      
    } else if (car_d > num_lanes*4) {
      
      // prevent segmentation errors
      curr_lane = num_lanes-1;
      
    }
  }
  
  /*
   Check how many lane shifts are possible each direction
   */
  int num_lanes_right = (num_lanes - 1) - curr_lane;
  int num_lanes_left = curr_lane;

  
  if (num_lanes_right == 0) {
    
    // right shift not possible (no lane exists)
	  cost[1] = 999;
    
  } else if (lane.status[curr_lane+1] == "blocked") {
    
    // right shift not possible (obstacle)
	  cost[1] = 999;
    
  } else {
    
    double loc_cost = 999;
    
    for (int i = 1; i <= num_lanes_right; ++i) {
      
      double sub_cost;
      sub_cost = (V_MAX - lane.speed[curr_lane+i])/V_MAX;
      sub_cost = fmax(sub_cost, 0.0);
      loc_cost = fmin(loc_cost, sub_cost);
      
    }
    
    cost[1] += loc_cost;
    
  }
	
  if (num_lanes_left == 0) {
    
    // efltt shift not possible (no lane exists)
	  cost[2] = 999;
    
  } else if (lane.status[curr_lane-1] == "blocked") {
    
    // left shift not possible (obstacle)
	  cost[2] = 999;
    
  } else {
    
    double loc_cost = 999;
    
    for (int i = 1; i <= num_lanes_left; ++i) {
      
      double sub_cost;
      sub_cost = (V_MAX - lane.speed[curr_lane-i])/V_MAX;
      sub_cost = fmax(sub_cost, 0.0);
      loc_cost = fmin(loc_cost, sub_cost);
    }
    
    cost[2] += loc_cost;
    
  }
	  
  cost[0] +=  fmax(0.0, (V_MAX - lane.speed[curr_lane]))/V_MAX;
  
  if ((cost[1] < cost[2]) && (cost[1] < cost[0])) {
    
    this->next = "LCR";
    this->target_lane = (lane.lines[curr_lane+1][0] + lane.lines[curr_lane+1][1])/2;
    
  } else if ((cost[2] < cost[0]) && (cost[2] < cost[1])) {
    
    this->next = "LCL";
    this->target_lane = (lane.lines[curr_lane-1][0] + lane.lines[curr_lane-1][1])/2;
    
  } else {
    
    this->next = "KL";
    this->target_lane = (lane.lines[curr_lane][0] + lane.lines[curr_lane][1])/2;
    
  }
  
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;
  
  /* 
   Instances of custom class Lane to describe road
   and Behavior to decide ego vehicle behavior
   */
  Lane lane;
  Behavior behavior;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  
  
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane,&behavior](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {

          /*
           Simulator data
           */
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];
          auto sensor_fusion = j[1]["sensor_fusion"];
          
          car_speed *= MPH2MPS; // change to m/sec
          
          /*
           Static variables to carry through to following loop
           */
          static int loop;
          static bool init = false;
          
          /*
           Initialize custom classes first time through
           */
          if (!init) {
            behavior.init(car_d);
            init = true;
            lane.init(3, 4); //3 lanes, 4m wide
          }
          
          int num_prev = previous_path_x.size();
          vector<double> ptsx;
          vector<double> ptsy;
          
          /*
           Set current x, y, yaw as reference points for coordinate system transforms
           */
          vector<double> ref = {car_x, car_y, deg2rad(car_yaw)};
          
          if (num_prev < 2) {
            
            // not enough previous path data available, use current position
            // look backward using the heading
            double prev_car_x = car_x-cos(car_yaw);
            double prev_car_y = car_y-sin(car_yaw);
            
            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);
            
            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
            
          } else {
            
            // enough previous path data is available
            
            car_s = end_path_s;
            car_d = end_path_d;
            
            ref[0] = previous_path_x[num_prev-1];
            ref[1] = previous_path_y[num_prev-1];
            
            
            vector<double> ref_prev = {previous_path_x[num_prev-2], previous_path_y[num_prev-2], 0.0};
            double ref_yaw_prev = atan2(ref[1]-ref_prev[1], ref[0]-ref_prev[0]);
            ref_prev[2] = ref_yaw_prev;
            
            ptsx.push_back(ref_prev[0]);
            ptsx.push_back(ref[0]);
            
            ptsy.push_back(ref_prev[1]);
            ptsy.push_back(ref[1]);
            
          }
          
          /*
           Evaluate target lane every 2 sec
           */
          if (loop >= int(2/TS)) {
            
            /*
             Update the road obstacles and behavior
             */
            double target_lane = behavior.target_lane;
            lane.update(sensor_fusion, car_s, car_d, num_prev);
            
            // don't change lanes at very low speed
            if (car_speed > 10*MPH2MPS) {
              
              behavior.update(lane, car_d);
              
            }

            
            if (target_lane != behavior.target_lane) {
              std::cout << "Target lane = " << behavior.target_lane << std::endl;
            }
            
            loop = 0;
            
          } else {
            
            ++loop;
            
          }
          
          /*
           Set future waypoints
           based on decided target lane
           */
          double s_ahead = 30;
          vector<double> wp_next;
          
          wp_next = getXY(car_s+1*s_ahead, behavior.target_lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          ptsx.push_back(wp_next[0]);
          ptsy.push_back(wp_next[1]);
          
          wp_next = getXY(car_s+2*s_ahead, behavior.target_lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          ptsx.push_back(wp_next[0]);
          ptsy.push_back(wp_next[1]);
          
          wp_next = getXY(car_s+3*s_ahead, behavior.target_lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          ptsx.push_back(wp_next[0]);
          ptsy.push_back(wp_next[1]);
          
        
          /*
           Transform to car x, y coordinates from global x, y coordinates
           */
          int num_spline = ptsx.size();
          
          for (int i = 0; i < num_spline; ++i) {
            
            // shift backward
            double x_delta = ptsx[i] - ref[0];
            double y_delta = ptsy[i] - ref[1];
            
            // rotate backward (negative yaw)
            ptsx[i] = x_delta*cos(0 - ref[2]) - y_delta*sin(0-ref[2]);
            ptsy[i] = x_delta*sin(0 - ref[2]) + y_delta*cos(0-ref[2]);
            
          }
          
          /*
           Set trajectory using spline and look ahead points
           */
          tk::spline trajectory;
          trajectory.set_points(ptsx, ptsy);
          

          /*
           Find objects that are "nearby"
           */
          double follow_dist = car_speed*1.1;   // 1 second minimum follow distance + hysteresis band
          vector<vector<double>> near_objs;     // list of nearby cars (irrespective of lane
          int num_objs = sensor_fusion.size();  // all possible objects
          
          for (int i = 0; i < num_objs; ++i) {
            
            // if objects are withing follow distance (meters), consider them "near"
            double obj_vx  = sensor_fusion[i][SF.vx];
            double obj_vy  = sensor_fusion[i][SF.vy];
            double obj_v = sqrt(obj_vx*obj_vx + obj_vy*obj_vy);
            double obj_s = sensor_fusion[i][SF.s];
            obj_s += obj_v*num_prev*TS;
            
            double dst = obj_s - car_s;
            if (dst < follow_dist) {
              
              vector<double> obj;
              for (int j = 0; j < SF.size; j++) {
                
                obj.push_back(sensor_fusion[i][j]);
                
              }
              
              obj.push_back(dst);
              near_objs.push_back(obj);
              
            }
          }

          /*
           Set the vehicle speed to follow the nearest object in the same lane
           */
          double car_v_s_des = V_MAX;    // desired speed without accounting for traffic
          
          for (int i = 0; i < near_objs.size(); ++i) {
            
            // set the desired to the minimum of
            // nearby vehicles in front of the car
            double obj_vx  = near_objs[i][SF.vx];
            double obj_vy  = near_objs[i][SF.vy];
            double obj_v = sqrt(obj_vx*obj_vx + obj_vy*obj_vy);
            double obj_s = near_objs[i][SF.s] + obj_v*num_prev*TS;
            double obj_d = near_objs[i][SF.d];

            if (obj_s - car_s > 0) {
              
              // current object is ahead of the ego vehicle in s direction
              
              if (fabs(car_d - obj_d) < 2) {
                
                // current object is within 1/2 lane width
                
                double obj_v_s = (obj_s - car_s)/1.0; // desired spee to achieve 1 second following distance
                car_v_s_des = (obj_v_s < car_v_s_des) ? obj_v_s : car_v_s_des; // minimum speed for all blocking objects

              }
            }
          }
          
          /*
           Acceleration handling
           only allow 1 m/sec change per loop
           */
          double car_speed_next;
          
          if ((car_speed + 1) < car_v_s_des) {
            
            car_speed_next = car_speed + 1;
            
          } else if ((car_speed - 1) > car_v_s_des) {
            
            car_speed_next = car_speed - 1;
            
          } else {
            
            car_speed_next = car_speed;
            
          }
          
          json msgJson;
          
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          
          /*
           Push back previous values
           */
          for(int i = 0; i < num_prev; i++)
          {

              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            
          }
          
          /*
           Break velocity up into orthogonal components to not exceed speed limit
           */
          double x_des = car_v_s_des; //30;
          double y_des = trajectory(x_des);
          double dist = sqrt(x_des*x_des + y_des*y_des);
          
          /*
           Fill in remaining points of trajectory
           */
          int num_remain = N_TRAJ_PTS - num_prev;
          
          for (int i = 1; i < num_remain; ++i) {
            
            double N = dist/(TS*car_speed_next);
            double x_new = i*x_des/N;
            double y_new = trajectory(x_new);
            
            double x_new_ = x_new; // dummy variable
            
            // switch back to global coordinates
            x_new = x_new*cos(ref[2]) - y_new*sin(ref[2]);
            y_new = x_new_*sin(ref[2]) + y_new*cos(ref[2]);
            
            x_new += ref[0];
            y_new += ref[1];
            
            next_x_vals.push_back(x_new);
            next_y_vals.push_back(y_new);
            
          }
          
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;
          
          auto msg = "42[\"control\","+ msgJson.dump()+"]";
          
          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
