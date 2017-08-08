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
#define MPH2MPS 0.44704                 // miles per hour to meters per second
const double V_MAX_MPH = 47.5;            // miles per hour
double V_MAX_MPS = V_MAX_MPH * MPH2MPS; // meters per second
double A_MAX_MPS2 = 5;                  // meters per second per second
double J_MAX_MPS2 = 10;                  // meters per second per second
const double TS = 0.02;                 // software timestep seconds
const int N_TRAJ_PTS = 50;              // number of trajectory points
const int NUM_PREV_MAX = 10;            // maximum number of previous trajectory points to use

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
  vector<string> status; // text description of the lane status
  vector<vector<double>> lines; // an array of pairs, containing lane borders in 'd' coordinates
  vector<double> speed;  // maximum speed of each lane
  void init(int dimension, double lane_width);
  void update(vector<vector<double>> sf, double s_ego, double d_ego);
  Lane();
};

Lane::Lane(void) {};

void Lane::init(int dimension, double lane_width){
  vector<string> ini;
  
  // initial state: assume no blockage, no speed limitation
  for (int i = 0; i < dimension; ++i){
    ini.push_back("open");
    double left = double(i)*lane_width;
    double right = double(i+1)*lane_width;
    vector<double> lr_boundary; // lane boundaries
    lr_boundary.push_back(left);
    lr_boundary.push_back(right);
    this->lines.push_back(lr_boundary);
    this->speed.push_back(99*MPH2MPS);
  }
  this->status = ini;

}

void Lane::update(vector<vector<double>> sf, double s_ego, double d_ego){
  
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
        if (sf[j][SF.d] <= this->lines[i][1]) {
          if (fabs(sf[j][SF.s] - s_ego) < 10) {
		  
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
    
    /* debug */
    //std::cout << "Lane " << i << " is " << this->status[i] << " - " << this->lines[i][0] << "m to " << this->lines[i][1] << "m" << std::endl;
    //std::cout << "     Speed is " << this->speed[i]/MPH2MPS << " MPH" << std::endl;
  }
}

/*
 A class to determine the optimal behavior of the vehicle
 */
class Behavior
{
public:
  string next; // text description of the next action to take
  void init(void);
  void update(Lane lane, double car_d);
  double target_lane;
  Behavior();
};

Behavior::Behavior(double car_d) {};

void Behavior::init(void) {
  this->next = "KL"; // keep current lane
  this->target_lane = car_d;
}

void Behavior::update(Lane lane, double car_d) {
  // calculate cost and decide on action
  /* 
   cost[0] = cost to maintain lane
   cost[1] = cost to increase lane (shift right)
   cost[2] = cost to decrease lane (shift left)
   */
  vector<double> cost = {0,0.03,0.03}; // biased against lane change

  // identify current lane
  int num_lanes = lane.speed.size();
  //std::cout << "Number of lanes: " << num_lanes << std::endl;
  int curr_lane;
  for (int i = 0; i < num_lanes; ++i) {
    if ((car_d >= lane.lines[i][0]) && (car_d <= lane.lines[i][1])) {
      curr_lane = i;
    } else if (car_d < 0.0) {
      curr_lane = 0;
    } else if (car_d > num_lanes*4) {
      curr_lane = num_lanes-1;
    }
  }
  
  int num_lanes_right = (num_lanes - 1) - curr_lane;
  int num_lanes_left = curr_lane;
  
  /* debug */
  //std::cout << "Current lane: " << curr_lane << std::endl;
  //std::cout << "Lanes to the right: " << num_lanes_right << std::endl;
  //std::cout << "Lanes to the left: " << num_lanes_left << std::endl;
  
  if (num_lanes_right == 0) {
	  cost[1] = 999;
  } else if (lane.status[curr_lane+1] == "blocked") {
	  cost[1] = 999;
  } else {
    
    double loc_cost = 999;
    for (int i = 1; i <= num_lanes_right; ++i) {
      double sub_cost;
      sub_cost = (V_MAX_MPS - lane.speed[curr_lane+i])/V_MAX_MPS;
      sub_cost = fmax(sub_cost, 0.0);
      loc_cost = fmin(loc_cost, sub_cost);
    }
    
    cost[1] += loc_cost;
    
  }
	
  if (num_lanes_left == 0) {
	  cost[2] = 999;
  } else if (lane.status[curr_lane-1] == "blocked") {
	  cost[2] = 999;
  } else {
    
    double loc_cost = 999;
    for (int i = 1; i <= num_lanes_left; ++i) {
      double sub_cost;
      sub_cost = (V_MAX_MPS - lane.speed[curr_lane-i])/V_MAX_MPS;
      sub_cost = fmax(sub_cost, 0.0);
      loc_cost = fmin(loc_cost, sub_cost);
    }
    
    cost[2] += loc_cost;
    
  }
	  
  cost[0] +=  fmax(0.0, (V_MAX_MPS - lane.speed[curr_lane]))/V_MAX_MPS;

  /* debug */
  //std::cout << "Cost to keep: " << cost[0] << std::endl;
  //std::cout << "Cost to change right: " << cost[1] << std::endl;
  //std::cout << "Cost to change left: " << cost[2] << std::endl;
  
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

// quintic polynomial solution (minimum jerk)
vector<double> JMT(vector< double> start, vector <double> end, double T)
{
  
  double a0 = start[0];
  double a1 = start[1];
  double a2 = start[2]*0.5;
  
  Eigen::MatrixXd A(3,3);
  A(0,0) = T*T*T;
  A(0,1) = T*T*T*T;
  A(0,2) = T*T*T*T*T;
  A(1,0) = 3*T*T;
  A(1,1) = 4*T*T*T;
  A(1,2) = 5*T*T*T*T;
  A(2,0) = 6*T;
  A(2,1) = 12*T*T;
  A(2,2) = 20*T*T*T;
  
  Eigen::VectorXd b(3);
  b(0) = end[0] - (start[0] + start[1]*T + 0.5*start[2]*T*T);
  b(1) = end[1] - (start[1] + start[2]*T);
  b(2) = end[2] - start[2];
  
  Eigen::MatrixXd Ainv = A.inverse();
  
  Eigen::VectorXd x = Ainv*b;
  double a3 = x(0);
  double a4 = x(1);
  double a5 = x(2);
  
  return {a0,a1,a2,a3,a4,a5};
  
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
   instance of custom class Lane to describe road
   and ego vehicle behavior
   */
  Lane lane;
  lane.init(3, 4); //3 lanes, 4m wide
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
  
  
  /*
   Define new waypoint maps with high resolution (every 0.5m in the s domain)
   */
  
  int num_upsampled_waypoints = int(max_s*2); // place a waypoint every 0.5m
  
  tk::spline upsample_waypoints_s_x;
  tk::spline upsample_waypoints_s_y;
  tk::spline upsample_waypoints_s_dx;
  tk::spline upsample_waypoints_s_dy;
  
  upsample_waypoints_s_x.set_points(map_waypoints_s, map_waypoints_x);
  upsample_waypoints_s_y.set_points(map_waypoints_s, map_waypoints_y);
  upsample_waypoints_s_dx.set_points(map_waypoints_s, map_waypoints_dx);
  upsample_waypoints_s_dy.set_points(map_waypoints_s, map_waypoints_dy);
  
  vector<double> upsample_waypoints_x;
  vector<double> upsample_waypoints_y;
  vector<double> upsample_waypoints_s;
  vector<double> upsample_waypoints_dx;
  vector<double> upsample_waypoints_dy;
  
  for (int i = 0; i < num_upsampled_waypoints; ++i) {
    double upsample_s = double(i)/2.0;
    upsample_waypoints_x.push_back(upsample_waypoints_s_x(upsample_s));
    upsample_waypoints_y.push_back(upsample_waypoints_s_y(upsample_s));
    upsample_waypoints_dx.push_back(upsample_waypoints_s_dx(upsample_s));
    upsample_waypoints_dy.push_back(upsample_waypoints_s_dy(upsample_s));
    upsample_waypoints_s.push_back(upsample_s);
  }
  
  h.onMessage([&upsample_waypoints_x,&upsample_waypoints_y,&upsample_waypoints_s,&upsample_waypoints_dx,&upsample_waypoints_dy,&upsample_waypoints_s_x,&upsample_waypoints_s_y,&upsample_waypoints_s_dx,&upsample_waypoints_s_dy,&lane,&behavior](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	//double car_s = j[1]["s"];
          	//double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
          
          // get better Frenet coordinates
          vector<double> frenet = getFrenet(car_x, car_y, car_yaw, upsample_waypoints_x, upsample_waypoints_y);
          double car_s = frenet[0];
          double car_d = frenet[1];
          double vel_ini = car_speed*MPH2MPS;
		
	  static bool init = false;
		
	  if (!init) {
	    behavior.init(car_d);	  
	  }
          
	  /*debug*/
	  //std::cout << vel_ini << std::endl;
          
          /* debug */
          //std::cout << sensor_fusion << std::endl;
          //std::cout << sensor_fusion[0] << std::endl;
          //std::cout << sensor_fusion.size() << std::endl;
          
          /* debug */
          /*
          for (int i = 0; i < sensor_fusion.size(); ++i) {
            std::cout << sensor_fusion[i] << std::endl;
          }
           */
          
          /* debug */
          
          /* low speed output to terminal */
          
          int num_prev = previous_path_x.size();
          
          std::cout << "Previous s" << std::endl;
          if (num_prev > 0) {
            for (int i = 0; i < num_prev; ++i) {
              vector<double> s_frenet_prev = getFrenet(previous_path_x[i], previous_path_y[i], car_yaw, upsample_waypoints_x, upsample_waypoints_y);
              double s_prev = s_frenet_prev[0];
              std::cout << s_prev << ", ";
            }
          }
          
          std::cout << std::endl;
          std::cout << std::endl;
          
          static int loop;
          
          static vector<double> jmt_d = {car_d, 0, 0, 0, 0, 0};
          vector<double> d_i;
          vector<double> d_f;
          
          double tgt_lane = car_d;
          static string behavior_old = "KL";
          
          if (loop >= int(0.5/TS)) {
            loop = 0;
            
            int curr_lane;
            // display current lane info
            if ((car_d > 12) || (car_d < 0)) {
              curr_lane = -1;
              std::cout << "Not in a valid lane!" << std::endl;
            } else if (car_d > 8) {
              curr_lane = 0;
              std::cout << "In right lane!" << std::endl;
            } else if (car_d > 4) {
              curr_lane = 1;
              std::cout << "In center lane!" << std::endl;
            } else {
              curr_lane = 2;
              std::cout << "In left lane!" << std::endl;
            }
            lane.update(sensor_fusion, car_s, car_d);
            
            
            behavior.update(lane, car_d);
            if (behavior_old != behavior.next) {
              std::cout << "Change lanes!" << std::endl;
              behavior_old = behavior.next;
            }
            

          } else {
            ++loop;
            behavior_old = "KL";
          }
          
          vector<double> d_trajectory(50);
          
          for (int i = 0; i < 50; ++i) {
            
          }
          
          
          // open loop d trajectory estimate
          double max_d_per_sec = 1; // 2 m/sec maximum lateral change
          
          double d_fin;
          double d_vel_fin;
          
          if ((behavior.target_lane - car_d) > max_d_per_sec) {
            d_fin = car_d + max_d_per_sec;
            d_vel_fin = max_d_per_sec;
          } else if (( behavior.target_lane - car_d) < -max_d_per_sec) {
            d_fin = car_d - max_d_per_sec;
            d_vel_fin = -max_d_per_sec;
          } else {
            d_fin = behavior.target_lane;
            d_vel_fin = 0.0;
          }
          
          double d_vel_ini;
          if (num_prev > 1) {
            vector<double> d_frenet_0 = getFrenet(previous_path_x[0], previous_path_y[0], car_yaw, upsample_waypoints_x, upsample_waypoints_y);
            double car_d_0 = d_frenet_0[1];
            vector<double> d_frenet_1 = getFrenet(previous_path_x[1], previous_path_y[1], car_yaw, upsample_waypoints_x, upsample_waypoints_y);
            double car_d_1 = d_frenet_1[1];
            d_vel_ini = (car_d_1-car_d_0)/TS;
          } else {
            d_vel_ini = 0.0;
          }
          
          d_i = {car_d, d_vel_ini, 0};
          d_f = {d_fin, d_vel_fin, 0};
          jmt_d = JMT(d_i, d_f, 1);
            
          

          
          
          double follow_dist = vel_ini*1.1; // 1 second minimum follow distance + hysteresis band
          
          vector<vector<double>> near_objs;
          int num_objs = sensor_fusion.size();
          for (int i = 0; i < num_objs; ++i) {
            // if objects are withing DIST_MAX meters, consider them "near"
            double obj_x = sensor_fusion[i][SF.x];
            double obj_y = sensor_fusion[i][SF.y];
            double dst_x = car_x - obj_x;
            double dst_y = car_y - obj_y;
            double dst = sqrt((dst_x*dst_x + dst_y*dst_y));
            if (dst < follow_dist) {
              // push back is broken up in case additional fields are added
              vector<double> obj;
              for (int j = 0; j < SF.size; j++) {
                obj.push_back(sensor_fusion[i][j]);
              }
              obj.push_back(dst);
              near_objs.push_back(obj);
            }
          }
          

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
          
          // follow the car in front
          double vel_des = V_MAX_MPS;
          
          for (int i = 0; i < near_objs.size(); ++i) {
            double obj_x = near_objs[i][SF.x];
            double obj_y = near_objs[i][SF.y];
            vector<double> obj_frenet = getFrenet(obj_x, obj_y, car_yaw, upsample_waypoints_x, upsample_waypoints_y);
            double obj_s = obj_frenet[0];
            double obj_d = obj_frenet[1];
            if (obj_s - car_s > 0) {
              if (fabs(car_d - obj_d) < 2) {
                double vel = (obj_s - car_s)/1.0; // 1 second following distance
                if (vel < vel_des) {
                  vel_des = vel;
                }
                //std::cout << "Follow car #" << near_objs[i][SF.id] << " at " << double(int(10*vel/MPH2MPS))/10 << "MPH" << std::endl;
              }
            }
          }
          
          // jerk minimization in s domain
          
          
          /* calculate s trajectory */
          
          vector<double> s_trajectory;
          
          double s_ini, s_v_ini, s_a_ini;
          double s_fin, s_v_fin, s_a_fin;
          double s_a_max = 10;
          int s_to_fill;
        
          if (num_prev > 0) {
            vector<double> s_frenet_0 = getFrenet(previous_path_x[0], previous_path_y[0], car_yaw, upsample_waypoints_x, upsample_waypoints_y);
            s_ini = s_frenet_0[0];
            s_v_ini = 0.0;
            s_a_ini = 0.0;
            s_to_fill = 49;
            s_trajectory.push_back(s_frenet_0[0]);
          } else {
            vector<double> s_frenet_0 = getFrenet(car_x, car_y, car_yaw, upsample_waypoints_x, upsample_waypoints_y);
            s_ini = s_frenet_0[0];
            s_v_ini = 0.0;
            s_a_ini = 0.0;
            s_to_fill = 49;
            s_trajectory.push_back(s_frenet_0[0]);
          }
          
          if (vel_des > (s_v_ini + s_a_max)) {
            s_v_fin = s_v_ini + s_a_max;
            s_fin = s_ini + (s_v_fin + s_v_ini)/2;
            s_a_fin = 0.1;
          } else if (vel_des < (s_v_ini - s_a_max)) {
            s_v_fin = s_v_ini - s_a_max;
            s_fin = s_ini + (s_v_fin + s_v_ini)/2;
            s_a_fin = -0.1;
          } else {
            s_v_fin = vel_des;
            s_fin = s_ini + (s_v_fin + s_v_ini)/2;
            s_a_fin = 0.1*(s_v_fin-s_v_ini)/fabs(s_v_fin-s_v_ini);
          }
          
          double s_f_T = 1;
          
          vector<double> jmt_s;
          vector<double> s_i;
          vector<double> s_f;
          
          s_i = {s_ini, s_v_ini, 0};
          s_f = {s_fin, s_v_fin, 0};
          jmt_s = JMT(s_i, s_f, s_f_T);
          
          for (int i = 1; i < 50; ++i) {
            double ts = (i)*TS;
            double s0 = jmt_s[0];
            double s1 = jmt_s[1]*ts;
            double s2 = jmt_s[2]*ts*ts;
            double s3 = jmt_s[3]*ts*ts*ts;
            double s4 = jmt_s[4]*ts*ts*ts*ts;
            double s5 = jmt_s[5]*ts*ts*ts*ts*ts;
            double next_s = s0 + s1 + s2 + s3 + s4 + s5;
            //s_trajectory.push_back(next_s);
          }
          
          double new_vel;
          if (fabs(vel_ini - vel_des) < 0.001) {
            new_vel = vel_des;
          } else if (vel_ini < vel_des){
            new_vel = vel_ini + 0.001;
          } else {
            new_vel = vel_ini - 0.001;
          }
          
          for (int i = 1; i < 50; ++i) {
            double last_s = s_trajectory[i-1];
            double new_s = last_s + new_vel*TS;
            s_trajectory.push_back(new_s);
          }

          static double next_d = 6;
          static double prev_d = 6;
          
          std::cout << "Next s" << std::endl;
          for(int i = 0; i < N_TRAJ_PTS; i++)
          {
            
            if ((i == 0)&&(num_prev > 0)) {
              next_x_vals.push_back(previous_path_x[0]);
              next_y_vals.push_back(previous_path_y[0]);
            } else if (i == 0) {
              next_x_vals.push_back(car_x);
              next_y_vals.push_back(car_y);
            } else{
          
              double next_s = s_trajectory[i];
              double prev_s = s_trajectory[i-1];
            
              std::cout << next_s << ", ";
            
              double next_x = upsample_waypoints_s_x(next_s) + upsample_waypoints_s_dx(next_s)*next_d;
              double next_y = upsample_waypoints_s_y(next_s) + upsample_waypoints_s_dy(next_s)*next_d;
            
              double prev_x = upsample_waypoints_s_x(prev_s) + upsample_waypoints_s_dx(prev_s)*prev_d;
              double prev_y = upsample_waypoints_s_y(prev_s) + upsample_waypoints_s_dy(prev_s)*prev_d;
              double next_x_delta = next_x_vals[i-1] + (next_x - prev_x);
              double next_y_delta = next_y_vals[i-1] + (next_y - prev_y);
              next_x_vals.push_back(next_x_delta);
              next_y_vals.push_back(next_y_delta);
            }
      

          }
          
          std::cout << std::endl;
          std::cout << std::endl;
          
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
















































































