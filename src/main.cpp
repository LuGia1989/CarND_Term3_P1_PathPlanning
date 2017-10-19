#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

#define PROJECT_STEP 50  // how far to project into the future

#define SAFETY_ZONE 20   // how close is too close

#define SPEED_LIMIT 49.5 // Speed limit in MPH, multiply by 2.24 to convert to m/s
#define DT 0.02          // 20 ms in s
#define LARGE_NUM 10000

#define ADJUST_VELOCITY 0.200

using namespace std;
using json = nlohmann::json; // for convenience

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

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
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
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

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
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

double CheckSpeed(double vx, double vy) {
	return (sqrt(vx*vx + vy*vy)*DT);
}


double Gap(double vx, double vy, int prev_path_size, double s, double end_path_s) {
	double check_speed = sqrt(vx*vx + vy*vy)*DT;
	return (s + ((double)prev_path_size * check_speed) - end_path_s);
}


double TargetSpeed (bool watchout_ahead, double speed) {
	double target_speed = speed;
	if (watchout_ahead) {
		target_speed -= ADJUST_VELOCITY;
	}
	else {
		target_speed += ADJUST_VELOCITY;
	}
	 if (target_speed > SPEED_LIMIT) {
		 target_speed = SPEED_LIMIT;
	 }
	if (target_speed < 0) {
		target_speed = 0.0;
	}
	return target_speed;
}

double AdjustCarVelocity (double target_speed, double car_ref_speed) {
	double car_velocity = car_ref_speed;
	if (car_velocity < ADJUST_VELOCITY) {car_velocity = ADJUST_VELOCITY;}
	else if ((car_velocity + ADJUST_VELOCITY) <= target_speed) {car_velocity += ADJUST_VELOCITY;}
	else {car_velocity -= ADJUST_VELOCITY;}
	return car_velocity;
}



double WayPointx(double x, double yaw) {
	return (x - cos(yaw));
}

double WayPointy(double y, double yaw) {
	return (y - sin(yaw));
}



int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;


  int lane = 1; // Lanes are numbered 0 - 2 with 0 adjacent to the centerline of the road
  int plan_steps = PROJECT_STEP;
  int safety_zone = SAFETY_ZONE;

  double car_ref_velocity = 0.0; // Car velocity in MPH, multiply by 2.24 to convert to m/s


  string map_file_ = "../data/highway_map.csv"; // Waypoint map to read from

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

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &lane, &safety_zone, &plan_steps, &car_ref_velocity](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
	  //double car_speed_adjust = ADJUST_VELOCITY;
	  double left_lane_open = LARGE_NUM;
	  double right_lane_open = LARGE_NUM;

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
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
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

            int previous_path_size = previous_path_x.size();

            if (previous_path_size == 0) {
            	 end_path_s = car_s;
            }

            double safety_zone = (car_ref_velocity / 2) + 5;
            bool watchout_ahead = false;
            bool change_left_lane = true;
            bool change_right_lane = true;

            if (lane == 0) {change_left_lane = false;}
            if (lane == 2) {change_right_lane = false;}

            for (int i=0; i<sensor_fusion.size(); i++) { // for each car
              float sf_d = sensor_fusion[i][6];
              float d_lane = 1;
              if (sf_d < 4) {d_lane = 0;}
              if (sf_d > 8) {d_lane = 2;}

              double sf_vx = sensor_fusion[i][3];
              double sf_vy = sensor_fusion[i][4];
              double sf_s = sensor_fusion[i][5];
              double check_speed = CheckSpeed(sf_vx, sf_vy);
              double project_s = sf_s + ((double)previous_path_size * check_speed);
              double gap = Gap(sf_vx, sf_vy, previous_path_size, sf_s, end_path_s);

              if (sf_s > car_s) {
                if (gap < safety_zone) {
                  if (d_lane == lane) {
                	  watchout_ahead = true;
                  }

                  if (d_lane == (lane-1)) {
                	  change_left_lane = false;
                      if (gap < left_lane_open) {left_lane_open = gap;}
                  }
                  if (d_lane == (lane+1)) {
                	  change_right_lane = false;
                    if (gap < right_lane_open) {right_lane_open = gap;}
                  }
                }
              }
              if ((sf_s > car_s) && (sf_s < end_path_s)) {
                if (d_lane == (lane-1)) {change_left_lane = false;}
                if (d_lane == (lane+1)) {change_right_lane = false;}
              }
              if ((project_s > car_s) && (project_s < end_path_s)) {
                if (d_lane == (lane-1)) {change_left_lane = false;}
                if (d_lane == (lane+1)) {change_right_lane = false;}
              }
            }

            if (watchout_ahead) {
              if (change_left_lane) {
                if (change_right_lane && (right_lane_open > left_lane_open)) {lane += 1;}
                else {lane -= 1;}
              }
              else if (change_right_lane) {lane += 1;}
            }

          	json msgJson;

            vector<double> ptsx; // Waypoint based anchor points
            vector<double> ptsy; // Waypoint based anchor points
	    
	    double ref_yaw = deg2rad(car_yaw);
            double ref_x = car_x;
            double ref_y = car_y;

            if (previous_path_size >= 2) { // use two previous path points
              ptsx.push_back(previous_path_x[previous_path_size-2]);
              ptsy.push_back(previous_path_y[previous_path_size-2]);


              ref_x = previous_path_x[previous_path_size-1];
              ref_y = previous_path_y[previous_path_size-1];
            }
            else { // if there aren't enough path points
              double ref_yaw = deg2rad(car_yaw);
              ptsx.push_back(WayPointx(car_x, car_yaw));
              ptsy.push_back(WayPointy(car_y, car_yaw));
            }
            ptsx.push_back(ref_x);
            ptsy.push_back(ref_y);

            // generate path points
            for (int i=safety_zone; i<=(3*safety_zone); i+=safety_zone) {
              vector<double> next_wp = getXY(end_path_s+i,(lane*4)+2, map_waypoints_s, map_waypoints_x, map_waypoints_y);
              ptsx.push_back(next_wp[0]);
              ptsy.push_back(next_wp[1]);
            }

            // move points to car coordinates
            for (int i=0; i<ptsx.size(); i++) {
              double shift_x = ptsx[i] - ref_x;
              double shift_y = ptsy[i] - ref_y;

              ptsx[i] = ((shift_x)*cos(0-ref_yaw) - (shift_y)*sin(0-ref_yaw));
              ptsy[i] = ((shift_x)*sin(0-ref_yaw) + (shift_y)*cos(0-ref_yaw));
            }

            // create and initialize spline
            tk::spline spl; // create spline
            spl.set_points(ptsx, ptsy); // add anchor points to spline

            vector<double> next_x_vals; // Planner path points
            vector<double> next_y_vals; // Planner path points

            for (int i=0; i<previous_path_size; i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            double target_x = safety_zone;
            double target_y = spl(target_x);
            double target_velocity = car_ref_velocity;
            double current_velocity = car_ref_velocity;
            double x_point = ptsx[1];
            double y_point = ptsy[1];

            int waypoints = plan_steps-previous_path_size;
            for (int i=1; i <= waypoints; i++) { // fill remaining points

              // Select target speed
              current_velocity = car_ref_velocity;
              double speed = TargetSpeed(watchout_ahead, target_velocity);
              target_velocity = speed;
              car_ref_velocity = AdjustCarVelocity(target_velocity, car_ref_velocity);

              x_point += target_x/((sqrt(target_x*target_x + target_y*target_y))/((car_ref_velocity/2.24)*DT));
              y_point = spl(x_point);

              next_x_vals.push_back((x_point * cos(ref_yaw) - y_point * sin(ref_yaw)) + ref_x);
              next_y_vals.push_back((x_point * sin(ref_yaw) + y_point * cos(ref_yaw)) + ref_y);
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
