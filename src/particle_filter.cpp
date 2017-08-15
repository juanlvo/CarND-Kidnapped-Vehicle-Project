/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

/**
 * init Initializes particle filter by initializing particles to Gaussian
 *   distribution around first position and all the weights to 1.
 * @param x Initial x position [m] (simulated estimate from GPS)
 * @param y Initial y position [m]
 * @param theta Initial orientation [rad]
 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 50;
	double initial_weight = 1.0;

	default_random_engine gen;

	// This line creates a normal (Gaussian) distribution
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// TODO: Sample  and from these normal distrubtions like this:
	//	 sample_x = dist_x(gen);
	//	 where "gen" is the random engine initialized earlier.


	for (int i = 0; i < num_particles; ++i) {

		Particle p = {
			i,
			dist_x(gen),
			dist_y(gen),
			dist_theta(gen),
			initial_weight
		};

	    weights.push_back(initial_weight);
	    particles.push_back(p);
	}

	is_initialized = true;
	return;
}

/**
 * prediction Predicts the state for the next time step
 *   using the process model.
 * @param delta_t Time between time step t and t+1 in measurements [s]
 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 * @param velocity Velocity of car from t to t+1 [m/s]
 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
 */
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	for (int i=0; i < num_particles; i++) {

		// This line creates a normal (Gaussian) distribution for sensor noise
		normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

		if (fabs(yaw_rate) < 0.00001) {
		      particles[i].x += velocity * delta_t * cos(particles[i].theta);
		      particles[i].y += velocity * delta_t * sin(particles[i].theta);
		} else {
		      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
		      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
		      particles[i].theta += yaw_rate * delta_t;
		}

		default_random_engine gen;

	    // add noise
	    particles[i].x += dist_x(gen);
	    particles[i].y += dist_y(gen);
	    particles[i].theta += dist_theta(gen);
	}

}

/**
 * dataAssociation Finds which observations correspond to which landmarks (likely by using
 *   a nearest-neighbors data association).
 * @param predicted Vector of predicted landmark observations
 * @param observations Vector of landmark observations
 */
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i=0; i < observations.size(); i++) {

		//current observation
		LandmarkObs lo = observations[i];

		// minimum distance to maximum possible
		double min_dist = numeric_limits<double>::max();

		// init id of landmark from map to be associated with the observation
		int map_id = -1;

		for (int j=0; j < predicted.size(); j++) {

			//grab current prediction
			LandmarkObs p = predicted[j];

			//get the distance between current/predicted landmarks
			double cur_dist = dist(lo.x, lo.y, p.x, p.y);

			//find the predicted landmark nearest the current observed landmark
			if (cur_dist < min_dist) {
				min_dist = cur_dist;
				map_id = p.id;
			}
		}

		//set the observation id to the nearest
		observations[i].id = map_id;
	}

	return;

}

/**
 * updateWeights Updates the weights for each particle based on the likelihood of the
 *   observed measurements.
 * @param sensor_range Range [m] of sensor
 * @param std_landmark[] Array of dimension 2 [standard deviation of range [m],
 *   standard deviation of bearing [rad]]
 * @param observations Vector of landmark observations
 * @param map Map class containing map landmarks
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int i=0; i < observations.size(); i++) {

		//get the coordinates of the particles
		double pX = particles[i].x;
		double pY = particles[i].y;
		double pTheta = particles[i].theta;

		//create a vector with landmarks locations predicted
		vector<LandmarkObs> predictions;

		//for each map landmark
		for (int j=0; j < map_landmarks.landmark_list.size(); j++) {

			//get id and x, y
			float lmX = map_landmarks.landmark_list[j].x_f;
			float lmY = map_landmarks.landmark_list[j].y_f;
			int lmId = map_landmarks.landmark_list[j].id_i;

			//landmarks in the sensor range of the particle
			if (fabs(lmX-pX) <= sensor_range && fabs(lmY-pY) <= sensor_range) {

				//add the prediction
				predictions.push_back(LandmarkObs{lmId, lmX, lmY});
			}
		}

		//create and populate of the list of observations transformed to map coordinates
		vector<LandmarkObs> transformedOs;

		for (int j = 0; j < observations.size(); j++) {

		      double tX = cos(pTheta)*observations[j].x - sin(pTheta)*observations[j].y + pX;
		      double tY = sin(pTheta)*observations[j].x + cos(pTheta)*observations[j].y + pY;
		      transformedOs.push_back(LandmarkObs{ observations[j].id, tX, tY });

		}

		//execute data association for the predictions and transformed  observations on current particle
		dataAssociation(predictions, transformedOs);

		//init again the weight
		particles[i].weight = 1.0;

		for (int j = 0; j <  transformedOs.size(); j++) {

		      //observation and associated prediction coordinates
		      double oX, oY, prX, prY;
		      oX = transformedOs[j].x;
		      oY = transformedOs[j].y;

		      int associatedPrediction = transformedOs[j].id;

		      // get the x,y coordinates of the prediction associated with the current observation
		      for (unsigned int k = 0; k < predictions.size(); k++) {
		        if (predictions[k].id == associatedPrediction) {
		          prX = predictions[k].x;
		          prY = predictions[k].y;
		        }
		      }

		      // calculate weight for this observation with multivariate Gaussian
		      double sX = std_landmark[0];
		      double sY = std_landmark[1];
		      double obsW = ( 1/(2*M_PI*sX*sY)) * exp( -( pow(prX-oX,2)/(2*pow(sX, 2)) + (pow(prY-oY,2)/(2*pow(sY, 2))) ) );

		      // product of this obersvation weight with total observations weight
		      particles[i].weight *= obsW;

		}
	}

	return ;
}

/**
 * resample Resamples from the updated set of particles to form
 *   the new set of particles.
 */
void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> particlesNew;

	default_random_engine gen;

	//read and assign all the weights
	for (int i=0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
	}

	//generate random starting index for resampling wheel
	uniform_int_distribution<int> uniintdist(0, num_particles-1);
	auto index = uniintdist(gen);

	// get max weight
	double maxWeight = *max_element(weights.begin(), weights.end());

	// uniform random distribution [0.0, max_weight)
	uniform_real_distribution<double> unirealdist(0.0, maxWeight);

	double beta = 0.0;

	// spin the resample wheel!
	for (int i = 0; i < num_particles; i++) {
	beta += unirealdist(gen) * 2.0;
	while (beta > weights[index]) {
	  beta -= weights[index];
	  index = (index + 1) % num_particles;
	}
	particlesNew.push_back(particles[index]);
	}

	particles = particlesNew;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
