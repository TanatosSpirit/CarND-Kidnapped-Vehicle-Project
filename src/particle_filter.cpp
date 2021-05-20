/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1;  // TODO: Set the number of particles

    std::default_random_engine gen;

    // This lines creates a normal (Gaussian) distributions
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_theta(theta, std[2]);

    double sample_x, sample_y, sample_theta;

    for (int i = 0; i < num_particles; ++i) {
        sample_x = dist_x(gen);
        sample_y = dist_y(gen);
        sample_theta = dist_theta(gen);
        Particle particle{i, sample_x, sample_y, sample_theta, 1.0/num_particles};
        particles.push_back(particle);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

    std::default_random_engine gen;

    for (int i = 0; i < num_particles; ++i) {

        double theta_0 = particles[i].theta;
        double x_0 = particles[i].x;
        double y_0 = particles[i].y;

        double x_new = x_0 + (velocity/yaw_rate) * (sin(theta_0 + yaw_rate * delta_t) - sin(theta_0));
        double y_new = y_0 + (velocity/yaw_rate) * (cos(theta_0) - cos(theta_0 + yaw_rate * delta_t));
        double theta_new = theta_0 + yaw_rate * delta_t;

        // This lines creates a normal (Gaussian) distributions
        std::normal_distribution<double> dist_x(x_new, std_pos[0]);
        std::normal_distribution<double> dist_y(y_new, std_pos[1]);
        std::normal_distribution<double> dist_theta(theta_new, std_pos[2]);

        double sample_x = dist_x(gen);
        double sample_y = dist_y(gen);
        double sample_theta = dist_theta(gen);

        particles[i].x = sample_x;
        particles[i].y = sample_y;
        particles[i].theta = sample_theta;
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    for(auto particle: particles)
    {
        // Transformation of observations into map coordinates given the position of the particle
        vector<LandmarkObs> obs_transformed;

        for (auto observation:observations)
        {
            LandmarkObs obs;
            obs.x = particle.x + (cos(particle.theta) * observation.x) - (sin(particle.theta) * observation.y);
            obs.y = particle.y + (sin(particle.theta) * observation.x) + (cos(particle.theta) * observation.y);

            obs_transformed.push_back(obs);
        }

        // Association each transformed observation with a landmark identifier
        for (auto observation:obs_transformed) {
            double min_distance = std::numeric_limits<double>::max();
            for (auto landmark:map_landmarks.landmark_list)
            {
                if(dist(particle.x, particle.y, landmark.x_f, landmark.y_f) <= sensor_range)
                {
                    double distance = dist(observation.x, observation.y, landmark.x_f, landmark.y_f);
                    if(distance < min_distance){
                        observation.id = landmark.id_i;
                        min_distance = distance;
                    }
                }
            }
        }
    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}