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
    num_particles = 10;  // TODO: Set the number of particles

    std::random_device rd;
    std::default_random_engine gen(rd());

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

    std::random_device rd;
    std::default_random_engine gen(rd());

    for (int i = 0; i < num_particles; ++i) {

        double theta_0 = particles[i].theta;
        double x_0 = particles[i].x;
        double y_0 = particles[i].y;

        double x_new ;
        double y_new;
        double theta_new;

        if(fabs(yaw_rate) < 1e-4)
        {
             x_new = x_0 + velocity * delta_t * cos(theta_0);
             y_new = y_0 + velocity * delta_t * sin(theta_0);
             theta_new = theta_0;
        }
        else {
             x_new = x_0 + (velocity / yaw_rate) * (sin(theta_0 + yaw_rate * delta_t) - sin(theta_0));
             y_new = y_0 + (velocity / yaw_rate) * (cos(theta_0) - cos(theta_0 + yaw_rate * delta_t));
             theta_new = theta_0 + yaw_rate * delta_t;
        }

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
    double weight_norm{0.0};
    weights.clear();

    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];

    // calculate normalization term
    double gauss_norm;
    gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

    for(auto &particle:particles)
    {
        // Transformation of observations into map coordinates given the position of the particle
        vector<LandmarkObs> obs_transformed;

//        double particle_x = particles[i].x;
//        double particle_y = particles[i].y;
//        double particle_theta = particles[i].theta;

        for (auto observation:observations)
        {
            LandmarkObs obs;
            obs.id = 0;
            obs.x = particle.x + (cos(particle.theta) * observation.x) - (sin(particle.theta) * observation.y);
            obs.y = particle.y + (sin(particle.theta) * observation.x) + (cos(particle.theta) * observation.y);

            obs_transformed.push_back(obs);
        }

        // Association each transformed observation with a landmark identifier
        particle.associations.clear();
        particle.sense_x.clear();
        particle.sense_y.clear();

        for (auto &observation:obs_transformed) {
            double min_distance = std::numeric_limits<double>::max();
            for (auto &landmark:map_landmarks.landmark_list)
            {
                if( dist(particle.x, particle.y, landmark.x_f, landmark.y_f) <= sensor_range)
                {
                    double distance = dist(observation.x, observation.y, landmark.x_f, landmark.y_f);
                    if(distance < min_distance){
                        observation.id = landmark.id_i;
                        min_distance = distance;
                    }
                }
            }

            particle.associations.push_back(observation.id);
            particle.sense_x.push_back(observation.x);
            particle.sense_y.push_back(observation.y);

        }

        // Calculating the Particle's Final Weight
        for (int j = 0; j < particle.associations.size(); ++j) {

            int landmark_id = particle.associations[j] - 1;

            double obs_x = particle.sense_x[j];
            double  obs_y = particle.sense_y[j];

            double landmark_x = map_landmarks.landmark_list[landmark_id].x_f;
            double landmark_y = map_landmarks.landmark_list[landmark_id].y_f;

            // calculate exponent
            double exponent;
            double x_diff = obs_x - landmark_x;
            double y_diff = obs_y - landmark_y;

            exponent = (pow(x_diff, 2) / (2 * pow(sig_x, 2)))
                       + (pow(y_diff, 2) / (2 * pow(sig_y, 2)));

            // calculate weight using normalization terms and exponent
            double weight = gauss_norm * exp(-exponent);

            if (weight == 0) {
                particle.weight *= 1e-6;
            } else {
                particle.weight *= weight;
            }
        }
        weight_norm += particle.weight;
        weights.push_back(particle.weight);
    }

    for (auto &particle:particles) {
        particle.weight /= weight_norm;

    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    std::random_device rd;
    std::default_random_engine gen(rd());

    std::discrete_distribution<int> dd_index(weights.begin(), weights.end());

    // create a bag of temp particles to transfer resampled particles
    std::vector<Particle> resampled_particles;

    for(int i=0; i<num_particles; i++){
        int index = dd_index(gen);
        resampled_particles.push_back(particles[index]);
    }

    // replace the original particles with survived particles
    particles = std::move(resampled_particles);

//   // find the particles max weight
//    double weight_norm{0};
//
//    double max_weight = std::numeric_limits<double>::min();
//    for (auto &particle:particles) {
//        if(particle.weight > max_weight) {
//            max_weight = particle.weight;
//        }
//    }
//
//    // Seed with a real random value, if available
//    std::random_device r;
//    std::default_random_engine gen(r());
//
//    std::uniform_real_distribution<double> uni_rand(0, 1);
//
//    std::uniform_int_distribution<int> particle_index(0, num_particles - 1);
//    int index = particle_index(gen);
//
//    double beta = 0.0;
//
//    std::vector<Particle> new_particles;
//
//    double random_value{1};
//
//    for(int i = 0; i < num_particles; ++i){
//        random_value = uni_rand(gen);
//        beta += random_value * 2.0 * max_weight;
//        while (beta > particles[index].weight)
//        {
//            beta -= particles[index].weight;
//            index = (index + 1) % num_particles;
//        }
//        new_particles.push_back(particles[index]);
//        weight_norm += particles[index].weight;
//    }
//
//    for (auto &particle:new_particles) {
//        particle.weight /= weight_norm;
//    }
//    particles = new_particles;
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