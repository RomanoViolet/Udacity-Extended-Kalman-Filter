#include "FusionEKF.h"
#include "Eigen/Dense"
#include "tools.h"
#include <cmath>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
	is_initialized_ = false;

	previous_timestamp_ = 0;

	// initializing matrices
	R_laser_ = MatrixXd(2, 2);
	R_radar_ = MatrixXd(3, 3);
	H_laser_ = MatrixXd(2, 4);
	Hj_ = MatrixXd(3, 4);

	// measurement covariance matrix - laser
	R_laser_ << 0.0225, 0, 0, 0.0225;

	// measurement covariance matrix - radar
	R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

	// measurement noise
	noise_ax = 9.0;
	noise_ay = 9.0;

	// initial precision matrix
	P_ = MatrixXd(4, 4);
	P_ << 10000, 0, 0, 0,
		  0, 10000, 0, 0,
		  0, 0, 10000, 0,
		  0, 0, 0, 10000;

	/**
	 TODO:
	 * Finish initializing the FusionEKF.
	 * Set the process and measurement noises
	 */
}



/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

	/*****************************************************************************
	 *  Initialization
	 ****************************************************************************/
	if (!is_initialized_) {

		// set the state with the initial location and zero velocity
		ekf_.x_ = VectorXd(4);

		// update the timestamp
		previous_timestamp_ = measurement_pack.timestamp_;

		/**
		 TODO:
		 * Initialize the state ekf_.x_ with the first measurement.
		 * Create the covariance matrix.
		 * Remember: you'll need to convert radar from polar to cartesian
		 coordinates.
		 */
		// first measurement
		// cout << "EKF: " << endl;
		// ekf_.x_ = VectorXd(4);
		// ekf_.x_ << 1, 1, 1, 1;
		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			 Convert radar from polar to cartesian coordinates and initialize state.
			 */

			float rho = measurement_pack.raw_measurements_[0];
			float theta = measurement_pack.raw_measurements_[1];
			float dtrho = measurement_pack.raw_measurements_[2];

			// The inferred px. The supplied angle is in radians.
			ekf_.x_(0) = rho * cos(theta);

			// The inferred py. The supplied angle is in radians.
			ekf_.x_(1) = rho * sin(theta);

			// The inferred vx is set to zero
			ekf_.x_(2) = dtrho * cos(theta);

			// The inferred vy is set to zero
			ekf_.x_(3) = dtrho * sin(theta);

		} else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			/**
			 Initialize state.
			 */
			// Lidar provides us with px, and py directly.
			ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
		}

		// instantiate the prediction matrix F
		ekf_.F_ = MatrixXd::Identity(4, 4);

		//  instantiate the process covariance matrix Q
		ekf_.Q_ = MatrixXd(4, 4);

		// instantiate and initialize the measurement matrix for Laser case.
		ekf_.H_ = MatrixXd(2, 4);
		ekf_.H_ << 	1, 0, 0, 0,
				  	  	0, 1, 0, 0;

		// fill in the initial precision matrix
		ekf_.P_ = P_;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/

	/**
	 TODO:
	 * Update the state transition matrix F according to the new elapsed time.
	 - Time is measured in seconds.
	 * Update the process noise covariance matrix.
	 * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	 */

	// Find the time elapsed between predictions.
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

	// update matrix F
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

	// update the process covariance matrix Q
	ekf_.Q_ << 	(std::pow(dt, 4.0) * noise_ax) / 4.0, 									0, 				(std::pow(dt, 3)* noise_ax) / 2.0, 		0,
													0, (std::pow(dt, 4) * noise_ay) / 4.0, 												0, 		(std::pow(dt, 3) * noise_ay) / 2.0,
				(std::pow(dt, 3) * noise_ax) / 2.0, 									0, 				(std::pow(dt, 2) * noise_ax), 			0,
													0, (std::pow(dt, 3) * noise_ay) / 2.0, 												0, 		(std::pow(dt, 2) * noise_ay);

	ekf_.Predict();

	// update the "old" timestamp
	previous_timestamp_ = measurement_pack.timestamp_;




	/*****************************************************************************
	 *  Update
	 ****************************************************************************/

	/**
	 TODO:
	 * Use the sensor type to perform the update step.
	 * Update the state and covariance matrices.
	 */

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		//update R matrix for radar case
		ekf_.R_ = R_radar_;
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	} else {
		// Laser updates
		// update R matrix for laser case
		ekf_.R_ = R_laser_;
		ekf_.Update(measurement_pack.raw_measurements_);
	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
