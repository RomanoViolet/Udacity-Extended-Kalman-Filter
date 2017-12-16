#include "kalman_filter.h"
#include <cmath>
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {
}

KalmanFilter::~KalmanFilter() {
}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
		MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}


double KalmanFilter::NormalizeYawAngle(double angle) {


/*	while (angle>  M_PI) angle-=(2.0 * M_PI);
	while (angle< -M_PI) angle+= (2.0 * M_PI);
	return (angle);
*/
	if (std::abs(angle) / (M_PI) > 1) {
		//https://stackoverflow.com/questions/24234609/standard-way-to-normalize-an-angle-to-%CF%80-radians-in-java
		return (std::atan2(std::sin(angle), std::cos(angle)));
	} else {
		return (angle);
	}

}/* KalmanFilter::NormalizeYawAngle */

void KalmanFilter::Predict() {

	/**
	 TODO:
	 * predict the state
	 */

	// The prediction uses linear model, therefore, linearization (e.g., via Taylor expansion) is not required.
	std::cout << "file: " << __FILE__ << std::endl;
	std::cout << "priori: x_: " << x_ << std::endl;
	std::cout << "priori: F_: " << F_ << std::endl;
	x_ = F_ * x_;
	std::cout << "post: x_: " << x_ << std::endl;
	MatrixXd Ft = F_.transpose();
	std::cout << "priori: P_: " << P_ << std::endl;
	P_ = F_ * P_ * Ft + Q_;
	std::cout << "post: P_: " << P_ << std::endl;
	std::cout << "post: Q_: " << Q_ << std::endl;


}

void KalmanFilter::Update(const VectorXd &z) {
	/**
	 TODO:
	 * update the state by using Kalman Filter equations
	 */

	// Compute y = z - Hx', where x' is the predicted state at the current timestep.
	VectorXd prediction = VectorXd(2);

	//H_ << 	1, 0, 0, 0,
	//		0, 1, 0, 0;

	// measurement matrix H has already been instantiated and initialized.
	prediction = H_ * x_;

	// now the error vector y
	VectorXd y = z - prediction;

	// Matrix S, which will be used to compute Kalman gain, K
	MatrixXd S = H_ * P_ * H_.transpose() + R_;

	// Kalman gain
	MatrixXd K = P_ * H_.transpose() * S.inverse();

	// new estimate: posterori
	x_ = x_ + (K * y);

	// precision update
	unsigned x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	 TODO:
	 * update the state by using Extended Kalman Filter equations
	 */
	// Compute y = z - Hx', where x' is the predicted state at the current timestep.
	// First compute Hx'
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	// inverse tangent value is computed using std::atan2 since the returned value is required to be between -pi and pi.
	VectorXd prediction = VectorXd(3);
	prediction << 	std::pow(std::pow(px, 2.0) + std::pow(py, 2.0), 0.5),
					std::atan2(py, px),
					(px * vx + py * vy) / (std::pow(std::pow(px, 2.0) + std::pow(py, 2.0), 0.5));

	std::cout << "prediction: rho" << prediction(0) << std::endl;
	std::cout << "prediction: phi" << prediction(1) << std::endl;
	std::cout << "prediction: rhodot" << prediction(2) << std::endl;


	std::cout << "measured: rho" << z(0) << std::endl;
	std::cout << "measured: phi" << z(1) << std::endl;
	std::cout << "measured: rhodot" << z(2) << std::endl;

	// now the error vector y
	VectorXd y = z - prediction;

	std::cout << "error: rho" << y(0) << std::endl;
	std::cout << "error: phi" << y(1) << std::endl;
	std::cout << "error: rhodot" << y(2) << std::endl;


	// correct the phi error, if necessary
/*	if(y(1)/(M_PI)>1)
	{
		y(1) = y(1)-(2*M_PI);
		std::cout << "corrected error: phi" << y(1) << std::endl;

	}*/

	y(1) = NormalizeYawAngle(y(1));

	// Compute the Jacobian
	MatrixXd Hj = tool.CalculateJacobian(x_);

	// Matrix S, which will be used to compute Kalman gain, K
	MatrixXd S = Hj * P_ * Hj.transpose() + R_;

	// Kalman gain
	MatrixXd K = P_ * Hj.transpose() * S.inverse();

	// new estimate: posterori
	x_ = x_ + (K * y);

	// precision update
	unsigned x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj) * P_;

}


