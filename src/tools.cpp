#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {
}

Tools::~Tools() {
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth) {
	/**
	 TODO:
	 * Calculate the RMSE here.
	 */

	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// TODO: YOUR CODE HERE
	eigen_assert(
					(estimations.size() != 0)
				&& 	(estimations.size() == ground_truth.size())
				&& 	"size of vectors do not match."
				);


	//accumulate squared residuals
	VectorXd totalError(4);
	totalError << 0, 0, 0, 0;

	for (unsigned i = 0; i < estimations.size(); ++i) {
		VectorXd currentEstimate = estimations[i];
		VectorXd currentGroundTruth = ground_truth[i];
		VectorXd currentError = currentEstimate - currentGroundTruth;
		totalError = totalError + static_cast<VectorXd>(currentError.array().pow(2.0));
	}

	//calculate the mean
	VectorXd meanError = totalError.array() / estimations.size();

	//calculate the squared root
	rmse = meanError.array().pow(0.5);

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	/**
	 TODO:
	 * Calculate a Jacobian here.
	 */
	Eigen::MatrixXd Hj(3, 4);
		// recover state parameters
		float px = x_state(0);
		float py = x_state(1);
		float vx = x_state(2);
		float vy = x_state(3);

		// Mark divide by zero error
		assert((std::pow(px, 2) + std::pow(py, 2)) > 0.0001);

		Hj(0, 0) = px / std::pow(std::pow(px, 2) + std::pow(py, 2), 0.5);
		Hj(0, 1) = py / std::pow(std::pow(px, 2) + std::pow(py, 2), 0.5);
		Hj(0, 2) = 0;
		Hj(0, 3) = 0;

		Hj(1, 0) = -py / (std::pow(px, 2) + std::pow(py, 2));
		Hj(1, 1) = px / (std::pow(px, 2) + std::pow(py, 2));
		Hj(1, 2) = 0;
		Hj(1, 3) = 0;

		Hj(2, 0) = (py * ((vx * py) - (vy * px)))
				/ std::pow(std::pow(px, 2) + std::pow(py, 2), 3 / 2.0);
		Hj(2, 1) = (px * ((vy * px) - (vx * py)))
				/ std::pow(std::pow(px, 2) + std::pow(py, 2), 3 / 2.0);
		Hj(2, 2) = Hj(0, 0);
		Hj(2, 3) = Hj(0, 1);

		return Hj;
}
