#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

    // TODO: YOUR CODE HERE

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	int est_size = estimations.size();
	int gt_size  = ground_truth.size();

	if (est_size == 0 or gt_size == 0)
	{
		std::cout << "Error vectors should not be zero!" << std::endl;
		return rmse;
	}
	//  * the estimation vector size should equal ground truth vector size
	if (est_size != gt_size)
	{
		std::cout << "Vector sizes should be the same!" << std::endl;
		return rmse;
	}

	//accumulate squared residuals
	
	for(int i=0; i < estimations.size(); ++i){
    	VectorXd temp(4);
    	temp = (estimations[i].array() - ground_truth[i].array());
       	rmse = rmse.array() + (temp.array() * temp.array());
 	}

	//calculate the mean
	rmse = rmse.array()/est_size;

	//calculate the squared root
	rmse  = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE 
	float px2_py2 = px*px + py*py;
	//check division by zero
	if (fabs(px2_py2) < 0.0001)
	{
		// https://discussions.udacity.com/t/action-on-divide-by-zero-in-jacobian/229082
		Hj << 0, 0, 0, 0,
		      1e+9, 1e+9, 0, 0,
		      0, 0, 0, 0;

		//std::cout << "Error division by zero!" << std::endl;
		return Hj;
	}
	//compute the Jacobian matrix
	float sqrt_px2_py2 = sqrt(px2_py2);
	float pow_3_2_px2_py2 = px2_py2 * sqrt_px2_py2;// pow(px2_py2,3/2);
	

	Hj << (px / sqrt_px2_py2), (py / sqrt_px2_py2), 0, 0,
	      -(py/px2_py2), (px/px2_py2), 0, 0,
	      py*(vx*py - vy*px) / pow_3_2_px2_py2, px*(vy*px - vx*py) / pow_3_2_px2_py2, px / sqrt_px2_py2, py / sqrt_px2_py2;

	return Hj;
}
