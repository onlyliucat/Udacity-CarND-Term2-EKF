#define _USE_MATH_DEFINES

#include "kalman_filter.h"
#include <math.h>
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
	x_ = F_*x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_*P_*Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
	VectorXd z_pred = H_*x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_*P_*Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_*Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K*y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size,x_size);
	P_ = (I - K*H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */

  // note that H_ was already set to Hj_ in FusionEKF.cpp
	// but we need to convert cartesion coordinates to polar


	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);
	float rho = sqrt(px*px + py*py);

	if(rho == 0) {
		std::cout << "MERP: Rho calculates to 0 in KalmanFilter::UpdateEKF" << std::endl;
		return;
	}
	double theta = atan2(py,px);
	float rhodot = px*vx+ py*vy;
	rhodot = rhodot / rho;
	VectorXd z_pred = VectorXd(3);
	z_pred << rho,theta,rhodot;
	VectorXd y = z - z_pred;

	double check = y(1);
	if(check > M_PI) {
		while(check > M_PI) {
			check -= 2*M_PI;
		}
	} else if(check < -1.0*M_PI) {
		while(check < -1.0*M_PI) {
			check += 2*M_PI;
		}
  }
  y(1) = check;

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_*P_*Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_*Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K*y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size,x_size);
	P_ = (I - K*H_) * P_;
}
