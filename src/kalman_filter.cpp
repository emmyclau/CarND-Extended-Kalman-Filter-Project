#include "kalman_filter.h"
#include "tools.h"
#include <math.h>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Constructor
 */
KalmanFilter::KalmanFilter() {}

/**
 * Destructor
 */
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

/**
 * Prediction
 */
void KalmanFilter::Predict() {
  
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
    
}

/**
 * Laser update
 */
void KalmanFilter::Update(const VectorXd &z) {

    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
    
    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

/**
 * Radar update
 */
void KalmanFilter::UpdateEKF(const VectorXd &z) {
    
    //map the state vector from cartesian to polar coordinates
    float px = x_[0];
    float py = x_[1];
    float vx = x_[2];
    float vy = x_[3];
    
    VectorXd z_pred = VectorXd(3);
    
    float rho = sqrt(px*px + py*py);
    float phi = atan2(py, px);
    float rho_dot = (px*vx + py*vy)/sqrt(px*px + py*py);
    
    z_pred << rho, phi, rho_dot;
    VectorXd y = z - z_pred;
    
    //normalizing angles
    if (y(1) > M_PI) {
        y(1) = fmod(y(1), 2*M_PI);
    } else if (y(1) < -M_PI) {
        y(1) = fmod(y(1), -2*M_PI);
    }    
    
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
    
    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
    
}
