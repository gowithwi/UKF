#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  
  // 1. is initialized?
  is_initialized_ = false;

  // 2. if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  
  // 3. if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // 4. initial state vector
  x_ = VectorXd(5);
    
  // 5.initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5,5);
  // 6. predicted sigma points matrix
  Xsig_pred = MatrixXd(5,15);
    
  // 7. time when the state is true, in us
  time_us_ = 0;

  // 8. Process noise standard deviation longitudinal acceleration in m/s^2
//  std_a_ = 0.8;
//  std_a_ = 4; // 0.075, 0.086, 0.369, 0.256 // 2,6 is worse,
  std_a_ = 3;

  // 9. Process noise standard deviation yaw acceleration in rad/s^2
//  std_yawdd_ = M_PI / 4; // 2,6 is worse
  std_yawdd_ = M_PI / 4;
    
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // 10a. Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // 10b. Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // 10c. Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // 10d. Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // 10e. Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // 13. Augmented state dimension
    n_aug = 7;
  
  // 11. Weights of sigma points
  weights = VectorXd(2*n_aug + 1);
    
  // 12. State dimension
  n_x = 5;
    
  // 14. Sigma point spreading parameter
//  lambda_ = 3 - n_aug;// -1 does not help, -2 make it worse, +1 is also worse
  lambda_ = 3 - n_aug;
    
  // parts that are interesting
  x_aug = VectorXd(n_aug);
  P_aug = MatrixXd(n_aug, n_aug);
  Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  x_pred = VectorXd (5);
  P_pred = MatrixXd(5,5);
    

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   Here we have three parts:
   1. initializatin;
   2. prediction;
   3. update.
  */
  // 1. initialization
    if (!is_initialized_) {
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

            float rho;
            float yaw;
            rho = meas_package.raw_measurements_[0];
            yaw = meas_package.raw_measurements_[1];
            x_ << rho*cos(yaw),rho*sin(yaw), 0, 0, 0;
            
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        }
        
        time_us_ = meas_package.timestamp_;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    
  // 2. generate sigma points-- class 7/18
    // 2.1 create augmented x_a (7x1), P_a (7x7), X_sig_a (7x15)


    // 2.2 input augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    
    // 2.3 input augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    
    // 2.4 create square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    
    // 2.5 input augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug) * L.col(i);
        Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda_+n_aug) * L.col(i);
    }
    
  // 3. predict sigma points-- 7/21
    const float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;    //dt - expressed in seconds
    time_us_ = meas_package.timestamp_;
    
    for (int i = 0; i< 2*n_aug+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }
    
    // 4. predict mean and covariance-- 7/24
    double weight_0 = lambda_/(lambda_+n_aug);
    weights(0) = weight_0;
    for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug+lambda_);
        weights(i) = weight;
    }
    
        //predicted state mean
    x_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
        x_pred = x_pred+ weights(i) * Xsig_pred.col(i);
    }
    
        //predicted state covariance matrix
    P_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
        
        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x_pred;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P_pred = P_pred + weights(i) * x_diff * x_diff.transpose() ;
    }
    
    // 5. predict measurement-- 7/27
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package); // update by radar
    }
    else {
        UpdateLidar(meas_package);// update by lidar
    }
    
}

//void UKF::Prediction(double delta_t) {
//
//    
//}


void UKF::UpdateLidar(MeasurementPackage meas_package) {
    VectorXd z = meas_package.raw_measurements_;
    
    n_z = 2; // 2 for lidar's px, py
    
    
    
    //innovation covariance matrix S
    MatrixXd H(2,5);
    H<< 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0;
    MatrixXd Ht = H.transpose();

    
    VectorXd z_pred = H*x_pred;
    VectorXd y = z - z_pred;
    
    MatrixXd R = MatrixXd(n_z,n_z);
    R << std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;
    
    MatrixXd S = H*P_pred*Ht + R;
    
    //Kalman gain K;
    MatrixXd K = P_pred*Ht*S.inverse();
        
    x_ = x_pred + K * y;
    P_ = P_pred - K*H*P_pred;
    
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    VectorXd z = meas_package.raw_measurements_;
    n_z = 3; // 3 for radar's rho, phi, rho dot
    MatrixXd Zsig(n_z, 2 * n_aug + 1);
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        double v  = Xsig_pred(2,i);
        double yaw = Xsig_pred(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        float eps  =0.001;
        float rho = sqrt(p_x*p_x+p_y*p_y);
        rho = std::max(eps, rho); // check rho is not equal to zero.
        
        float MyArctan;
        if(p_y==0&& p_x==0)
            MyArctan = 0;
        else
            MyArctan = atan2(p_y,p_x);
        
        Zsig(0,i) = rho;                        //rho
        Zsig(1,i) = MyArctan;                   //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / rho;   //rho_dot
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug+1; i++) {
        z_pred = z_pred + weights(i) * Zsig.col(i);
    }
    
    //innovation covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        S = S + weights(i) * z_diff * z_diff.transpose();
    }
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;
    S = S + R;
    
    // 6. update state-- 7/30 (radar only, just for now)
    MatrixXd Tc = MatrixXd(n_x, n_z);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x_pred;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    //update state mean and covariance matrix
    x_ = x_pred + K * z_diff;
    P_ = P_pred - K*S*K.transpose();
    
}
