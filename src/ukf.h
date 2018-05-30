#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {

public:

  ///* 1. initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* 2. if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* 3. if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* 4. state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* 5. state covariance matrix
  MatrixXd P_;

  ///* 6. predicted sigma points matrix
  MatrixXd Xsig_pred;

  ///* 7. time when the state is true, in us
  long long time_us_;

  ///* 8. Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* 9. Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* 10a. Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* 10b. Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* 10c. Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* 10d. Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* 10e. Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* 11. Weights of sigma points
  VectorXd weights;

  ///* 12. State dimension
  int n_x;

  ///* 13. Augmented state dimension
  int n_aug;

  ///* 14. Sigma point spreading parameter
  double lambda_;

  // parts that are interesting
  VectorXd x_aug; // 7x1
  MatrixXd P_aug; // 7x7
  MatrixXd Xsig_aug; //7x15
  VectorXd x_pred; //5x1
  MatrixXd P_pred;  //5x5
  int n_z; // radar: 3, lidar: 2
//  MatrixXd R;
  
    
    
    
    
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
//  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

};

#endif /* UKF_H */
