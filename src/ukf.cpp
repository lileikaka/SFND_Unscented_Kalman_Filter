#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	//std_a_ = 0.2;
	//std_a_ = 2.5;
	std_a_ = 5.0;

	// Process noise standard deviation yaw acceleration in rad/s^2
	//std_yawdd_ = 0.2;
	std_yawdd_ = 0.5;

	/**
	 * DO NOT MODIFY measurement noise values below.
	 * These are provided by the sensor manufacturer.
	 */

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	/**
	 * End DO NOT MODIFY section for measurement noise values 
	 */


	/**
	 * TODO: Complete the initialization. See ukf.h for other member properties.
	 * Hint: one or more values initialized above might be wildly off...
	 */

	// set state dimension
	n_x_ = 5;
	// set augmented dimension
	n_aug_ = 7;    

	n_z_ = 3;
	n_z_l = 2;
	lambda_ = 3 - n_aug_;

	// set weights
	weights = VectorXd(2 * n_aug_ + 1);
	weights(0) = lambda_ / (lambda_ + n_aug_);
	weights.tail(weights.rows() - 1) = VectorXd::Ones(weights.rows() - 1)*0.5 / (lambda_ + n_aug_);
	Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);


	Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
	z_pred = VectorXd(n_z_);
	S = MatrixXd(n_z_, n_z_);

	Zsig_l = MatrixXd(n_z_l, 2 * n_aug_ + 1);
	z_pred_l = VectorXd(n_z_l);
	S_l = MatrixXd(n_z_l, n_z_l);
	x_ << 0,
		0,
		0,
		0,
		0;
	P_ << 1.0, 0.000, 0.000, 0.000,0.000,
		0.000, 1.0, 0.000, 0.000,0.000,
		0.000, 0.000, 1.00, 0.000,0.000,
		0.000, 0.000, 0.000, 1.00,0.000,
		0.000, 0.000, 0.000, 0.000,1.00;
	is_initialized_=false;
	use_fmod = true;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

	if (!is_initialized_)
	{
		//cout << "Kalman Filter Initialization " << endl;

		// set the state with the initial location and zero velocity

		if(meas_package.sensor_type_ ==MeasurementPackage::LASER)
		{

			x_ << meas_package.raw_measurements_(0),
				meas_package.raw_measurements_(1),
				0,
				0,
				0;
			P_(0,0) = std_laspx_*std_laspx_;
			P_(1,1) = std_laspy_*std_laspy_;
			P_(2,2) = 0.01;
			P_(3,3) = 0.01;
			P_(4,4) = 0.01;
			
			time_us_ = meas_package.timestamp_;
			//is_initialized_ = true;
			return;			
		} 
		else
		{
			double rho = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			double rho_ = meas_package.raw_measurements_(2);
			double x = rho*cos(phi);
			double y = rho*sin(phi);
			double vx = rho_*cos(phi);
			double vy = rho_*sin(phi);
			double v = sqrt(vx*vx + vy*vy);
			//x_(0) = x;
			//x_(1) = y;
			//x_(2) = rho_;
			std::cout << "x = " << x << std::endl << "y = " << y << std::endl;
			std::cout << "rho_ = " << rho_ << std::endl << "v = " << v << std::endl;
			std::cout << "phi = " << phi << std::endl << "atan = " << atan2(x_(1),x_(0)) << std::endl;
			if(x_(0) < 0.0)
			{
				x_(2) = -rho_;
			}
			else
			{
				x_(2) = rho_;
			}
			//x_(2) = v;
			//x_(3) = phi;
			//P_(2,2) = std_radrd_*std_radrd_;
			//P_(3,3) = std_radphi_*std_radphi_;

			time_us_ = meas_package.timestamp_;				
			is_initialized_ = true;
			return;
		}

	}

	/*
  if(time_us_ != meas_package.timestamp_)
  {
	AugmentedSigmaPoints(&Xsig,x_,P_);
	SigmaPointPrediction(&Xsig_pred,Xsig,meas_package.timestamp_-time_us_);
	PredictMeanAndCovariance(&x_,&P_,Xsig_pred);	
	PredictLidarMeasurement(&z_pred_l,&Zsig_l,&S_l,Xsig_pred);
	time_us_ = meas_package.timestamp_;
  }
  */
  

  if(meas_package.sensor_type_ ==MeasurementPackage::LASER)
  {
	/*
	++cnt_lidar;
	if(cnt_lidar<3)
	{
		double diffx = meas_package.raw_measurements_[0] - x_(0);
		double diffy = meas_package.raw_measurements_[1] - x_(1);
		x_ << meas_package.raw_measurements_[0],
			meas_package.raw_measurements_[1],
			sqrt(diffx*diffx+diffy*diffy)*(meas_package.timestamp_-time_us_)/1000000,
			0,
      		0;
	}
	*/
	///*
	AugmentedSigmaPoints(&Xsig,x_,P_);
	SigmaPointPrediction(&Xsig_pred,Xsig,meas_package.timestamp_-time_us_);
	PredictMeanAndCovariance(&x_,&P_,Xsig_pred);	
	PredictLidarMeasurement(&z_pred_l,&Zsig_l,&S_l,Xsig_pred);
	time_us_ = meas_package.timestamp_;
    UpdateLidar(meas_package);
	//*/
  }
  else
  {
	///*
	AugmentedSigmaPoints(&Xsig,x_,P_);
	//SigmaPointPrediction(&Xsig_pred,Xsig,meas_package.timestamp_-time_us_);
	//PredictMeanAndCovariance(&x_,&P_,Xsig_pred);	
	//PredictRadarMeasurement(&z_pred,&Zsig,&S,Xsig_pred);
	PredictRadarMeasurement(&z_pred,&Zsig,&S,Xsig);
	time_us_ = meas_package.timestamp_;
    UpdateRadar(meas_package);
	//*/
  }  
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
	VectorXd z = VectorXd(n_z_l);
	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
	UpdateLidarState(&x_,&P_,Xsig_pred,x_,P_,Zsig_l,z_pred_l,S_l,z);
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  VectorXd z = VectorXd(n_z_);
  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);
  UpdateRadarState(&x_,&P_,Xsig_pred,x_,P_,Zsig,z_pred,S,z);

}

void UKF::UpdateLidarState(VectorXd* x_out, MatrixXd* P_out,MatrixXd& Xsig_pred,VectorXd& x,MatrixXd& P,MatrixXd& Zsig,VectorXd& z_pred,MatrixXd& S,VectorXd& z) {

  // create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_l);

	/**
	* Student part begin
	*/
	///*
	// calculate cross correlation matrix
	MatrixXd diff_x = MatrixXd(n_x_, 2 * n_aug_ + 1);
	MatrixXd weighted_diff_x = MatrixXd(n_x_, 2 * n_aug_ + 1);
	MatrixXd diff_z = MatrixXd(n_z_l, 2 * n_aug_ + 1);
	for (int i = 0;i < weighted_diff_x.cols();i++)
	{
		diff_x.col(i) = Xsig_pred.col(i) - x;
		
		if(use_fmod)
		{
			diff_x.col(i)(3) = fmod(diff_x.col(i)(3),2.*M_PI);
		}
		else
		{
			while (diff_x.col(i)(3)> M_PI) diff_x.col(i)(3) -= 2.*M_PI;
			while (diff_x.col(i)(3)<-M_PI) diff_x.col(i)(3) += 2.*M_PI;
		}
		
		weighted_diff_x.col(i) = weights(i)*diff_x.col(i);
		diff_z.col(i) = Zsig.col(i) - z_pred;
		//while (diff_z.col(i)(1)> M_PI) diff_z.col(i)(1) -= 2.*M_PI;
		//while (diff_z.col(i)(1)<-M_PI) diff_z.col(i)(1) += 2.*M_PI;
		
	}
	Tc = weighted_diff_x * diff_z.transpose();
	//std::cout << "TC: " << std::endl << Tc << std::endl;
	// calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();
	// update state mean and covariance matrix
	x = x + K * (z - z_pred);
	P = P - K * S * K.transpose();
	//*/
	/**
	* Student part end
	*/

	// print result
	//std::cout << "Updated state x: " << std::endl << x << std::endl;
	//std::cout << "Updated state covariance P: " << std::endl << P << std::endl;

	// write result
	*x_out = x;
	*P_out = P;
}

void UKF::UpdateRadarState(VectorXd* x_out, MatrixXd* P_out,MatrixXd& Xsig_pred,VectorXd& x,MatrixXd& P,MatrixXd& Zsig,VectorXd& z_pred,MatrixXd& S,VectorXd& z) {

  // create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z_);

	/**
	* Student part begin
	*/
	///*
	// calculate cross correlation matrix
	MatrixXd diff_x = MatrixXd(n_x_, 2 * n_aug_ + 1);
	MatrixXd weighted_diff_x = MatrixXd(n_x_, 2 * n_aug_ + 1);
	MatrixXd diff_z = MatrixXd(n_z_, 2 * n_aug_ + 1);
	for (int i = 0;i < weighted_diff_x.cols();i++)
	{
		diff_x.col(i) = Xsig_pred.col(i) - x;
		
		if(use_fmod)
		{
			diff_x.col(i)(3) = fmod(diff_x.col(i)(3),2.*M_PI);
		}
		else
		{
			while (diff_x.col(i)(3)> M_PI) diff_x.col(i)(3) -= 2.*M_PI;
			while (diff_x.col(i)(3)<-M_PI) diff_x.col(i)(3) += 2.*M_PI;
		}
		
		weighted_diff_x.col(i) = weights(i)*diff_x.col(i);
		diff_z.col(i) = Zsig.col(i) - z_pred;
		
		if(use_fmod)
		{
			diff_z.col(i)(1) = fmod(diff_z.col(i)(1),2.*M_PI);
		}
		else
		{
			while (diff_z.col(i)(1)> M_PI) diff_z.col(i)(1) -= 2.*M_PI;
			while (diff_z.col(i)(1)<-M_PI) diff_z.col(i)(1) += 2.*M_PI;
		}
		
	}
	Tc = weighted_diff_x * diff_z.transpose();
	//std::cout << "TC: " << std::endl << Tc << std::endl;
	// calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();
	// update state mean and covariance matrix
	x = x + K * (z - z_pred);
	P = P - K * S * K.transpose();
	//*/
	/**
	* Student part end
	*/

	// print result
	//std::cout << "Updated state x: " << std::endl << x << std::endl;
	//std::cout << "Updated state covariance P: " << std::endl << P << std::endl;

	// write result
	*x_out = x;
	*P_out = P;
}

void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* z_sig,MatrixXd* S_out,MatrixXd& Xsig_pred) {


	// create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z_l, 2 * n_aug_ + 1);

	// mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_l);

	// measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_l, n_z_l);

	/**
	* Student part begin
	*/

	Zsig = Xsig_pred.topLeftCorner(n_z_l,2 * n_aug_ + 1);
	// calculate mean predicted measurement
	z_pred =Zsig*weights;
	// calculate innovation covariance matrix S
	///*
	MatrixXd z_diff = MatrixXd(n_z_l, 2 * n_aug_ + 1);
	MatrixXd Buffer = MatrixXd(n_z_l, 2 * n_aug_ + 1);	
	
	for (int i = 0;i < z_diff.cols();i++)
	{
		z_diff.col(i) = Zsig.col(i) - z_pred;
		Buffer.col(i) = z_diff.col(i) * weights(i);
	}
	MatrixXd R = MatrixXd::Zero(n_z_l, n_z_l);
	R(0, 0) = std_radr_ * std_radr_;
	R(1, 1) = std_radphi_ * std_radphi_;	

	S = Buffer * z_diff.transpose() + R;
	//*/
	/**
	* Student part end
	*/

	// print result
	//std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	//std::cout << "S: " << std::endl << S << std::endl;

	// write result
	*z_out = z_pred;
	*S_out = S;
	*z_sig = Zsig;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* z_sig,MatrixXd* S_out,MatrixXd& Xsig_pred) {


	// create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

	// mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);

	// measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_, n_z_);

	/**
	* Student part begin
	*/

	// transform sigma points into measurement space
	for (int i = 0; i < Xsig_pred.cols();i++)
	{
		double px = Xsig_pred.col(i)(0);
		double py = Xsig_pred.col(i)(1);
		double v = Xsig_pred.col(i)(2);
		double yaw = Xsig_pred.col(i)(3);
		Zsig.col(i)(0) = sqrt(px*px + py * py);
		Zsig.col(i)(1) = atan2(py, px);
		Zsig.col(i)(2) = (px*cos(yaw)*v + py * sin(yaw)*v) / Zsig.col(i)(0);
	}	
	// calculate mean predicted measurement
	z_pred =Zsig*weights;
	// calculate innovation covariance matrix S
	///*
	MatrixXd z_diff = MatrixXd(n_z_, 2 * n_aug_ + 1);
	MatrixXd Buffer = MatrixXd(n_z_, 2 * n_aug_ + 1);	
	
	for (int i = 0;i < z_diff.cols();i++)
	{
		z_diff.col(i) = Zsig.col(i) - z_pred;
		
		if(use_fmod)
		{
			z_diff.col(i)(1) = fmod(z_diff.col(i)(1),2.*M_PI);
		}
		else
		{
			while (z_diff.col(i)(1)> M_PI) z_diff.col(i)(1) -= 2.*M_PI;
			while (z_diff.col(i)(1)<-M_PI) z_diff.col(i)(1) += 2.*M_PI;
		}		
		
		Buffer.col(i) = z_diff.col(i) * weights(i);
	}
	MatrixXd R = MatrixXd::Zero(n_z_, n_z_);
	R(0, 0) = std_radr_ * std_radr_;
	R(1, 1) = std_radphi_ * std_radphi_;
	R(2, 2) = std_radrd_ * std_radrd_;

	S = Buffer * z_diff.transpose() + R;
	//*/
	/**
	* Student part end
	*/

	// print result
	//std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	//std::cout << "S: " << std::endl << S << std::endl;

	// write result
	*z_out = z_pred;
	*S_out = S;
	*z_sig = Zsig;
}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out, MatrixXd& Xsig_pred) {
	
	// create vector for predicted state
	VectorXd x = VectorXd(n_x_);
	
	// create covariance matrix for prediction
	MatrixXd P = MatrixXd(n_x_, n_x_);

	/**
	* Student part begin
	*/
	// predict state mean
	x = Xsig_pred * weights;
	//std::cout << "Predicted multiply weights state" << std::endl;
	//std::cout << x << std::endl;
	// predict state covariance matrix
	/* my solution	
	MatrixXd x_diff = MatrixXd(n_x_, 2 * n_aug_ + 1);
	MatrixXd Buffer = MatrixXd(n_x_, 2 * n_aug_ + 1);
	for (int i = 0;i < x_diff.cols();i++)
	{
		x_diff.col(i) = Xsig_pred.col(i)-x;
		while (x_diff.col(i)(3)> M_PI) x_diff.col(i)(3)-=2.*M_PI;
		while (x_diff.col(i)(3)<-M_PI) x_diff.col(i)(3)+=2.*M_PI;
		Buffer.col(i) = x_diff.col(i) * weights(i);
	}
	P = Buffer * x_diff.transpose();
	*/

	///* udacity solution
	P.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; ++i)// iterate over sigma points
	{  
		// state difference
		VectorXd x_diff = Xsig_pred.col(i) - x;
		// angle normalization
		if(use_fmod)
		{
			x_diff(3) = fmod(x_diff(3),2.*M_PI);
		}
		else
		{
			while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;		
			while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		}

		P = P + weights(i) * x_diff * x_diff.transpose() ;
	}
	//*/
	/**
	* Student part end
	*/

	// print result
  /*
	std::cout << "Predicted state" << std::endl;
	std::cout << x << std::endl;
	std::cout << "Predicted covariance matrix" << std::endl;
	std::cout << P << std::endl;
  */
	// write result
	*x_out = x;
	*P_out = P;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out,MatrixXd& Xsig_aug,double delta_t) {



	// create matrix with predicted sigma points as columns
	MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

	//double delta_t = 0.1; // time diff in sec
	delta_t = delta_t /1000000;
	/**
	* Student part begin
	*/

	// predict sigma points
	double sqr_delta_t = delta_t * delta_t;
	// avoid division by zero
	for (int i = 0; i< Xsig_aug.cols(); i++)
	{
		double px = Xsig_aug.col(i)(0);
		double py = Xsig_aug.col(i)(1);
		double v = Xsig_aug.col(i)(2);
		double yaw = Xsig_aug.col(i)(3);
		double yaw_ = Xsig_aug.col(i)(4);
		double va = Xsig_aug.col(i)(5);
		double v_yaw_ = Xsig_aug.col(i)(6);
		if (fabs(yaw_) < 0.001)
		{
			px = px + v * cos(yaw)*delta_t + 0.5*sqr_delta_t*cos(yaw)*va;
			py = py + v * sin(yaw)*delta_t + 0.5*sqr_delta_t*sin(yaw)*va;
		}
		else
		{
			px = px + v / yaw_ * (sin(yaw + delta_t*yaw_) - sin(yaw)) + 0.5*sqr_delta_t*cos(yaw)*va;
			py = py + v / yaw_ * (-cos(yaw + delta_t * yaw_) + cos(yaw)) + 0.5*sqr_delta_t*sin(yaw)*va;			
		}
		v = v + delta_t * va;
		yaw = yaw + yaw_ * delta_t + 0.5*sqr_delta_t*v_yaw_;
		yaw_ = yaw_ + delta_t * v_yaw_;

		Xsig_pred.col(i)(0) = px;
		Xsig_pred.col(i)(1) = py;
		Xsig_pred.col(i)(2) = v;
		Xsig_pred.col(i)(3) = yaw;
		Xsig_pred.col(i)(4) = yaw_;
		//Xsig_pred.col(i)(5) = va;
		//Xsig_pred.col(i)(6) = v_xi_;		
	}
	// write predicted sigma points into right column

	/**
	* Student part end
	*/

	// print result
	//std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

	// write result
	*Xsig_out = Xsig_pred;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out,VectorXd& x,MatrixXd& P) {

	// create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);

	// create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	// create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	/**
	* Student part begin
	*/

	// create augmented mean state
	x_aug.head(n_x_) = x;
	x_aug(5) = 0;
  	x_aug(6) = 0;
	//x_aug(n_x) = std_a;
	//x_aug(n_x+1) = std_yawdd;
	// create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_,n_x_) = P;
	P_aug(5,5) = std_a_*std_a_;
	P_aug(6,6) = std_yawdd_*std_yawdd_;

	// create square root matrix
	// create square root matrix
	///*
	MatrixXd L = P_aug.llt().matrixL();

	// create augmented sigma points
	Xsig_aug.col(0)  = x_aug;
	for (int i = 0; i< n_aug_; ++i) {
		Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
	}
	//*/
	/*
	MatrixXd A = P_aug.llt().matrixL();
	double sqrt_lamda_x = sqrt(lambda_ + n_aug_);
	// create augmented sigma points
	for (int i = 0;i<2 * n_aug_ + 1;i++)
	{
		Xsig_aug.col(i) = x_aug;
	}
	Xsig_aug.block(0, 1, n_aug_, n_aug_) = Xsig_aug.block(0, 1, n_aug_, n_aug_) + sqrt_lamda_x * A;
	Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) = Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) - sqrt_lamda_x * A;
	*/
	/**
	* Student part end
	*/

	// print result
	//std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

	// write result
	*Xsig_out = Xsig_aug;
}