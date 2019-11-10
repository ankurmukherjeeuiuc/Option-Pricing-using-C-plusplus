// pricing an european call option using fast fourier transform (fft)
// the technical details can be found in p. carr and d. madan, 
// "option valuation using the fast fourier transform," journal of 
// computational finance, volume 2, number 4, summer, 1999.
// this paper is available on the compass website.  it require some
// mathematical heavy-lifting.  come see me if you have questions,
// as the nitty-gritty of the paper are a little beyond the scope of
// the class. the method, however, is not...
// 
// the grid size must be of the form 2^k for some integer k for the 
// fft algorithm to work efficiently.   you will need to link with the
// newmat library for the fft algorithm etc.
// 

// i have only done the evaluation of a put option... i leave the call
// option pricing as a (challenging) exercise to you!

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include "newmat.h"
#include "newmatap.h"  // need this for the fft algorithm
#include "normdist.h"  // need this for odegaard's black-scholes routine
#define pi 3.141592654 // declaring the constant pi for future use

using namespace std;

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility;
long grid_size;

double max(double a, double b) {
	return (b < a) ? a : b;
}
// prof. odegaard's implementation of the black-scholes call price formula
double option_price_call_black_scholes(const double& s,       // spot (underlying) price
	const double& k,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time)	  // time to maturity 
{
	double time_sqrt = sqrt(time);
	double d1 = (log(s / k) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	double d2 = d1 - (sigma * time_sqrt);
	return s * n(d1) - k * exp(-r * time) * n(d2);
};

// odegaard's black-scholes put c++ code
double option_price_put_black_scholes(const double& s,      // spot price
	const double& k,      // strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(s / k) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	double d2 = d1 - (sigma * time_sqrt);
	return k * exp(-r * time) * n(-d2) - s * n(-d1);
};


double n(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-z) * (z / 2.0));
	double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

// my implemenation of the carr-madan approach
double call_price_from_fft()
{
	double log_of_initial_stock_price = log(initial_stock_price);
	double log_of_strike_price = log(strike_price);

	// equation 20 of carr & madan's paper
	double b = log_of_strike_price;
	double lambda = (2 * b) / ((double)grid_size);

	// equation 23 of carr & madan's paper
	double eta = (2 * pi) / (lambda * ((double)grid_size));

	// we are going to use an alpha of 2
	double alpha = 5.0;

	// characteristic function 6 and the expression for phi_t in my notes
	// for this we have to create a grid and do the necessary incantations
	//complex <double> final_series[grid_size];

	complex <double>* final_series = new complex<double>[grid_size];

	for (int i = 0; i < grid_size; i++) {

		complex <double> numerator;

		complex <double> phi(((double)i) * eta, -(alpha + 1));

		// simpson's rule weights
		double weight = 3.0 + pow(-1.0, (double)i + 1);
		if (i == 0)
			weight -= 1.0;
		weight = weight / 3.0;

		complex <double> temp1((volatility * volatility * expiration_time) / 2, 0);
		temp1 = temp1 * phi * phi;

		complex <double> temp2(0.0, log_of_initial_stock_price + (risk_free_rate - (volatility * volatility * 0.5)) * expiration_time);
		temp2 = temp2 * phi;

		temp1 = exp(temp2 - temp1);

		numerator = weight * exp(-risk_free_rate * expiration_time) * temp1;


		complex <double> denominator((alpha * alpha) + alpha - (((double)i * i) * eta * eta), (2 * alpha + 1) * (((double)i) * eta));

		complex <double> temp3(0.0, log_of_strike_price * eta * ((double)i));
		temp3 = exp(temp3);

		final_series[i] = (eta * temp3 * numerator) / denominator;
	}

	// have to get the real and imaginary parts for fft
	ColumnVector real_part((int)grid_size), imaginary_part((int)grid_size);
	for (int i = 1; i <= grid_size; i++) {
		real_part(i) = final_series[i - 1].real();
		imaginary_part(i) = final_series[i - 1].imag();
	}
	ColumnVector fft_real_part((int)grid_size), fft_imaginary_part((int)grid_size);
	FFT(real_part, imaginary_part, fft_real_part, fft_imaginary_part);
	return(fft_real_part(1) * exp(-alpha * log_of_strike_price) / pi);
}

// my implemenation of the carr-madan approach
double put_price_from_fft()
{
	double log_of_initial_stock_price = log(initial_stock_price);
	double log_of_strike_price = log(strike_price);

	// equation 20 of carr & madan's paper
	double b = log_of_strike_price;
	double lambda = (2 * b) / ((double)grid_size);

	// equation 23 of carr & madan's paper
	double eta = (2 * pi) / (lambda * ((double)grid_size));

	// we are going to use an alpha of 2
	double alpha = -5.0;

	// characteristic function 6 and the expression for phi_t in my notes
	// for this we have to create a grid and do the necessary incantations
	complex <double>* final_series = new complex<double>[grid_size];

	for (int i = 0; i < grid_size; i++) {

		complex <double> numerator;

		complex <double> phi(((double)i) * eta, -(alpha + 1));

		// simpson's rule weights
		double weight = 3.0 + pow(-1.0, (double)i + 1);
		if (i == 0)
			weight -= 1.0;
		weight = weight / 3.0;

		complex <double> temp1((volatility * volatility * expiration_time) / 2, 0.0);
		temp1 = temp1 * phi * phi;

		complex <double> temp2(0.0, log_of_initial_stock_price + (risk_free_rate - (volatility * volatility * 0.5)) * expiration_time);
		temp2 = temp2 * phi;

		temp1 = exp(temp2 - temp1);

		numerator = weight * exp(-risk_free_rate * expiration_time) * temp1;

		complex <double> denominator((alpha * alpha) + alpha - (((double)i * i) * eta * eta), (2 * alpha + 1) * (((double)i) * eta));

		complex <double> temp3(0.0, log_of_strike_price * eta * ((double)i));
		temp3 = exp(temp3);

		final_series[i] = (eta * temp3 * numerator) / denominator;
	}

	// have to get the real and imaginary parts for fft
	ColumnVector real_part((int)grid_size), imaginary_part((int)grid_size);
	for (int i = 1; i <= grid_size; i++) {
		real_part(i) = final_series[i - 1].real();
		imaginary_part(i) = final_series[i - 1].imag();
	}
	ColumnVector fft_real_part((int)grid_size), fft_imaginary_part((int)grid_size);
	FFT(real_part, imaginary_part, fft_real_part, fft_imaginary_part);
	return(fft_real_part(1) * exp(-alpha * log_of_strike_price) / pi);
}


int main(int argc, char* argv[])
{
	/*sscanf(argv[1], "%lf", &expiration_time);
	sscanf(argv[2], "%ld", &grid_size);
	sscanf(argv[3], "%lf", &risk_free_rate);
	sscanf(argv[4], "%lf", &volatility);
	sscanf(argv[5], "%lf", &initial_stock_price);
	sscanf(argv[6], "%lf", &strike_price);*/

	cout << "enter t,grid sizerf,volatility(in decimal), s0, k: " << endl;
	cin >> expiration_time >> grid_size>> risk_free_rate >> volatility >> initial_stock_price >> strike_price;

	// checking if the grid size is 2^k for some k
	long temp = grid_size;
	while ((temp >= 2) && (temp % 2 == 0))
		temp = temp / 2;

	if (temp == 1) {
		cout << "european put option pricing by fft" << endl;
		cout << "expiration time (years) = " << expiration_time << endl;
		cout << "grid size = " << grid_size << endl;
		cout << "risk free interest rate = " << risk_free_rate << endl;
		cout << "volatility (%age of stock value) = " << volatility * 100 << endl;
		cout << "initial stock price = " << initial_stock_price << endl;
		cout << "strike price = " << strike_price << endl;
		cout << "--------------------------------------" << endl;
		cout << "price of an european call option via carr & madan's fft method = " <<
			call_price_from_fft() << endl;
		cout << "price of an european call from the black-scholes formula = " <<
			option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time) << endl;
		cout << "--------------------------------------" << endl;
		cout << "price of an european put option via carr & madan's fft method = " <<
			put_price_from_fft() << endl;
		cout << "price of an european put from the black-scholes formula = " <<
			option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time) << endl;
		cout << "--------------------------------------" << endl;
	}
	else {
		cout << "your grid size is " << grid_size << endl;
		cout << "retry with something that is 2^k for some k" << endl;
		cout << "exiting..." << endl;
	}

}
