// calculating the price of an european option using dynamic 
// programming

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "newmat.h"
#include "normdist.h"          // this defines the normal distribution from odegaard's files
using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, r;
int no_of_divisions;

double max(double a, double b) {
	return (b < a) ? a : b;
}

matrix repeated_squaring(matrix a, int exponent, int no_rows)
{
	if (exponent == 0) {
		identitymatrix i(no_rows);
		return (i);
	}
	else if (exponent % 2 == 0)
		return (repeated_squaring(a * a, exponent / 2, no_rows));
	else
		return (a * repeated_squaring(a * a, (exponent - 1) / 2, no_rows));
}

double option_price_put_black_scholes(const double& s,      // spot price
	const double& k,      // strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time)
{
	double time_sqrt = sqrt(time);
	double d1 = (log(s / k) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	double d2 = d1 - (sigma * time_sqrt);
	return k * exp(-r * time) * n(-d2) - s * n(-d1);
};

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

double european_call_option_dyn_prog()
{
	// create the probability matrix
	matrix transition_probability(2 * no_of_divisions - 1, 2 * no_of_divisions - 1);
	matrix final_matrix(2 * no_of_divisions - 1, 2 * no_of_divisions - 1);

	for (int i = 1; i <= 2 * no_of_divisions - 1; i++)
		for (int j = 1; j <= 2 * no_of_divisions - 1; j++)
			transition_probability(i, j) = 0.0;

	// boundary values of the probabilities need to be entered
	transition_probability(1, 2) = 1.0;
	transition_probability(2 * no_of_divisions - 1, 2 * no_of_divisions - 1) = 1.0;

	for (int i = 2; i <= 2 * no_of_divisions - 2; i++)
		for (int j = 1; j <= 2 * no_of_divisions - 1; j++) {
			if (j == (i - 1))
				transition_probability(i, j) = 1 - uptick_prob;
			if (j == (i + 1))
				transition_probability(i, j) = uptick_prob;
		}

	// compute pi^no_of_divisions by method of repeated squaring
	// keep in mind that we are computing pi^(t-1) not pi^t 
	// because the pi matrix in itself represents the 1-step state
	// transition probability... for k-step you need to multiply pi
	// with itself (k-1) times 
	final_matrix = repeated_squaring(transition_probability, no_of_divisions - 1, 2 * no_of_divisions - 1);

	// value at expiration
	matrix v_t(2 * no_of_divisions - 1, 1);
	for (int i = 1; i <= 2 * no_of_divisions - 1; i++)
		v_t(i, 1) = (1.0 / pow(r, no_of_divisions)) *
		max(0.0, (initial_stock_price * pow(up_factor, i - no_of_divisions)) - strike_price);

	double final_result = 0;
	for (int i = 1; i < 2 * no_of_divisions - 1; i++)
		final_result += final_matrix(no_of_divisions, i) * v_t(i, 1);

	return (final_result);
}

double european_put_option_dyn_prog()
{
	// create the probability matrix
	matrix transition_probability(2 * no_of_divisions - 1, 2 * no_of_divisions - 1);
	matrix final_matrix(2 * no_of_divisions - 1, 2 * no_of_divisions - 1);

	for (int i = 1; i <= 2 * no_of_divisions - 1; i++)
		for (int j = 1; j <= 2 * no_of_divisions - 1; j++)
			transition_probability(i, j) = 0.0;

	// boundary values of the probabilities need to be entered
	transition_probability(1, 2) = 1.0;
	transition_probability(2 * no_of_divisions - 1, 2 * no_of_divisions - 1) = 1.0;

	for (int i = 2; i <= 2 * no_of_divisions - 2; i++)
		for (int j = 1; j <= 2 * no_of_divisions - 1; j++) {
			if (j == (i - 1))
				transition_probability(i, j) = 1 - uptick_prob;
			if (j == (i + 1))
				transition_probability(i, j) = uptick_prob;
		}

	// compute pi^no_of_divisions by method of repeated squaring
	// keep in mind that we are computing pi^(t-1) not pi^t 
	// because the pi matrix in itself represents the 1-step state
	// transition probability... for k-step you need to multiply pi
	// with itself (k-1) times 
	final_matrix = repeated_squaring(transition_probability, no_of_divisions - 1, 2 * no_of_divisions - 1);

	// value at expiration
	matrix v_t(2 * no_of_divisions - 1, 1);
	for (int i = 1; i <= 2 * no_of_divisions - 1; i++)
		v_t(i, 1) = (1.0 / pow(r, no_of_divisions)) *
		max(0.0, strike_price - (initial_stock_price * pow(up_factor, i - no_of_divisions)));

	double final_result = 0;
	for (int i = 1; i < 2 * no_of_divisions - 1; i++)
		final_result += final_matrix(no_of_divisions, i) * v_t(i, 1);

	return (final_result);
}

int main(int argc, char* argv[])
{

	cout << "enter t,n,rf,volatility(in decimal), s0, k: " << endl;
	cin >> expiration_time >> no_of_divisions >> risk_free_rate >> volatility >> initial_stock_price >> strike_price;

	up_factor = exp(volatility * sqrt(expiration_time / ((float)no_of_divisions)));
	r = exp(risk_free_rate * expiration_time / ((float)no_of_divisions));
	uptick_prob = (r - (1 / up_factor)) / (up_factor - (1 / up_factor));
	cout << "european put option pricing by dynamic programming" << endl;
	cout << "expiration time (years) = " << expiration_time << endl;
	cout << "number of divisions = " << no_of_divisions << endl;
	cout << "risk free interest rate = " << risk_free_rate << endl;
	cout << "volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "initial stock price = " << initial_stock_price << endl;
	cout << "strike price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "up factor = " << up_factor << endl;
	cout << "uptick probability = " << uptick_prob << endl;
	cout << "--------------------------------------" << endl;
	cout << "price of an european call option = " << european_call_option_dyn_prog() << endl;
	cout << "call price according to black-scholes = " <<
		option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	cout << "price of an european put option = " << european_put_option_dyn_prog() << endl;
	cout << "put price according to black-scholes = " <<
		option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time) << endl;
}

