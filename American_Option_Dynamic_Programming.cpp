// calculating the price of an american option using dynamic 
// programming where the probability Matrix is derived from the
// trinomial lattice.  

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "newmat.h"
using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, r;
int no_of_divisions;

double max(double a, double b) {
	return (b < a) ? a : b;
}

double american_call_option_dyn_prog()
{
	// create the probability Matrix
	Matrix transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions + 1);

	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			transition_probability(i, j) = 0.0;

	// boundary values of the probabilities need to be entered
	transition_probability(1, 1) = 1.0 - uptick_prob;
	transition_probability(1, 2) = uptick_prob;
	transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions + 1) = uptick_prob;
	transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions) = 1 - uptick_prob;

	for (int i = 2; i <= 2 * no_of_divisions; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++) {
			if (j == (i - 1))
				transition_probability(i, j) = 1 - uptick_prob;
			if (j == (i + 1))
				transition_probability(i, j) = uptick_prob;
		}

	Matrix v_t(2 * no_of_divisions + 1, 1);
	// value at expiration
	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		v_t(i, 1) = max(0.0, (initial_stock_price * pow(up_factor, i - 1 - no_of_divisions)) - strike_price);

	// value at intermediate stages
	Matrix discounted_one_step_forward_value(2 * no_of_divisions + 1, 1);
	Matrix value_if_option_is_exercised_now(2 * no_of_divisions + 1, 1);
	for (int i = no_of_divisions; i > 0; i--) {
		// going backwards from expiration to zero-time
		discounted_one_step_forward_value = (transition_probability * v_t) / r;
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			value_if_option_is_exercised_now(j, 1) = max(0, (initial_stock_price * pow(up_factor, j - 1 - no_of_divisions) - strike_price));
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			v_t(j, 1) = max(value_if_option_is_exercised_now(j, 1), discounted_one_step_forward_value(j, 1));
	}

	// at the end of this backward movement, the middle-entry is the value of the american option

	return (v_t(no_of_divisions + 1, 1));
}

double american_put_option_dyn_prog()
{
	// create the probability Matrix
	Matrix transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions + 1);

	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			transition_probability(i, j) = 0.0;

	// boundary values of the probabilities need to be entered
	transition_probability(1, 1) = 1.0 - uptick_prob;
	transition_probability(1, 2) = uptick_prob;
	transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions + 1) = uptick_prob;
	transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions) = 1 - uptick_prob;

	for (int i = 2; i <= 2 * no_of_divisions; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++) {
			if (j == (i - 1))
				transition_probability(i, j) = 1 - uptick_prob;
			if (j == (i + 1))
				transition_probability(i, j) = uptick_prob;
		}

	Matrix v_t(2 * no_of_divisions + 1, 1);
	// value at expiration
	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		v_t(i, 1) = max(0.0, strike_price - (initial_stock_price * pow(up_factor, i - 1 - no_of_divisions)));
	// value at intermediate stages
	Matrix discounted_one_step_forward_value(2 * no_of_divisions + 1, 1);
	Matrix value_if_option_is_exercised_now(2 * no_of_divisions + 1, 1);
	for (int i = no_of_divisions; i > 0; i--) {
		// going backwards from expiration to zero-time
		discounted_one_step_forward_value = (transition_probability * v_t) / r;
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			value_if_option_is_exercised_now(j, 1) = max(0, strike_price - (initial_stock_price * pow(up_factor, j - 1 - no_of_divisions)));
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			v_t(j, 1) = max(value_if_option_is_exercised_now(j, 1), discounted_one_step_forward_value(j, 1));
	}

	// at the end of this backward movement, the middle-entry is the value of the american option

	return (v_t(no_of_divisions + 1, 1));
}

int main(int argc, char* argv[])
{

	cout << "Enter T,N,Rf,volatility(in decimal), S0, K: " << endl;
		cin >> expiration_time >> no_of_divisions >> risk_free_rate >> volatility >> initial_stock_price >> strike_price;

	up_factor = exp(volatility * sqrt(expiration_time / ((float)no_of_divisions)));
	r = exp(risk_free_rate * expiration_time / ((float)no_of_divisions));
	uptick_prob = (r - (1 / up_factor)) / (up_factor - (1 / up_factor));

	cout << "american option pricing by trinomial-model-inspired dynamic programming" << endl;
	cout << "expiration time (years) = " << expiration_time << endl;
	cout << "number of divisions = " << no_of_divisions << endl;
	cout << "risk free interest rate = " << risk_free_rate << endl;
	cout << "volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "initial stock price = " << initial_stock_price << endl;
	cout << "strike price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "r = " << r << endl;
	//cout << "up factor = " << up_factor << endl;
	//cout << "uptick probability = " << uptick_prob << endl;
	cout << "--------------------------------------" << endl;
	cout << "price of an american call option = " << american_call_option_dyn_prog() << endl;
	cout << "price of an american put option = " << american_put_option_dyn_prog() << endl;
	cout << "--------------------------------------" << endl;
}