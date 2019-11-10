#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

double american_call_option(int k, int i, double current_stock_price)
{
	if (k == no_of_divisions)
		return max(0.0, (current_stock_price - strike_price));
	else
		return max((current_stock_price - strike_price),
		(uptick_prob * american_call_option(k + 1, i + 1, current_stock_price * up_factor) +
			(1 - uptick_prob) * american_call_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
}

double american_put_option(int k, int i, double current_stock_price)
{
	if (k == no_of_divisions)
		return max(0.0, (strike_price - current_stock_price));
	else
		return max((strike_price - current_stock_price),
		(uptick_prob * american_put_option(k + 1, i + 1, current_stock_price * up_factor) +
			(1 - uptick_prob) * american_put_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
}

int main(int argc, char* argv[])
{

	cout << "Enter T,N,Rf,volatility(in decimal), S0, K: " << endl;
    cin >> expiration_time >> no_of_divisions >> risk_free_rate >> volatility >> initial_stock_price >> strike_price;

	up_factor = exp(volatility * sqrt(expiration_time / ((float)no_of_divisions)));
	R = exp(risk_free_rate * expiration_time / ((float)no_of_divisions));
	uptick_prob = (R - (1 / up_factor)) / (up_factor - (1 / up_factor));
	cout << "Recursive Binomial American Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "--------------------------------------" << endl;
	double call_price = american_call_option(0, 0, initial_stock_price);
	cout << "Binomial Price of an American Call Option = " << call_price << endl;
	double put_price = american_put_option(0, 0, initial_stock_price);
	cout << "Binomial Price of an American Put Option = " << put_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Let us verify the Put-Call Parity: S+P-C = Kexp(-r*T) for American Options" << endl;
	cout << initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price - call_price << " ?=? " << strike_price * exp(-risk_free_rate * expiration_time) << endl;
	if (abs(initial_stock_price + put_price - call_price - strike_price * exp(-risk_free_rate * expiration_time)) <= 1e-3)
		cout << "Looks like Put-Call Parity holds within three decimal places" << endl;
	else
		cout << "Looks like Put-Call Parity does NOT hold" << endl;
	cout << "--------------------------------------" << endl;

}
