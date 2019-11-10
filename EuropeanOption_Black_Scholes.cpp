//// Black-Scholes European Option Pricing Code
//// Adapted from Prof. Odegaard's Notes
//
//#include <iostream>
//#include <iomanip>
//#include <cmath>
//#include <fstream>
//#include <cstdlib>
//#include "./normdist.h"          // this defines the normal distribution from Odegaard's files
//using namespace std;
//
//double risk_free_rate, strike_price;
//double initial_stock_price, expiration_time, volatility;
//
//double option_price_put_black_scholes(const double& S,      // spot price
//	const double& K,      // Strike (exercise) price,
//	const double& r,      // interest rate
//	const double& sigma,  // volatility
//	const double& time) {
//	double time_sqrt = sqrt(time);
//	double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
//	double d2 = d1 - (sigma * time_sqrt);
//	return K * exp(-r * time) * N(-d2) - S * N(-d1);
//}
//
//double option_price_call_black_scholes(const double& S,       // spot (underlying) price
//	const double& K,       // strike (exercise) price,
//	const double& r,       // interest rate
//	const double& sigma,   // volatility 
//	const double& time) {  // time to maturity 
//	double time_sqrt = sqrt(time);
//	double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
//	double d2 = d1 - (sigma * time_sqrt);
//	return S * N(d1) - K * exp(-r * time) * N(d2);
//}
//
//double N(const double& z) {
//	if (z > 6.0) { return 1.0; }; // this guards against overflow 
//	if (z < -6.0) { return 0.0; };
//	double b1 = 0.31938153;
//	double b2 = -0.356563782;
//	double b3 = 1.781477937;
//	double b4 = -1.821255978;
//	double b5 = 1.330274429;
//	double p = 0.2316419;
//	double c2 = 0.3989423;
//	double a = fabs(z);
//	double t = 1.0 / (1.0 + a * p);
//	double b = c2 * exp((-z) * (z / 2.0));
//	double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
//	n = 1.0 - b * n;
//	if (z < 0.0) n = 1.0 - n;
//	return n;
//}
//
//int main(int argc, char* argv[])
//{
//
//	cout << "Enter T,Rf,volatility(in decimal), S0, K: " << endl;
//	cin >> expiration_time >> risk_free_rate >> volatility >> initial_stock_price >> strike_price;
//
//
//	cout << "Black-Scholes European Option Pricing" << endl;
//	cout << "Expiration Time (Years) = " << expiration_time << endl;
//	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
//	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
//	cout << "Initial Stock Price = " << initial_stock_price << endl;
//	cout << "Strike Price = " << strike_price << endl;
//	cout << "--------------------------------------" << endl;
//	double call_price = option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time);
//	double put_price = option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time);
//	cout << "Call Price according to Black-Scholes = " << call_price << endl;
//	cout << "Put Price according to Black-Scholes = " << put_price << endl;
//	cout << "--------------------------------------" << endl;
//	cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
//	cout << initial_stock_price << " + " << put_price << " - " << call_price;
//	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
//	cout << initial_stock_price + put_price - call_price << " = " << strike_price * exp(-risk_free_rate * expiration_time) << endl;
//	cout << "--------------------------------------" << endl;
//}
