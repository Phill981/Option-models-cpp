#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <chrono>
#include "gnuplot-iostream.h"  // Gnuplot-iostream for plotting
#include <functional>

// Integration function using the Trapezoidal rule
double integrate(std::function<double(double)> func, double a, double b, int n = 1000) {
    double h = (b - a) / n;  // Step size
    double sum = 0.5 * (func(a) + func(b));  // End points contribution

    for (int i = 1; i < n; ++i) {
        sum += func(a + i * h);  // Midpoint contributions
    }

    return sum * h;
}

// Black-Scholes Call Price Function
double black_scholes_call(double S, double K, double T, double r, double sigma) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    
    auto norm_cdf = [](double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); };
    
    double call_price = S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
    return call_price;
}

// Merton Jump Diffusion Model Call Price Function
double merton_jump_diffusion_call(double S, double K, double T, double r, double sigma, double lambda, 
                                  double muJ, double sigmaJ) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);

    auto norm_cdf = [](double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); };
    double BS_price = S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);

    double jump_part = std::exp(lambda * (std::exp(muJ + 0.5 * sigmaJ * sigmaJ) - 1)) * 
                       (std::exp(muJ + 0.5 * sigmaJ * sigmaJ) * norm_cdf(d1) - norm_cdf(d2));
    
    return BS_price + jump_part;
}

// Function to measure run time of Black-Scholes
double measure_runtime(double (*func)(double, double, double, double, double), double S, double K, double T, double r, double sigma) {
    auto start = std::chrono::high_resolution_clock::now();
    func(S, K, T, r, sigma);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return duration.count();
}

// Function to measure run time of Merton
double measure_runtime(double (*func)(double, double, double, double, double, double, double, double), double S, double K, double T, double r, double sigma, double lambda, double muJ, double sigmaJ) {
    auto start = std::chrono::high_resolution_clock::now();
    func(S, K, T, r, sigma, lambda, muJ, sigmaJ);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return duration.count();
}

// Plot Sensitivity and Performance Analysis Function
void plot_sensitivity_and_performance(double S, double K, double T, double r, double sigma, double lambda, 
                                       double muJ, double sigmaJ) {
    Gnuplot gp;

    std::vector<std::pair<double, double>> bs_vol_prices, merton_vol_prices;
    std::vector<std::pair<double, double>> bs_maturities, merton_maturities;
    std::vector<std::pair<double, double>> bs_strikes, merton_strikes;
    std::vector<std::pair<double, double>> price_discrepancy, run_times;

    // Volatility sensitivity
    for (double vol = 0.01; vol <= 0.6; vol += 0.01) {
        bs_vol_prices.emplace_back(vol, black_scholes_call(S, K, T, r, vol));
        merton_vol_prices.emplace_back(vol, merton_jump_diffusion_call(S, K, T, r, vol, lambda, muJ, sigmaJ));
    }
    // Time to maturity sensitivity
    for (double t = 0.01; t <= 2.0; t += 0.01) {
        bs_maturities.emplace_back(t, black_scholes_call(S, K, t, r, sigma));
        merton_maturities.emplace_back(t, merton_jump_diffusion_call(S, K, t, r, sigma, lambda, muJ, sigmaJ));
    }
    // Strike price sensitivity
    for (double strike = 50; strike <= 150; strike += 1.0) {
        bs_strikes.emplace_back(strike, black_scholes_call(S, strike, T, r, sigma));
        merton_strikes.emplace_back(strike, merton_jump_diffusion_call(S, strike, T, r, sigma, lambda, muJ, sigmaJ));
    }
    
    // Calculate price discrepancies
    for (size_t i = 0; i < bs_vol_prices.size(); ++i) {
        double discrepancy = std::fabs(bs_vol_prices[i].second - merton_vol_prices[i].second);
        price_discrepancy.emplace_back(bs_vol_prices[i].first, discrepancy);
    }

    // Measure run times
    double avg_bs_runtime = measure_runtime(black_scholes_call, S, K, T, r, sigma);
    double avg_merton_runtime = measure_runtime(merton_jump_diffusion_call, S, K, T, r, sigma, lambda, muJ, sigmaJ);
    
    run_times.emplace_back(0, avg_bs_runtime);  // Placeholder for Black-Scholes
    run_times.emplace_back(1, avg_merton_runtime); // Placeholder for Merton
    
    // Set multiplot for 3x2 layout
    gp << "set multiplot layout 3,2 title 'Sensitivity and Performance Analysis'\n";

    // Plot Call Price vs Volatility
    gp << "set xlabel 'Volatility'\n";
    gp << "set ylabel 'Call Price'\n";
    gp << "$data_vol_bs << EOD\n";
    for (const auto& point : bs_vol_prices) {
        gp << point.first << " " << point.second << "\n";
    }
    gp << "EOD\n";

    gp << "$data_vol_merton << EOD\n";
    for (const auto& point : merton_vol_prices) {
        gp << point.first << " " << point.second << "\n";
    }
    gp << "EOD\n";

    gp << "set title 'Call Price vs Volatility'\n";
    gp << "plot $data_vol_bs with lines title 'Black-Scholes', $data_vol_merton with lines title 'Merton'\n";

    // Plot Call Price vs Time to Maturity
    gp << "set xlabel 'Time to Maturity'\n";
    gp << "set ylabel 'Call Price'\n";
    gp << "$data_time_bs << EOD\n";
    for (const auto& point : bs_maturities) {
        gp << point.first << " " << point.second << "\n";
    }
    gp << "EOD\n";

    gp << "$data_time_merton << EOD\n";
    for (const auto& point : merton_maturities) {
        gp << point.first << " " << point.second << "\n";
    }
    gp << "EOD\n";

    gp << "set title 'Call Price vs Time to Maturity'\n";
    gp << "plot $data_time_bs with lines title 'Black-Scholes', $data_time_merton with lines title 'Merton'\n";

    // Plot Call Price vs Strike Price for Black-Scholes
    gp << "set xlabel 'Strike Price'\n";
    gp << "set ylabel 'Call Price'\n";
    gp << "$data_strike_bs << EOD\n";
    for (const auto& point : bs_strikes) {
        gp << point.first << " " << point.second << "\n";
    }
    gp << "EOD\n";

    gp << "set title 'Call Price vs Strike Price - Black-Scholes'\n";
    gp << "plot $data_strike_bs with lines title 'Black-Scholes'\n";

    // Plot Call Price vs Strike Price for Merton
    gp << "set xlabel 'Strike Price'\n";
    gp << "set ylabel 'Call Price'\n";
    gp << "$data_strike_merton << EOD\n";
    for (const auto& point : merton_strikes) {
        gp << point.first << " " << point.second << "\n";
    }
    gp << "EOD\n";

    gp << "set title 'Call Price vs Strike Price - Merton'\n";
    gp << "plot $data_strike_merton with lines title 'Merton'\n";

    // Plot Price Discrepancy
    gp << "set xlabel 'Volatility'\n";
    gp << "set ylabel 'Price Discrepancy'\n";
    gp << "$data_discrepancy << EOD\n";
    for (const auto& point : price_discrepancy) {
        gp << point.first << " " << point.second << "\n";
    }
    gp << "EOD\n";

    gp << "set title 'Price Discrepancy Between Black-Scholes and Merton'\n";
    gp << "plot $data_discrepancy with lines title 'Discrepancy'\n";

    // Plot Run Time Comparison
    gp << "set xlabel 'Model'\n";
    gp << "set ylabel 'Run Time (seconds)'\n";
    gp << "$data_runtime << EOD\n";
    for (const auto& point : run_times) {
        gp << point.first << " " << point.second << "\n";
    }
    gp << "EOD\n";

    gp << "set title 'Run Time Comparison'\n";
    gp << "set xtics ('Black-Scholes' 0, 'Merton' 1)\n";
    gp << "plot $data_runtime with boxes title 'Run Time'\n";

    // Finish multiplot
    gp << "unset multiplot\n";
}


int main() {
    double S = 100;             // Stock price
    double K = 100;             // Strike price
    double T = 1;               // Time to maturity (in years)
    double r = 0.05;            // Risk-free rate
    double sigma = 0.2;         // Volatility for Black-Scholes and Merton
    double lambda = 0.1;        // Jump intensity for Merton
    double muJ = 0.0;           // Mean jump size for Merton
    double sigmaJ = 0.2;        // Jump size volatility for Merton
    
    plot_sensitivity_and_performance(S, K, T, r, sigma, lambda, muJ, sigmaJ);
    return 0;
}
