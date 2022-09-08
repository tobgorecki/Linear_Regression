#ifndef FITTING_LINEAR_MODEL_H
#define FITTING_LINEAR_MODEL_H
#include <vector>

//############ A comparison of three linear regression methods. ######################

// class cointains three( maybe more in future) methods to fitting linaer regeresion
class Fitting_linear_model
{
private :
    int Number_Data_Points;
     // initial conditions, required for only GD and MCMC methods.
    double init_value_A;
    double init_value_B;

    // Vectors with data to fit. Data will be upload by user.
    // VERY IMPORTANT - x values must be monotonically increase, otherwise there will be problem
    std::vector<double> x_value;
    std::vector<double> y_value;
    // error for y values
    std::vector<double> y_error;

    void check_input_data( std::vector<double> x_value,std::vector<double> y_value,
                           std::vector<double> y_value_error );
public:
    // parameters with results
    double *parameter_A;
    double *parameter_B;
    double *parameter_A_error;
    double *parameter_B_error;
    double *correlation_coefficient;

    Fitting_linear_model(double init_A, double init_B);

    // These functions are using to fitting linear model.
    void Least_sqaure_classic(int Number_Data_Points, std::vector<double> x_value,std::vector<double> y_value,
                            double *parameter_A, double *parameter_B, double *parameter_A_error, double *parameter_B_error, double *correlation_coefficient );
    void Least_sqaure_weighted(int Number_Data_Points, std::vector<double> x_value,std::vector<double> y_value,
                            std::vector<double> y_value_error,
                            double *parameter_A, double *parameter_B, double *parameter_A_error, double *parameter_B_error, double *correlation_coefficient );

    void Gradinet_descent(int Number_Data_Points,std::vector<double> x_value,std::vector<double> y_value,
                          std::vector<double> y_value_error,
                          double*parameter_A, double *parameter_B,
                          double init_value_A,  double init_value_B );

    // computing mcmc in parallel mode (thread c++)
    void Monte_Carto_Markov_Chains(int Number_Data_Points, std::vector<double> x_value,std::vector<double> y_value,
                                   std::vector<double> y_value_error,
                                   double  * parameter_A,double * parameter_B,
                                   double init_value_A,  double init_value_B );
    static void One_chain_mcmc(int Number_Data_Points, std::vector<double> x_value,std::vector<double> y_value,
                                              std::vector<double> y_value_error,
                                              double *parameter_A, double *parameter_B,
                                              double init_value_A,  double init_value_B );
    static double log_likehood(int Number_Data_Points, std::vector<double> x_value,std::vector<double> y_value,
                            std::vector<double> y_error,
                         double parameter_A, double parameter_B);

};

#endif // FITTING_LINEAR_MODEL_H
