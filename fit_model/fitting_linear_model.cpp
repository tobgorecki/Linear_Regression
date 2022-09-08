#include "fitting_linear_model.h"
#include <math.h>
#include <qmath.h>
#include <QtMath>
#include <QDebug>
#include <QString>
#include <random>
#include <QtCore>
#include <thread>
#include <memory>
Fitting_linear_model::Fitting_linear_model(double init_A, double init_B)
: init_value_A(init_A), init_value_B(init_B)  {};


void Fitting_linear_model::Least_sqaure_classic(int Number_Data_Points, std::vector<double> x_value,std::vector<double> y_value,
                            double *parameter_A, double *parameter_B,
                            double *parameter_A_error, double *parameter_B_error,
                            double *correlation_coefficient ){

    double sum_value_error_x = 0.0;   // x
    double sum_value_error_y = 0.0;   // y
    double sum_value_error_xy = 0.0;  // x*y
    double sum_value_error_xx = 0.0;  // x*x
    double sum_value_error_yy = 0.0;  // x*x


    for (int i=0 ; i < Number_Data_Points; ++i){
            sum_value_error_x += x_value[i];
            sum_value_error_y += y_value[i];
            sum_value_error_xy += x_value[i] *  y_value[i];
            sum_value_error_xx += x_value[i] *  x_value[i];
            sum_value_error_yy += y_value[i] *  y_value[i];
    }

     *parameter_A = ( Number_Data_Points * sum_value_error_xy - sum_value_error_x *sum_value_error_y )/
            ( Number_Data_Points  *sum_value_error_xx - qPow(sum_value_error_x, 2) );

     *parameter_B = ( sum_value_error_y -(*parameter_A) *sum_value_error_x )/ Number_Data_Points;

    if ( Number_Data_Points > 2 ){
        *parameter_A_error= qSqrt( Number_Data_Points/(Number_Data_Points - 2) *
                             ( sum_value_error_yy  -  (*parameter_A )*sum_value_error_xy - (*parameter_B)*sum_value_error_y )/
                             ( Number_Data_Points  *sum_value_error_xx - qPow(sum_value_error_x, 2) )
                              );
    }else{

        *parameter_A_error=0.0;
    }

    *parameter_B_error =  (*parameter_A_error) *  qSqrt( sum_value_error_xx/Number_Data_Points);


     *correlation_coefficient= ( Number_Data_Points* sum_value_error_xy - sum_value_error_x *  sum_value_error_y )/
                              qSqrt(
                                    ( Number_Data_Points  *sum_value_error_xx - qPow(sum_value_error_x, 2) )*
                                     ( Number_Data_Points  *sum_value_error_yy - qPow(sum_value_error_y, 2) )
                                    );
}


void Fitting_linear_model::Least_sqaure_weighted(int Number_Data_Points,std::vector<double> x_value,std::vector<double> y_value,
                        std::vector<double> y_value_error,
                        double *parameter_A, double *parameter_B,
                        double *parameter_A_error, double *parameter_B_error,
                        double *correlation_coefficient)
{


    std::vector<double> y_value_error_weighted;
    double sum_y_value_error_weighted = 0.0; //  w=(1/y)^2
    double sum_value_error_weighted_and_x = 0.0;   // w *x
    double sum_value_error_weighted_and_y = 0.0;  // w *y
    double sum_value_error_weighted_and_xy = 0.0;  // w *x*y
    double sum_value_error_weighted_and_xx = 0.0;  // w *x*x
    double sum_value_error_weighted_and_yy = 0.0;  // w *x*x


    for (auto it =y_value_error.begin(); it != y_value_error.end(); ++it)
     y_value_error_weighted.push_back(qSqrt(1.0/ *it));



    for (auto i: y_value_error_weighted)
            sum_y_value_error_weighted += i;

    for (int i=0 ; i < Number_Data_Points; ++i){
            sum_value_error_weighted_and_x += x_value[i] *y_value_error_weighted[i];
            sum_value_error_weighted_and_y += y_value[i] *y_value_error_weighted[i];
            sum_value_error_weighted_and_xy += x_value[i] *  y_value[i]* y_value_error_weighted[i];
            sum_value_error_weighted_and_xx += x_value[i] *  x_value[i]* y_value_error_weighted[i];
            sum_value_error_weighted_and_yy += y_value[i] *  y_value[i]* y_value_error_weighted[i];
    }
    // Y= A*x + B

    *parameter_A = ( sum_y_value_error_weighted * sum_value_error_weighted_and_xy  - sum_value_error_weighted_and_x * sum_value_error_weighted_and_y )/
            ( sum_y_value_error_weighted * sum_value_error_weighted_and_xx - qPow(sum_value_error_weighted_and_x,2) );

    *parameter_B = (sum_value_error_weighted_and_y  - *parameter_A * sum_value_error_weighted_and_x) / sum_y_value_error_weighted ;

    // computing error for parameters A i B
    *parameter_A_error = qSqrt(  sum_y_value_error_weighted /( sum_y_value_error_weighted * sum_value_error_weighted_and_xx - qPow( sum_value_error_weighted_and_x,2) ) );
    *parameter_B_error = (*parameter_A_error) * qSqrt( sum_value_error_weighted_and_xx /sum_y_value_error_weighted  );


    *correlation_coefficient =( sum_y_value_error_weighted * sum_value_error_weighted_and_xy  - sum_value_error_weighted_and_x * sum_value_error_weighted_and_y )/
            qSqrt(
                ( sum_y_value_error_weighted * sum_value_error_weighted_and_xx - qPow( sum_value_error_weighted_and_x,2) ) *
                ( sum_y_value_error_weighted * sum_value_error_weighted_and_yy - qPow( sum_value_error_weighted_and_y,2) )
                 );


}


void  Fitting_linear_model::Gradinet_descent (int Number_Data_Points,std::vector<double> x_value,std::vector<double> y_value,
                                             std::vector<double> y_value_error,
                                             double*parameter_A, double *parameter_B,
                                             double init_value_A,  double init_value_B )
{


    double lambda =1.0e-6;  //for good accuracy.



    // the partial derivative of the loos function, respect to parameter A and B
    // derivative_A = -2* sum x_n(y_n -  Ax_n -B)/sigma
    // derivative_A = -2* sum (y_n -  Ax_n -B)/sigma


    *parameter_A = init_value_A;
    *parameter_B = init_value_B;

    double derivative_A = 0.0;
    double derivative_B = 0.0;


    std::vector<double>loss_function;
    double model_old=0;
    double model_new=0;
    int iter=0;
  while(true)
    {



      for (int n = 0; n <Number_Data_Points; ++n){
                derivative_A +=  - 2* x_value[n] * (y_value[n] - (*parameter_A) *x_value[n] - (*parameter_B) )/Number_Data_Points ;
                derivative_B +=  - 2* (y_value[n] - (*parameter_A) *x_value[n] - (*parameter_B) )/Number_Data_Points;
           }

        *parameter_A = *parameter_A - lambda * derivative_A ;
        *parameter_B = *parameter_B - lambda * derivative_B ;



       model_new=0;
      for (int n = 0; n <Number_Data_Points; ++n){
             model_new+=  qPow(y_value[n] - (*parameter_A) *x_value[n] - (*parameter_B) ,2)/ Number_Data_Points ;
        }


        loss_function.push_back(qPow( model_old -  model_new,2));
        if (iter> 200){
            if (abs(loss_function[iter] -  loss_function[iter-100])<lambda){
                break;
                }


        }

         model_old = model_new;



        iter+=1;
    }
}


 void Fitting_linear_model::One_chain_mcmc(int Number_Data_Points, std::vector<double> x_value,std::vector<double> y_value,
                                          std::vector<double> y_value_error,
                                            double *parameter_A, double *parameter_B,
                                          double init_value_A,  double init_value_B )
{
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    std::normal_distribution<>norm(0, .01);
    int Number_iteration_mcmc=10000;
    std::vector<double> par_A(Number_iteration_mcmc);
    std::vector<double> par_B(Number_iteration_mcmc);

    double prob=0.0;
    double accept=0.0;
    par_A[0] = init_value_A;
    par_B[0] = init_value_B;

    double A =  par_A[0]+ norm(e2);
    double B =  par_B[0]+ norm(e2);

    int iter=0;
    int cou=0;
    while (iter < Number_iteration_mcmc )
    {

        while(true){

         prob =exp(log_likehood( Number_Data_Points,  x_value, y_value, y_value_error, A,  B)-
                log_likehood( Number_Data_Points,  x_value, y_value, y_value_error, par_A[iter],  par_B[iter]));

        accept  = fmin(1,prob);
        if ( dist(e2) <= accept ){
            par_A[iter+1] = A;
            par_B[iter+1] = B ;
            A =  par_A[iter+1]+ norm(e2);
            B =  par_B[iter+1]+ norm(e2);
           cou+=1;

            break;
         }
        else{
               A =  par_A[iter]+ norm(e2);
               B =  par_B[iter]+ norm(e2);
            }

       }

        iter +=1 ;
    }


 *parameter_A = accumulate(par_A.begin(), par_A.end(), 0.0)/par_A.size();
 *parameter_B = accumulate(par_B.begin(), par_B.end(), 0.0)/par_B.size();

}

void Fitting_linear_model::Monte_Carto_Markov_Chains(int Number_Data_Points, std::vector<double> x_value,std::vector<double> y_value,
                                                     std::vector<double> y_value_error,
                                                     double *parameter_A, double *parameter_B,
                                                     double init_value_A,  double init_value_B )
{



    //number of available processors minus 2
    const auto processor_count =std::thread::hardware_concurrency() - 2;

   // table contains smart pointers
   std::unique_ptr<double[]> par_A =  std::make_unique <double[]>(processor_count) ;
   std::unique_ptr<double[]> par_B =  std::make_unique <double[]>(processor_count) ;
   std::unique_ptr<std::thread[]> chain  =  std::make_unique <std::thread[]>(processor_count) ;

    *parameter_A = 0.0;
    *parameter_B = 0.0;


   for (size_t  i=0; i <  processor_count;  ++i ){

     chain[i] = std::thread(One_chain_mcmc, Number_Data_Points, x_value, y_value,y_value_error,  &par_A[i], &par_B[i], init_value_A,  init_value_B );

   }

   for (size_t  i=0; i <processor_count; ++i ){
        chain[i].join();
   }

   for (size_t i=0 ; i < processor_count; ++i){
       *parameter_A = *parameter_A + par_A[i];
       *parameter_B = *parameter_B + par_B[i];
    }
   *parameter_A = *parameter_A/processor_count;
   *parameter_B = *parameter_B/processor_count;



}

double Fitting_linear_model::log_likehood(int Number_Data_Points, std::vector<double> x_value,std::vector<double> y_value,std::vector<double> y_error,double parameter_A, double parameter_B)
{
    double log_likehood_result=0;


    for (int i=0; i <Number_Data_Points; ++i){
        log_likehood_result +=   fabs (pow(y_value[i] - parameter_A*  x_value[i] - parameter_B,2)/ pow(y_error[i],2) + log(2.0*M_PI*pow(y_error[i],2)));
    }


    return  -0.5 *log_likehood_result;
}
