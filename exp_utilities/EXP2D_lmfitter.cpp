#include <EXP2D_lmfitter.h>


// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*

    This is an example illustrating the use the general purpose non-linear 
    least squares optimization routines from the dlib C++ Library.

    This example program will demonstrate how these routines can be used for data fitting.
    In particular, we will generate a set of data and then use the least squares  
    routines to infer the parameters of the model which generated the data.
*/



// ----------------------------------------------------------------------------------------


// using namespace dlib;


// ----------------------------------------------------------------------------------------

// We will use this function to generate data.  It represents a function of 2 variables
// and 3 parameters.   The least squares procedure will be used to infer the values of 
// the 3 parameters based on a set of input/output pairs.
double lmfitter::model ( const input_vector& input, const parameter_vector& params )
{
    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);
    const double p3 = params(3);

    const double i0 = input(0);
    const double i1 = input(1);

    // const double temp = p0*i0 + p1*i1 + p2;
    // const double temp = p0 * exp(- i0*i0 / p1 - i1 * i1 / p2);
    // const double temp = p0 * exp(- p1 * i0 * i0 - p2 * i0 * i1 - p3 * i1 * i1);
    // const double temp = p0 - (p1 * i0 * i0 + p2 * i0 * i1 + p3 * i1 * i1);
    double temp = 0;

    double value = 2.0 * 10000 / (M_PI * p1 * p3) * (1 - (i0*i0)/(p1*p1) - (i1*i1)/(p3*p3) - p2 * i0 * i1);
    // double value = p0 * (1 - (i0*i0)/(p1*p1) - (i1*i1)/(p3*p3) - p2 * i0 * i1) ;
    if(value > 0){
        temp = sqrt(value);
    }

    // return temp*temp;
    return temp;
}

// ----------------------------------------------------------------------------------------

// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double lmfitter::residual ( const std::pair<input_vector, double>& data, const parameter_vector& params )
{
    return model(data.first, params) - data.second;
}

// ----------------------------------------------------------------------------------------

// This function is the derivative of the residual() function with respect to the parameters.
parameter_vector lmfitter::residual_derivative ( const std::pair<input_vector, double>& data, const parameter_vector& params )
{
    parameter_vector der;

    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);
    const double p3 = params(3);

    const double i0 = data.first(0);
    const double i1 = data.first(1);

    // const double temp = p0*i0 + p1*i1 + p2;
    // const double temp = exp(- p1 * i0 * i0 - p2 * i0 * i1 - p3 * i1 * i1);
    const double temp = (1 - (i0*i0)/(p1*p1) - (i1*i1)/(p3*p3) - p2 * i0 * i1);

	// der(0) = temp;
	// der(1) = - p0 * temp * ( i0*i0);
	// der(2) = - p0 * temp * ( i1*i0);
	// der(3) = - p0 * temp * ( i1*i1);

	der(0) = temp;
	der(1) = p0 * 2.0 * ( i0*i0) / (p1 * p1 * p1);
	der(2) = - p0 * ( i1*i0);
	der(3) = p0 * 2.0 * ( i1*i1) / (p3 * p3 * p3);


    // der(0) = i0*2*temp;
    // der(1) = i1*2*temp;
    // der(2) = 2*temp;

    return der;
}

parameter_vector lmfitter::set_initial_parameters()
{
		parameter_vector params;
	    int n = density.cols();
        // int m = density.rows();
        double coordinate_axis = meta.coord[0];
        // double delta_x = 2 * coordinate_axis / n;
        // double a,b,c;
        double peak = density(meta.grid[0]/2,meta.grid[1]/2);


        for(int i = 0; i < meta.grid[0]/2; ++i){
        	if( density(i,meta.grid[1]/2) >= peak * 0.005){
        		double tmp = ((meta.grid[0]/2) - i ) * meta.spacing[0];
        		// a = 0.5 * 2.35 / (tmp * tmp);
        		params(1) = tmp;
        		// cerr << "Rx: " << tmp << endl;
        		break;
        	}
        }

        for(int i = 0; i < meta.grid[1]/2; ++i){
        	if( density(i,meta.grid[0]/2) >= peak * 0.005){
        		double tmp = ((meta.grid[1]/2) - i ) * meta.spacing[0];
        		// c = 0.5 * 2.35 / (tmp * tmp);
        		params(3) = tmp;
        		// cerr << "Ry: " << tmp << endl;
        		break;
        	}
        }

        // double n0 = 2 * (10000.0 / M_PI) * (1 / (a * c));
        params(0) = 10000;
        // double n0 = peak;

        // b = (a + c) /2.0;
        params(2) = 0.0;
        return params;
}

// ----------------------------------------------------------------------------------------

vector<double> lmfitter::optimize()
{
    try
    {
        // randomly pick a set of parameters to use in this example
    	// parameter_vector params;
    	// params(0) = 1;
     //    params(1) = 10*randm(3,1);
     //    params(2) = 10*randm(3,1);
     //    params(3) = 10*randm(3,1);
     //    cout << "params: " << trans(params) << endl;


        // Now let's generate a bunch of input/output pairs according to our model.
        std::vector<std::pair<input_vector, double> > data_samples;
        input_vector input;
        // for (int i = 0; i < 160000; ++i)
        // {
        //     input = 20*randm(2,1) - 10;
        //     const double output = model(input, params);

        //     // save the pair
        //     data_samples.push_back(make_pair(input, output));
        // }


        

        for(int i = 0; i < meta.grid[0]; ++i){
        	for(int j = 0; j < meta.grid[1]; ++j){
        		if(density(i,j) > 0){
        			input(0) = - meta.coord[0] + meta.spacing[0] * i;
        			input(1) = - meta.coord[1] + meta.spacing[1] * j;

        			// const double output = model(input,params) + 0.1 * randm(1,1);
        			const double output = density(i,j);

        			data_samples.push_back(make_pair(input, output));
        		}

        	}
        }

        // plotPair(data_samples);

        // Before we do anything, let's make sure that our derivative function defined above matches
        // the approximate derivative computed using central differences (via derivative()).  
        // If this value is big then it means we probably typed the derivative function incorrectly.
        // cout << "derivative error: " << dlib::length(lmfitter::residual_derivative(data_samples[0], params) - dlib::derivative(lmfitter::residual)(data_samples[0], params) ) << endl;





        // Now let's use the solve_least_squares_lm() routine to figure out what the
        // parameters are based on just the data_samples.
        parameter_vector x = set_initial_parameters();
        parameter_vector params = x;

        // cerr << "Use Levenberg-Marquardt" << endl;
        // Use the Levenberg-Marquardt method to determine the parameters which
        // minimize the sum of all squared residuals.
        // solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-4,2e3).be_verbose(), residual, residual_derivative, data_samples, x);
        solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7,5e2)/*.be_verbose()*/, residual, derivative(residual), data_samples, x);

        // Now x contains the solution.  If everything worked it will be equal to params.
        cerr << "initial parameters: " << trans(params) << endl;
        cerr << "inferred parameters: "<< trans(x) << endl;
        // cerr << "solution error:      "<< dlib::length(x - params) << endl;
        cerr << endl;

        vector<double> paramter_results {x(0),x(1),x(2),x(3)};


        // plotGauss(params,n,coordinate_axis);
        plotPairAndGauss(data_samples,x,density.cols(),meta.coord[0]);






        // x = 1;
        // // cout << "Use Levenberg-Marquardt, approximate derivatives" << endl;
        // // // If we didn't create the residual_derivative function then we could
        // // // have used this method which numerically approximates the derivatives for you.
        // 

        // // // Now x contains the solution.  If everything worked it will be equal to params.
        // cout << "inferred parameters: "<< dlib::trans(x) << endl;
        // // // cout << "solution error:      "<< dlib::length(x - params) << endl;
        // cout << endl;




        // x = 1;
        // cerr << "Use Levenberg-Marquardt/quasi-newton hybrid" << endl;
        // // This version of the solver uses a method which is appropriate for problems
        // // where the residuals don't go to zero at the solution.  So in these cases
        // // it may provide a better answer.
        // dlib::solve_least_squares(dlib::objective_delta_stop_strategy(1e-20).be_verbose(), 
        //                     residual,
        //                     residual_derivative,
        //                     data_samples,
        //                     x);

        // // // Now x contains the solution.  If everything worked it will be equal to params.
        // cerr << "inferred parameters: "<< dlib::trans(x) << endl;
        // // cout << "solution error:      "<< dlib::length(x - params) << endl;


        // cerr << "press [Enter] to continue.." << endl;
        // std::cin.get();

        return paramter_results;

    }
    catch (std::exception& e)
    {
        cout << e.what() << endl;
    }
}

// Example output:
/*
params: 8.40188 3.94383 7.83099 

derivative error: 9.78267e-06
Use Levenberg-Marquardt
iteration: 0   objective: 2.14455e+10
iteration: 1   objective: 1.96248e+10
iteration: 2   objective: 1.39172e+10
iteration: 3   objective: 1.57036e+09
iteration: 4   objective: 2.66917e+07
iteration: 5   objective: 4741.9
iteration: 6   objective: 0.000238674
iteration: 7   objective: 7.8815e-19
iteration: 8   objective: 0
inferred parameters: 8.40188 3.94383 7.83099 

solution error:      0

Use Levenberg-Marquardt, approximate derivatives
iteration: 0   objective: 2.14455e+10
iteration: 1   objective: 1.96248e+10
iteration: 2   objective: 1.39172e+10
iteration: 3   objective: 1.57036e+09
iteration: 4   objective: 2.66917e+07
iteration: 5   objective: 4741.87
iteration: 6   objective: 0.000238701
iteration: 7   objective: 1.0571e-18
iteration: 8   objective: 4.12469e-22
inferred parameters: 8.40188 3.94383 7.83099 

solution error:      5.34754e-15

Use Levenberg-Marquardt/quasi-newton hybrid
iteration: 0   objective: 2.14455e+10
iteration: 1   objective: 1.96248e+10
iteration: 2   objective: 1.3917e+10
iteration: 3   objective: 1.5572e+09
iteration: 4   objective: 2.74139e+07
iteration: 5   objective: 5135.98
iteration: 6   objective: 0.000285539
iteration: 7   objective: 1.15441e-18
iteration: 8   objective: 3.38834e-23
inferred parameters: 8.40188 3.94383 7.83099 

solution error:      1.77636e-15
*/
