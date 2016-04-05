#include "lmfitter.h"


double lmfitter::gauss_dist ( const input_vector& input, const parameter_vector& params )
{
    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);
    const double p3 = params(3);

    const double i0 = input(0);
    const double i1 = input(1);

    const double temp = p0 * exp(- p1 * i0 * i0 - p2 * i0 * i1 - p3 * i1 * i1);
    return temp;
}

// This function is the derivative of the residual() function with respect to the parameters.
parameter_vector lmfitter::gauss_residual_derivative ( const std::pair<input_vector, double>& data, const parameter_vector& params )
{
    parameter_vector der;

    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);
    const double p3 = params(3);

    const double i0 = data.first(0);
    const double i1 = data.first(1);

    const double temp = exp(- p1 * i0 * i0 - p2 * i0 * i1 - p3 * i1 * i1);

    der(0) = temp;
    der(1) = - p0 * temp * ( i0*i0);
    der(2) = - p0 * temp * ( i1*i0);
    der(3) = - p0 * temp * ( i1*i1);

    return der;
}

double lmfitter::gauss_residual ( const std::pair<input_vector, double>& data, const parameter_vector& params )
{
    return gauss_dist(data.first, params) - data.second;
}

double lmfitter::tf_dist ( const input_vector& input, const parameter_vector& params )
{
    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);
    const double p3 = params(3);

    const double i0 = input(0);
    const double i1 = input(1);

    double temp = 0;

    double value = p0 * (1 - (i0*i0)/(p1*p1) - (i1*i1)/(p3*p3) - p2 * i0 * i1);
    if(value > 0){
        temp = value;
    }

    return temp;
}

parameter_vector lmfitter::tf_residual_derivative ( const std::pair<input_vector, double>& data, const parameter_vector& params )
{
    parameter_vector der;

    const double p0 = params(0);
    const double p1 = params(1);
    const double p2 = params(2);
    const double p3 = params(3);

    const double i0 = data.first(0);
    const double i1 = data.first(1);

    const double temp = (1 - (i0*i0)/(p1*p1) - (i1*i1)/(p3*p3) - p2 * i0 * i1);

	if(temp > 0){

		der(0) = temp;
		der(1) = p0 * 2.0 * ( i0*i0) / (p1 * p1 * p1);
		der(2) = - p0 * ( i1*i0);
		der(3) = p0 * 2.0 * ( i1*i1) / (p3 * p3 * p3);
	} else {
		der = 0;
	}

    return der;
}

double lmfitter::tf_residual ( const std::pair<input_vector, double>& data, const parameter_vector& params )
{
    return tf_dist(data.first, params) - data.second;
}

parameter_vector lmfitter::set_initial_parameters()
{
		parameter_vector params;
	    int n = density.cols();
        double coordinate_axis = meta.coord[0];

        double last_largest_element = 0;
        std::queue<double> largest_elements;
        for(int i = 0; i < meta.grid[0]; ++i){
        	for(int j = 0; j < meta.grid[1]; ++j){
        		const double val = density(i,j);
        		if(val > last_largest_element){
        			largest_elements.push(val);
        			last_largest_element = val;
        		}
        		if(largest_elements.size() > /*meta.grid[0]*meta.grid[1]*0.001*/ 10){
        			largest_elements.pop();
        		}
        	}
        }
        double n0 = 0;
        int le_size = largest_elements.size();
        while(!largest_elements.empty()){
        	n0 += largest_elements.front();
        	largest_elements.pop();
        }
        n0 /= le_size;



        for(int i = 0; i < meta.grid[0]/2; ++i){
        	if( density(i,meta.grid[1]/2) >= n0 * 0.05){
        		double tmp = ((meta.grid[0]/2) - i ) * meta.spacing[0];
        		params(1) = tmp;
        		break;
        	}
        }

        for(int i = 0; i < meta.grid[1]/2; ++i){
        	if( density(i,meta.grid[0]/2) >= n0 * 0.05){
        		double tmp = ((meta.grid[1]/2) - i ) * meta.spacing[0];
        		params(3) = tmp;
        		break;
        	}
        }

        params(0) = n0;
        params(2) = 0.0;
        return params;
}


void lmfitter::plotQuerschnitte(const string& str, const parameter_vector& params_tf)
{
	    ArrayXd querschnitt_x(meta.grid[0]);
		for(int i = 0; i < meta.grid[0]; ++i){
			double value =  density(i,meta.grid[1]/2);
			if(value <= 100.0)
				querschnitt_x(i) = value;
			else 
				querschnitt_x(i) = 100.0;
		}

		ArrayXd fitschnitt_x(meta.grid[0]);
		for(int i = 0; i < meta.grid[0]; ++i){
			double i0 = - meta.coord[0] + meta.spacing[0] * i;
			double i1 = 0.0;
			double value = /*2 * (params_tf(0) / M_PI) * (1 / (params_tf(1) * params_tf(3)))*/ params_tf(0) * (1 - (i0*i0)/(params_tf(1)*params_tf(1)) - (i1*i1)/(params_tf(3)*params_tf(3)) - params_tf(2) * i0 * i1) ;
			if(value < 0) value = 0.0;
			if(value <= 100.0)
				fitschnitt_x(i) = value;
			else 
				fitschnitt_x(i) = 100.0;
			

		}
		

		ArrayXd querschnitt_y(meta.grid[0]);
		for(int i = 0; i < meta.grid[0]; ++i){
			double value =  density(meta.grid[0]/2,i);
			if(value <= 100.0)
				querschnitt_y(i) = value;
			else 
				querschnitt_y(i) = 100.0;
		}

		ArrayXd fitschnitt_y(meta.grid[0]);
		for(int i = 0; i < meta.grid[0]; ++i){
			double i1 = - meta.coord[0] + meta.spacing[0] * i;
			double i0 = 0.0;
			double value = /*2 * (params_tf(0) / M_PI) * (1 / (params_tf(1) * params_tf(3)))*/ params_tf(0) * (1 - (i0*i0)/(params_tf(1)*params_tf(1)) - (i1*i1)/(params_tf(3)*params_tf(3)) - params_tf(2) * i0 * i1) ;
			if(value < 0) value = 0.0;
			if(value <= 100.0)
				fitschnitt_y(i) = value;
			else 
				fitschnitt_y(i) = 100.0;
			

		}
		plotTwoVectors("Querschnitt_" + to_string(meta.steps) + str,"Querschnitt " +str,querschnitt_x,fitschnitt_x,querschnitt_y,fitschnitt_y);
}

// ----------------------------------------------------------------------------------------

vector<double> lmfitter::optimize()
{
    try
    {
        std::vector<std::pair<input_vector, double> > data_samples;
        input_vector input;

        parameter_vector x = set_initial_parameters();
        parameter_vector params = x;
        

        for(int i = 0; i < meta.grid[0]; ++i){
        	for(int j = 0; j < meta.grid[1]; ++j){
        		// if(density(i,j) < 0.5 * params(0)){
        			input(0) = - meta.coord[0] + meta.spacing[0] * i;
        			input(1) = - meta.coord[1] + meta.spacing[1] * j;

        			// const double output = model(input,params) + 0.1 * randm(1,1);
        			const double output = density(i,j);

        			data_samples.push_back(make_pair(input, output));
        		// }

        	}
        }
 
        solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-20,5e2)/*.be_verbose()*/, tf_residual, tf_residual_derivative, data_samples, x);
        // solve_least_squares_lm(dlib::gradient_norm_stop_strategy(1e-7,5e2)/*dlib::objective_delta_stop_strategy(1e-7,5e2)*//*.be_verbose()*/, residual, residual_derivative, data_samples, x);
        // solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-7,5e2)/*.be_verbose()*/, residual, derivative(residual), data_samples, x);

        cout << currentTime() << " Fit density with least squares Levenberg-Marquardt" << endl;
        cout << "\t initial parameters: " << "\t" << trans(params);
        cout << "\t fitted parameters: " << trans(x)/* << endl*/;
        // cerr << "solution error:      "<< dlib::length(x - params) << endl;


        vector<double> paramter_results {x(0),x(1),x(2),x(3)};

        return paramter_results;

    }
    catch (std::exception& e)
    {
        cout << e.what() << endl;
    }
}