#include <EXP2D_rk4.hpp>

RungeKutta::RungeKutta(/*MatrixData* &d,*/ Options &extOpt) : /*w(d),*/ opt(extOpt),
											laplacian_coefficient_x(3),
											laplacian_coefficient_y(3),
											gradient_coefficient_x(3),
											gradient_coefficient_y(3) {}

void RungeKutta::setVariables(){
	wavefctcp = MatrixXcd::Zero(w->meta.grid[0],w->meta.grid[1]);
	k0 = MatrixXcd::Zero(w->meta.grid[0],w->meta.grid[1]);
	k1 = MatrixXcd::Zero(w->meta.grid[0],w->meta.grid[1]);
	k2 = MatrixXcd::Zero(w->meta.grid[0],w->meta.grid[1]);
	k3 = MatrixXcd::Zero(w->meta.grid[0],w->meta.grid[1]);

	threads = omp_get_max_threads();
	subx = w->meta.grid[0]-4;
	suby = w->meta.grid[1]-4;
	frontx = vector<int32_t>(threads);
	endx = vector<int32_t>(threads);
	partx = w->meta.grid[0] / threads;

	for(int i = 0; i < threads; i++){
		if(i == 0){ frontx[i] = 2;}
		else{ frontx[i] = (i *partx);}
		if((i == threads-1) || (i == 0)){ endx[i] = partx-2;}
		else{endx[i] = partx;}
	}

	// Coordinate vectors/arrays in different forms etc.
	// FIXME in dimensionless coordinates, this should be removed, as the grid indices are the coordinates!!
  	x_axis.resize(w->meta.grid[0]);
  	y_axis.resize(w->meta.grid[1]);
  	for(int i=0;i<w->meta.grid[0];i++){x_axis[i]=-opt.min_x+i*real(w->meta.initSpacing[0]);}
  	for(int j=0;j<w->meta.grid[1];j++){y_axis[j]=-opt.min_y+j*real(w->meta.initSpacing[1]);}

  	X = VectorXcd(w->meta.grid[0]); Y = VectorXcd(w->meta.grid[1]);
	for(int i = 0;i<w->meta.grid[0];i++){X(i) = complex<double>(x_axis[i],0.0);}
	for(int j = 0;j<w->meta.grid[1];j++){Y(j) = complex<double>(y_axis[j],0.0);}

	Xmatrix = MatrixXcd(w->meta.grid[0],w->meta.grid[1]); Ymatrix = MatrixXcd(w->meta.grid[0],w->meta.grid[1]);
	for( int i = 0; i < w->meta.grid[0]; i++){ Xmatrix.col(i) = X;	}
	for( int i = 0; i < w->meta.grid[0]; i++){ Ymatrix.row(i) = Y;	}
	// END FIXME

	PotentialGrid = MatrixXcd::Zero(w->meta.grid[0],w->meta.grid[1]);
   	for(int i = 0; i< w->meta.grid[0]; i++){
   		for(int j = 0; j < w->meta.grid[1]; j++){
			PotentialGrid(i,j) = ( half * opt.omega_x * opt.omega_x * X(i) * X(i) +  half * opt.omega_y * opt.omega_y * Y(j) * Y(j) );
		}
	}

}

void RungeKutta::assignMatrixData(MatrixData* &d) {
	 w = d;
	 setVariables();
}


void RungeKutta::timeStep(double delta_T){

	int32_t t = 0;

	computeCoefficients(delta_T);

	for(int s = 0; s < w->wavefunction.size(); s++){

		#pragma omp parallel
		{	
				#pragma omp for
				for(int i = 0; i < threads; ++i){
					wavefctcp.block(i * partx,0,partx,w->meta.grid[1]) = w->wavefunction[s].block(i * partx,0,partx,w->meta.grid[1]);
				}
	
				#pragma omp for
				for(int i = 0; i < threads; ++i){
					singleK(k0,frontx[i],endx[i],subx,suby,t);
				}
	
				#pragma omp for
				for(int i = 0; i < threads; ++i){
					wavefctcp.block(i * partx,0,partx,w->meta.grid[1]) = w->wavefunction[s].block(i * partx,0,partx,w->meta.grid[1]) + half * complex<double>(delta_T,0.0) * k0.block(i * partx,0,partx,w->meta.grid[1]);
				}
	
				#pragma omp barrier
				#pragma omp single
				{
					t++;
				}
	
				#pragma omp for
				for(int i = 0; i < threads; ++i){
					singleK(k1,frontx[i],endx[i],subx,suby,t);
				}
				#pragma omp for
				for(int i = 0; i < threads; ++i){
					wavefctcp.block(i * partx,0,partx,w->meta.grid[1]) = w->wavefunction[s].block(i * partx,0,partx,w->meta.grid[1]) + half * complex<double>(delta_T,0.0) * k1.block(i * partx,0,partx,w->meta.grid[1]);
				}
			
				#pragma omp for
				for(int i = 0; i < threads; ++i){
					singleK(k2,frontx[i],endx[i],subx,suby,t);
				}
				#pragma omp for
				for(int i = 0; i < threads; ++i){
					wavefctcp.block(i * partx,0,partx,w->meta.grid[1]) = w->wavefunction[s].block(i * partx,0,partx,w->meta.grid[1]) + complex<double>(delta_T,0.0) * k2.block(i * partx,0,partx,w->meta.grid[1]);
				}	
				
				#pragma omp barrier
				#pragma omp single
				{
					t++;
				}
			
				#pragma omp for
				for(int i = 0; i < threads; ++i){
					singleK(k3,frontx[i],endx[i],subx,suby,t);
				}
			
				#pragma omp for
				for(int i = 0; i < threads; ++i){
					wavefctcp.block(i * partx,0,partx,w->meta.grid[1]) = (one/six) * ( k0.block(i * partx,0,partx,w->meta.grid[1]) + two * ( k1.block(i * partx,0,partx,w->meta.grid[1]) + k2.block(i * partx,0,partx,w->meta.grid[1]) ) + k3.block(i * partx,0,partx,w->meta.grid[1]));
				}
	
				#pragma omp barrier
				#pragma omp single
				{
					MSDBoundaries(w->wavefunction[s],wavefctcp);
				}
	
			#pragma omp for
			for(int i = 0; i < threads; ++i){
				w->wavefunction[s].block(i * partx,0,partx,w->meta.grid[1]) += complex<double>(delta_T,0.0) * wavefctcp.block(i * partx,0,partx,w->meta.grid[1]);
			}
		}
		
		w->increment(delta_T, real(lambda_x(complex<double>(w->meta.time,0.0))), real(lambda_y(complex<double>(w->meta.time,0.0))));
	}	
}

void RungeKutta::MSDBoundaries(MatrixXcd &U,MatrixXcd &Ut){
	// #pragma omp parallel for
	for(int x = 1; x < w->meta.grid[0]-1; x++){		
		if(U(x,1).imag() != 0.0){
			Ut(x,0) = i_unit * U(x,0) * Ut(x,1).imag() / U(x,1).imag();
		}else{
			Ut(x,0) = zero;
		}
		if(U(x,w->meta.grid[1]-2).imag() != 0.0){
			Ut(x,w->meta.grid[1]-1) = i_unit * U(x,w->meta.grid[1]-1) * Ut(x,w->meta.grid[1]-2).imag() / U(x,w->meta.grid[1]-2).imag();
		}else{
			Ut(x,w->meta.grid[1]-1) = zero;
		}
	}
	// #pragma omp parallel for
	for(int y = 1; y < w->meta.grid[1]-1; y++){
		if(U(1,y).imag() != 0.0){
			Ut(0,y) = i_unit * U(0,y) * Ut(1,y).imag() / U(1,y).imag();
		}else{
			Ut(0,y) = zero;
		}
		if(U(w->meta.grid[0]-2,y).imag() != 0.0){
			Ut(w->meta.grid[0]-1,y) = i_unit * U(w->meta.grid[0]-1,y) * Ut(w->meta.grid[0]-2,y).imag() / U(w->meta.grid[0]-2,y).imag();
		}else{
			Ut(w->meta.grid[0]-1,y) = zero;
		}
	}
	Ut(0,0) = (Ut(0,1) + Ut(1,0)) / two;
	Ut(0,w->meta.grid[1]-1) = (Ut(0,w->meta.grid[1]-2) + Ut(1,w->meta.grid[1]-1)) / two;
	Ut(w->meta.grid[0]-1,0) = (Ut(w->meta.grid[0]-2,0) + Ut(w->meta.grid[0]-1,1)) / two;
	Ut(w->meta.grid[0]-1,w->meta.grid[1]-1) = (Ut(w->meta.grid[0]-2,w->meta.grid[1]-1) + Ut(w->meta.grid[0]-1,w->meta.grid[1]-2)) / two;

}

void Expansion::singleK(MatrixXcd &k, int32_t &front, int32_t &end,int32_t &subx,int32_t &suby, int &t){
	k.block(front,2,end,suby).noalias() = ( - wavefctcp.block(front-2,2,end,suby) - wavefctcp.block(front+2,2,end,suby) + sixteen * ( wavefctcp.block(front-1,2,end,suby) + wavefctcp.block(front+1,2,end,suby)) - thirty * wavefctcp.block(front  ,2,end,suby)) * laplacian_coefficient_x[t]
										+ ( - wavefctcp.block(front  ,0,end,suby) - wavefctcp.block(front  ,4,end,suby) + sixteen * ( wavefctcp.block(front  ,1,end,suby) + wavefctcp.block(front  ,3,end,suby)) - thirty * wavefctcp.block(front  ,2,end,suby)) * laplacian_coefficient_y[t];

	k.block(front,2,end,suby).array() += ( - wavefctcp.block(front+2,2,end,suby).array() + eight * ( wavefctcp.block(front+1,2,end,suby).array() - wavefctcp.block(front-1,2,end,suby).array() ) + wavefctcp.block(front-2,2,end,suby).array()) * Xmatrix.block(front,2,end,suby).array() * gradient_coefficient_x[t]
									   + ( - wavefctcp.block(front  ,4,end,suby).array() + eight * ( wavefctcp.block(front  ,3,end,suby).array() - wavefctcp.block(front  ,1,end,suby).array() ) + wavefctcp.block(front  ,0,end,suby).array()) * Ymatrix.block(front,2,end,suby).array() * gradient_coefficient_y[t];

	k.block(front,2,end,suby).array() -= i_unit * ( /*PotentialGrid.block(front,2,end,suby).array() +*/ (complex<double>(opt.g,0.0) * ( wavefctcp.block(front,2,end,suby).conjugate().array() * wavefctcp.block(front,2,end,suby).array() ))) * wavefctcp.block(front,2,end,suby).array();
}

void RotatingTrap::singleK(MatrixXcd &k, int32_t &front, int32_t &end,int32_t &subx,int32_t &suby, int &t){
	k.block(front,2,end,suby).noalias() = ( - wavefctcp.block(front-2,2,end,suby) - wavefctcp.block(front+2,2,end,suby) + sixteen * ( wavefctcp.block(front-1,2,end,suby) + wavefctcp.block(front+1,2,end,suby)) - thirty * wavefctcp.block(front  ,2,end,suby)) * laplacian_coefficient_x[t]
										+ ( - wavefctcp.block(front  ,0,end,suby) - wavefctcp.block(front  ,4,end,suby) + sixteen * ( wavefctcp.block(front  ,1,end,suby) + wavefctcp.block(front  ,3,end,suby)) - thirty * wavefctcp.block(front  ,2,end,suby)) * laplacian_coefficient_y[t];

	k.block(front,2,end,suby).array() += - (-wavefctcp.block(front+2,2,end,suby).array() + eight * ( wavefctcp.block(front+1,2,end,suby).array() - wavefctcp.block(front-1,2,end,suby).array() ) + wavefctcp.block(front-2,2,end,suby).array()) * Ymatrix.block(front,2,end,suby).array() * gradient_coefficient_y[t]
									     + (-wavefctcp.block(front  ,4,end,suby).array() + eight * ( wavefctcp.block(front  ,3,end,suby).array() - wavefctcp.block(front  ,1,end,suby).array() ) + wavefctcp.block(front  ,0,end,suby).array()) * Xmatrix.block(front,2,end,suby).array() * gradient_coefficient_x[t];

	k.block(front,2,end,suby).array() -= i_unit * ( PotentialGrid.block(front,2,end,suby).array() + (complex<double>(opt.g,0.0) * ( wavefctcp.block(front,2,end,suby).conjugate().array() * wavefctcp.block(front,2,end,suby).array() ))) * wavefctcp.block(front,2,end,suby).array();
}

void Expansion::computeCoefficients(double &delta_T){
	complex<double> absTime(w->meta.time,0.0);
	for(int32_t t = 0; t < 3; ++t){
		laplacian_coefficient_x[t] = i_unit / ( twelve * w->meta.initSpacing[0] * w->meta.initSpacing[0] * lambda_x(absTime) * lambda_x(absTime) );
   		laplacian_coefficient_y[t] = i_unit / ( twelve * w->meta.initSpacing[1] * w->meta.initSpacing[1] * lambda_y(absTime) * lambda_y(absTime) );
   		gradient_coefficient_x[t] = lambda_x_dot(absTime) / (twelve * w->meta.initSpacing[0] * lambda_x(absTime));
   		gradient_coefficient_y[t] = lambda_y_dot(absTime) / (twelve * w->meta.initSpacing[1] * lambda_y(absTime));
   		absTime += ( half * complex<double>(delta_T,0.0) );
   	}
}

void RotatingTrap::computeCoefficients(double &delta_T){
	complex<double> absTime(w->meta.time,0.0);
	for(int32_t t = 0; t < 3; ++t){
		laplacian_coefficient_x[t] = i_unit / ( twelve * w->meta.initSpacing[0] * w->meta.initSpacing[0] );
   		laplacian_coefficient_y[t] = i_unit / ( twelve * w->meta.initSpacing[1] * w->meta.initSpacing[1] );
   		gradient_coefficient_x[t] = opt.omega_w / (twelve * w->meta.initSpacing[0] );
   		gradient_coefficient_y[t] = opt.omega_w / (twelve * w->meta.initSpacing[1] );
   		absTime += ( half * complex<double>(delta_T,0.0) );
   	}
}   		


