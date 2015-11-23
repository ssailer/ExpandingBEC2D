	void Bh3Evaluation::assert_space(Space s)
{
	if(space != s)
	{
		for(int i = 0; i < data.size(); i++)
			ComplexGrid::fft(data[i], data[i], s == KSpace);
		
		space = s;
	}
}

		typedef struct {
			int32_t n;
			Coordinate<double> x;
			vector<Vector<double> > velocity;
			list<Coordinate<int32_t> > points;
			int32_t num_points;
			double pair_distance;
		} VortexData;

		vector<vector<double> > kspace;
		vector<vector<complex<double> > > kfspace, kbspace;

		enum Space { RSpace, KSpace };

		struct PathResults {
			list<VortexData> vlist;
		};

	class Averages {
		public:
			double nm3, nm2, nm1, n0, n1, n2, n3, mass_zeros, num_vortices;
			complex<double> meanfield;
			ArrayXXf vortex_velocities_x, vortex_velocities_y, vortex_velocities_norm;
	               ArrayXd velocities_x, velocities_y, velocities_z, velocities_norm;
			ArrayXd av_vortex_velocity;
			double av_velocity;
			double pair_distance_all, pair_distance_nonzero, max_pair_distance;
			ArrayXf /*pd_histogram_all, pd_histogram_closest,*/ g2, g2_av, g2_vv, g2_closest;
			ArrayXcd g1;
			double kinetick_int, pressure_int, interaction_int, particle_count, ckinetic_int, ikinetic_int,
					j_int, cj_int, ij_int, Ekin;
			ArrayXd number, omega, ikinetick, ckinetick, kinetick, pressure, j, cj, ij, pressure_kin_mixture, phase;
			ArrayXd ikinetick_wo_phase, ckinetick_wo_phase, kinetick_wo_phase, pressure_wo_phase;
			ArrayXd k;
			double time;
			
			Averages() {}
			Averages(int avgrid, int deltas);
			Averages operator+ (const Averages &a) const;
			Averages &operator+= (const Averages &a);
			Averages operator- (const Averages &a) const;
			Averages &operator-= (const Averages &a);
			Averages operator* (const Averages &a) const;
			Averages &operator*= (const Averages &a);
			Averages operator/ (double d) const;
			Averages &operator/= (double d);
			Averages operator* (double d) const;
			Averages &operator*= (double d);
		};
inline Bh3Evaluation::Averages::Averages(int avgrid, int deltas) :
		vortex_velocities_x(deltas, 4*avgrid),
		vortex_velocities_y(deltas, 4*avgrid),
		vortex_velocities_norm(deltas, 4*avgrid),
		av_vortex_velocity(deltas),
		velocities_x(avgrid),
		velocities_y(avgrid),
		velocities_z(avgrid),
		velocities_norm(avgrid),
		g1(avgrid),
		g2(avgrid),
		g2_av(avgrid),
		g2_vv(avgrid),
		g2_closest(avgrid),
		number(avgrid),
		ikinetick(avgrid),
		ckinetick(avgrid),
		kinetick(avgrid),
		pressure(avgrid),
		pressure_kin_mixture(avgrid),
		j(avgrid),
		cj(avgrid),
		ij(avgrid),
		omega(avgrid),
		k(avgrid),
		ikinetick_wo_phase(avgrid),
		ckinetick_wo_phase(avgrid),
		kinetick_wo_phase(avgrid),
		pressure_wo_phase(avgrid),
		phase(avgrid)
{
	time = nm3 = nm2 = nm1 = n0 = n1 = n2 = n3 = mass_zeros = num_vortices = 0.0;
	pair_distance_all = pair_distance_nonzero = max_pair_distance = 0.0;
	Ekin = kinetick_int = pressure_int = interaction_int = particle_count = ckinetic_int = ikinetic_int = j_int = cj_int = ij_int = 0.0;
	meanfield = 0.0;
	av_velocity = 0.0;

    vortex_velocities_x.setZero();
    vortex_velocities_y.setZero();
    vortex_velocities_norm.setZero();
    av_vortex_velocity.setZero();
    velocities_x.setZero();
    velocities_y.setZero();
    velocities_z.setZero();
    velocities_norm.setZero();
    g1.setZero();
    g2.setZero();
    g2_av.setZero();
    g2_vv.setZero();
    g2_closest.setZero();
    number.setZero();
    ikinetick.setZero();
    ckinetick.setZero();
    kinetick.setZero();
    pressure.setZero();
    pressure_kin_mixture.setZero();
    j.setZero();
    cj.setZero();
    ij.setZero();
    omega.setZero();
    k.setZero();
    ikinetick_wo_phase.setZero();
    ckinetick_wo_phase.setZero();
    kinetick_wo_phase.setZero();
    pressure_wo_phase.setZero();
    phase.setZero();
}

double kwidth2[3];
	for(int i = 0; i < 3; i++)
		kwidth2[i] = (options.grid[i+1] == 1) ? 0 : options.klength[i]*options.klength[i];
double index_factor = (ares.number.size() - 1) / sqrt(kwidth2[0] + kwidth2[1] + kwidth2[2]);

		PathResults pres;
		Averages ares;
		Space space;


void Bh3Evaluation::init_k_space()
{
	kspace.resize(3);
	kfspace.resize(3);
	kbspace.resize(3);
	for(int d = 0; d < 3; d++)
	{
		// set k-space
		kspace[d].resize(options.grid[d+1]);
		for (int i=0; i<options.grid[d+1]/2; i++)
			kspace[d][i] = options.klength[d]*sin( M_PI*((double)i)/((double)options.grid[d+1]) );

		for (int i=options.grid[d]/2; i<options.grid[d+1]; i++)
			kspace[d][i] = options.klength[d]*sin( M_PI*((double)(-options.grid[d+1]+i))/((double)options.grid[d+1]) );

		//Set k-forward
		kfspace[d].resize(options.grid[d+1]);
		for (int i=0; i<options.grid[d+1]; i++)
			kfspace[d][i] = - complex<double>(0,1) * (polar( 1.0, 2*M_PI*((double)i)/((double)options.grid[d+1]) ) - 1.0);
		//Set k-backward
		kbspace[d].resize(options.grid[d+1]);
		for (int i=0; i<options.grid[d+1]; i++)
			kbspace[d][i] = complex<double>(0,1) * (polar( 1.0, - 2*M_PI*((double)i)/((double)options.grid[d+1]) ) - 1.0);
	}
}


complex<double> Bh3Evaluation::calc_compressible_part_x(int x, int y, int z, complex<double> kx, complex<double> ky, complex<double> kz) const
{
	return (kx*kfspace[0][x] + ky*kfspace[1][y] + kz*kfspace[2][z]) * kbspace[0][x] /
			(kbspace[0][x]*kfspace[0][x] + kbspace[1][y]*kfspace[1][y] + kbspace[2][z]*kfspace[2][z]);
}

	// temporaere Variable
	vector<ComplexGrid> field(3, ComplexGrid(1,options.grid[1], options.grid[2], options.grid[3]));
	//Vector<int32_t> up = data->make_vector(0,1,0);
	Vector<int32_t> down = data[0].make_vector(0,-1,0);
	//Vector<int32_t> right = data->make_vector(1,0,0);
	Vector<int32_t> left = data[0].make_vector(-1,0,0);
	//Vector<int32_t> forward = data->make_vector(0,0,1);
	Vector<int32_t> back = data[0].make_vector(0,0,-1);

	inline Coordinate<int32_t> make_coord(int x, int y, int z) const {return Coordinate<int32_t>(x,y,z,width(),height(),depth());}
	inline Vector<int32_t> make_vector(int x, int y, int z) const {return Vector<int32_t>(x,y,z,width(),height(),depth());}


	// Im Ortsraum Geschwindigkeitsfeld berechnen
	#pragma omp parallel for schedule(guided, 1)
	for(int x = 0; x < data[0].width(); x++)
	{
		for (int y = 0; y < data[0].height(); y++)
		{
			for (int z = 0; z < data[0].depth(); z++)
			{
				Coordinate<int32_t> c = data[0].make_coord(x,y,z);
				
				double sqrt_n_back[3] = {abs(data[0](0,c+left)),
										 abs(data[0](0,c+down)),
										 abs(data[0](0,c+back))};
				double sqrt_n = abs(data[0](0,c));
				complex<double> exp_i_phi_back[3] = {data[0](0,c+left) / sqrt_n_back[0],
                                                     data[0](0,c+down) / sqrt_n_back[1],
                                                     data[0](0,c+back) / sqrt_n_back[2]};
				complex<double> exp_i_phi = data[0](0,c) / sqrt_n;
				
				field[0](0,c) = (sqrt_n_back[0] + sqrt_n) / 2.0 * (exp_i_phi - exp_i_phi_back[0]) / exp_i_phi;
				field[1](0,c) = (sqrt_n_back[1] + sqrt_n) / 2.0 * (exp_i_phi - exp_i_phi_back[1]) / exp_i_phi;
				field[2](0,c) = (sqrt_n_back[2] + sqrt_n) / 2.0 * (exp_i_phi - exp_i_phi_back[2]) / exp_i_phi;
			}
		}
	}

	ComplexGrid::fft(field[0],field[0]);
	ComplexGrid::fft(field[1],field[1]);
	ComplexGrid::fft(field[2],field[2]);
	assert_space(KSpace);
	
	// Im Impulsraum Werte aufaddieren
	for(int x = 0; x < data[0].width(); x++)
	{
		for (int y = 0; y < data[0].height(); y++)
		{
			for (int z = 0; z < data[0].depth(); z++)
			{
				double k = sqrt(kspace[0][x]*kspace[0][x] + kspace[1][y]*kspace[1][y] + kspace[2][z]*kspace[2][z]);

				int index = index_factor * k;
				Coordinate<int32_t> c = data[0].make_coord(x,y,z);
				
				// occupation number
				ares.k(index) += k;
				double number = abs2(data[0](0,c));
				ares.number(index) += number;
				ares.particle_count += number;
				ares.Ekin += number * k * k;
				
				if(data.size() > 1)
					ares.omega(index) += abs((data[0](0,c) - data[1](0,c)) / data[0](0,c));
				
				// Calculating the incompressible:
				// compressible and incompressible velocity:
				complex<double> cvelofield[3];
				if((x==0)&&(y==0)&&(z==0))
				{
					cvelofield[0] = field[0](0,0,0,0);
					cvelofield[1] = field[1](0,0,0,0);
					cvelofield[2] = field[2](0,0,0,0);
				}
				else
				{
					cvelofield[0] = calc_compressible_part_x(x,y,z,field[0](0,c),field[1](0,c), field[2](0,c));
					cvelofield[1] = calc_compressible_part_y(x,y,z,field[0](0,c),field[1](0,c), field[2](0,c));
					cvelofield[2] = calc_compressible_part_z(x,y,z,field[0](0,c),field[1](0,c), field[2](0,c));
				}
				
				double ckinetick = abs2(cvelofield[0]) + abs2(cvelofield[1]) + abs2(cvelofield[2]);
				double kinetick = abs2(field[0](0,c)) + abs2(field[1](0,c)) + abs2(field[2](0,c));
				double ikinetick = kinetick - ckinetick;
				
				ares.ikinetick_wo_phase(index) += ikinetick;
				ares.ckinetick_wo_phase(index) += ckinetick;
				ares.kinetick_wo_phase(index) += kinetick;
				
				ares.ikinetic_int += ikinetick;
				ares.ckinetic_int += ckinetick;
				ares.kinetick_int += kinetick;
				
				field[0](0,c) = cvelofield[0];
				field[1](0,c) = cvelofield[1];
				field[2](0,c) = cvelofield[2];
			}
		}
	}