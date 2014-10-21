#include <EXP2D_Contour.h>

using namespace std;

Contour::Contour(Options &external_opt){
	opt = external_opt;
	v_down = Vector<int32_t>(0,-1,0,opt.grid[1],opt.grid[2],opt.grid[3]);
	v_right = Vector<int32_t>(1,0,0,opt.grid[1],opt.grid[2],opt.grid[3]);
	v_up = Vector<int32_t>(0,1,0,opt.grid[1],opt.grid[2],opt.grid[3]);
	v_left = Vector<int32_t>(-1,0,0,opt.grid[1],opt.grid[2],opt.grid[3]);
}

inline void Contour::setDirection(int32_t &direction){
		switch(direction){
			case 0:
				direction = 6;
				break;
			case 1:
				direction = 0;
				break;
			case 2:
				direction = 0;
				break;
			case 3:
				direction = 2;
				break;
			case 4:
				direction = 2;
				break;
			case 5:
				direction = 4;
				break;
			case 6:
				direction = 4;
				break;
			case 7:
				direction = 6;
				break;
			default:
				cout << "Function Contour::setDirection() ERROR!" << endl;
				break;
		}
}

inline Coordinate<int32_t> Contour::nextClockwise(Coordinate<int32_t> &s, int32_t &direction){
	Coordinate<int32_t> c;

	switch(direction){
		case 0 :
			c = s + v_up;
			break;
		case 1 :
			c = s + v_right;
			break;
		case 2 :
			c = s + v_right;
			break;
		case 3 :
			c = s + v_down;
			break;
		case 4 :
			c = s + v_down;
			break;
		case 5 :
			c = s + v_left;
			break;
		case 6 :
			c = s + v_left;
			break;
		case 7 :
			c = s + v_up;
			break;
		default :
			cout << "Function Contour::nextClockwise() ERROR!" << endl;
	}
	return c; 
}

inline void Contour::findInitialP(RealGrid &data,Coordinate<int32_t> &p,Coordinate<int32_t> &s){

	Coordinate<int32_t> tmp = p;


	for(int x = tmp.x(); x < data.width(); x++){

		if(data(0,x,tmp.y(),0) > 0){

			p = data.make_coord(x,tmp.y(),0);
			s = p + v_left;
			// cout <<  p << " | " << s << endl;
			// initial[0] = p;
			// initial[1] = s;			
			break;
		}
	}
}

inline void Contour::findSecondP(RealGrid &data,Coordinate<int32_t> &p,Coordinate<int32_t> &s){

	Coordinate<int32_t> tmp = p;
	bool end = true;

	for(int x = tmp.x(); x < data.width()-1; x++){

		if(data(0,x,tmp.y(),0) > 0){

			p = data.make_coord(x,tmp.y(),0);
			s = p + v_left;
			end = true;
			// cout <<  p << " | " << s << endl;
			// initial[0] = p;
			// initial[1] = s;			
			break;
		}

	}
	if(end == false){
		p = data.make_coord(data.width()-1,tmp.y(),0);
	}
}

inline void Contour::findMostRightP(c_set &contour, Coordinate<int32_t> &p){
	c_set::iterator max_it = contour.begin();
	for(c_set::iterator it = contour.begin(); it != contour.end(); ++it){
		if(it->x() >= max_it->x())
			max_it = it;
	}
	p = *max_it;
}

c_set Contour::trackContour(RealGrid &data){

	c_set contour;
	c_set wholeContour;
	// c_set::iterator it;
	// std::pair<c_set::iterator,bool> ret;

	Coordinate<int32_t> s;
	Coordinate<int32_t> p = data.make_coord(0,data.height()/2,0);
	Coordinate<int32_t> initial[2];
	initial[0] = data.make_coord(0,0,0);
	initial[1] = data.make_coord(0,0,0);

	findInitialP(data,p,s/*,initial*/);
	contour.insert(p);

	int32_t counter = 0;
	int direction = 0;
	int insert_counter = 1;
	bool singlepoint = false;
	bool stop = false;
	// cout << "Test1: " << endl;
	// string name2 = "Test1";
	// string name3 = "Test2";
	// plotdatatopng(name3,data,opt);
	// plotContourSurround(name2,data,contour,opt);

	do{ 
		// cout << "\r" << flush;
		// cout <<  p << " | " << s << " | " << insert_counter << " " << endl;
		counter = 0;
		while(counter < 8){
			Coordinate<int32_t> c = nextClockwise(s,direction);
			if(data(0,c) > 0){
				p = c;
				setDirection(direction);
				contour.insert(p);
				insert_counter++;
				break;
			}
			
			counter++;
			if(counter == 8){
				singlepoint = true;
				break;
			}
			s = c;

			direction = (direction+1)%8;
		}

		if(singlepoint == true){
			// std::cout << "Found single point, continuing the search. " << p << endl;
			// string name = "ERROR_3-SinglePoint_" + to_string(insert_counter) + "_"+ to_string(p.x()) + "_" + to_string(p.y());
			// plotContourSurround(name, data,contour,opt);
			s = p;
			p = p + v_right;			
			contour.clear();
			findInitialP(data,p,s/*,initial*/);
			contour.insert(p);
			direction = 0;
			insert_counter = 1;
			singlepoint = false;
		}
		
		if(contour.size() >= 2){
			double scalingFromRatio;
			scalingFromRatio = (opt.omega_x.real() > opt.omega_y.real()) ? opt.omega_y.real()/opt.omega_x.real() : opt.omega_x.real()/opt.omega_y.real();
			int size_condition = (data.width()/2 - initial[0].x()) * 2 * M_PI * scalingFromRatio * 0.5; // Circumference of a circle going through p, 90%
			if((initial[0] == p) && (initial[1] == s)){
				if(contour.size() > size_condition){				
					findMostRightP(contour,p);
					s = p;
					p = p+v_right;
					findSecondP(data,p,s);
					if(p.x() == opt.grid[1]-1){
						stop = true;
					} else {
						wholeContour.insert(contour.begin(),contour.end());
						contour.clear();
						contour.insert(p);
						direction = 0;
						insert_counter =1;
					}

					// cout << "Found initial conditions with big enough contour. Size: " << contour.size() << endl;
				}
			// }else if((initial[0] == p) && (initial[1] == s)){
			// 	// cout << "Found initial conditions with small contour. Size:" << contour.size() << " Searching new contour. "<< p << " with initial " << initial[0] << " | " << initial[1] << endl;
			// 	// string name = "ERROR_2-ContourTooSmall_" + to_string(insert_counter) + "_" + to_string(p.x()) + "_" + to_string(p.y());
			// 	// plotContourSurround(name, data,contour,opt);
			// 	findMostRightP(contour,p);
			// 	s = p;
			// 	p = p + v_right;
			// 	contour.clear();
			// 	findInitialP(data,p,s/*,initial*/);
			// 	contour.insert(p);
			// 	direction = 0;
			// 	insert_counter = 1;	
			}
		}

		if(insert_counter == 2){
			initial[0] = p;
			initial[1] = s;
		}

		if(insert_counter >= (2 * contour.size() + 1)){

			findMostRightP(contour,p);
			s = p;
			p = p+v_right;
			findSecondP(data,p,s);
			if(p.x() == opt.grid[1]-1){
				stop = true;
			} else {
				// wholeContour.insert(contour.begin(),contour.end());
				contour.clear();
				contour.insert(p);
				direction = 0;
				insert_counter =1;
			}
			

			// throw std::string("Contour::trackContour error #1: Surrounded the contour two times.");



			// cout << "Surrounded the contour two times. Size:" << contour.size() << " Coordinates: " << p << " with initial " << initial[0] << " | " << initial[1] << endl;
			// string name = "ERROR_1-SurroundedTwoTimes_" + to_string(insert_counter) + "_"+ to_string(p.x()) + "_" + to_string(p.y());
			// plotContourSurround(name, data,contour,opt);

			 
			// findMostRightP(contour,p);
			// s = p;
			// p = p + v_right;
			// contour.clear();
			// findInitialP(data,p,s/*,initial*/);
			// contour.insert(p);
			// insert_counter = 1;			
		}


	}while(stop == false);
	
	wholeContour.insert(contour.begin(),contour.end());
	return wholeContour;
}