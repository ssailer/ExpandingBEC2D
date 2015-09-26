#include <contour.h>

using namespace std;

Contour::Contour(MatrixData::MetaData &external_meta){
	meta = external_meta;
	v_down = Vector<int32_t>(0,-1,0,meta.grid[0],meta.grid[1],1);
	v_right = Vector<int32_t>(1,0,0,meta.grid[0],meta.grid[1],1);
	v_up = Vector<int32_t>(0,1,0,meta.grid[0],meta.grid[1],1);
	v_left = Vector<int32_t>(-1,0,0,meta.grid[0],meta.grid[1],1);
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

inline void Contour::findInitialP(MatrixXi &data,Coordinate<int32_t> &p,Coordinate<int32_t> &s){

	Coordinate<int32_t> tmp = p;


	for(int x = tmp.x(); x < meta.grid[0]; x++){

		if(data(x,tmp.y()) > 0){

			p = Coordinate<int32_t>(x,tmp.y(),0,meta.grid[0],meta.grid[1],1);
			s = p + v_left;
			// cout <<  p << " | " << s << endl;
			// initial[0] = p;
			// initial[1] = s;			
			break;
		}
	}
	// cout << "findInitialP Size:" << " Coordinates: " << p << endl;
}

inline void Contour::findSecondP(MatrixXi &data,Coordinate<int32_t> &p,Coordinate<int32_t> &s){

	Coordinate<int32_t> tmp = p;
	bool found = false;

	for(int x = tmp.x(); x < meta.grid[0]-1; x++){

		if(data(x,tmp.y()) > 0){

			p = Coordinate<int32_t>(x,tmp.y(),0,meta.grid[0],meta.grid[1],1);
			s = p + v_left;
			found = true;
			// cout <<  p << " | " << s << endl;
			// initial[0] = p;
			// initial[1] = s;			
			break;
		}

	}
	if(found == false){
		// cout << "found end" << endl;
		p = Coordinate<int32_t>(meta.grid[0]-1,tmp.y(),0,meta.grid[0],meta.grid[1],1);
	}
	// cout << "findSecondP Size:" << " Coordinates: " << p << " density " << data(0,p.x(),p.y(),0) << endl;
}

inline void Contour::findMostRightP(c_set &contour, Coordinate<int32_t> &p){
	c_set::iterator max_it = contour.begin();
	for(c_set::iterator it = contour.begin(); it != contour.end(); ++it){
		if(it->x() >= max_it->x())
			max_it = it;
	}
	p = *max_it;
	// cout << "findMostRightP Size:" << " Coordinates: " << p << endl;
}

c_set Contour::trackContour(MatrixXi &data){

	c_set contour;
	c_set wholeContour;
	// c_set::iterator it;
	// std::pair<c_set::iterator,bool> ret;

	// double scalingFromRatio = (opt.omega_x.real() > opt.omega_y.real()) ? opt.omega_y.real()/opt.omega_x.real() : opt.omega_x.real()/opt.omega_y.real();

	Coordinate<int32_t> s;
	Coordinate<int32_t> p = Coordinate<int32_t>(0,meta.grid[1]/2,0,meta.grid[0],meta.grid[1],1);
	Coordinate<int32_t> initial[2];
	initial[0] = Coordinate<int32_t>(0,0,0,meta.grid[0],meta.grid[1],1);
	initial[1] = Coordinate<int32_t>(0,0,0,meta.grid[0],meta.grid[1],1);

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
			if(data(c.x(),c.y()) > 0){
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

			findMostRightP(contour,p);
			s = p;
			p = p+v_right;

			findSecondP(data,p,s);
			if(p.x() == meta.grid[0]-1){
				stop = true;
			} else {
				contour.clear();
				contour.insert(p);
				direction = 0;
				insert_counter =1;
			}


			// s = p;
			// p = p + v_right;			
			// contour.clear();
			// findInitialP(data,p,s/*,initial*/);
			// contour.insert(p);
			// direction = 0;
			// insert_counter = 1;
			// singlepoint = false;
		}
		
		if(contour.size() >= 100){
			
			// int size_condition = (data.width()/2 - initial[0].x()) * 2 * M_PI * scalingFromRatio * 0.5; // Circumference of a circle going through p, 90%
			if((initial[0] == p) && (initial[1] == s)){
				// if(contour.size() > size_condition){				

					wholeContour.insert(contour.begin(),contour.end());
					// cout << "Contour saved to wholeContour." << endl;

					findMostRightP(contour,p);
					s = p;
					p = p+v_right;
					findSecondP(data,p,s);
					if(p.x() == meta.grid[0]-1){
						stop = true;
					} else {
						contour.clear();
						contour.insert(p);
						direction = 0;
						insert_counter =1;
					}

					// cout << "Found initial conditions with big enough contour. Size: " << contour.size() << endl;
				// }
				// cout << "Found initial conditions with too small contour. Size: " << contour.size() << endl;
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

			// cout << "Surrounded the contour two times. Size:" << contour.size() << " Coordinates: " << p << " with initial " << initial[0] << " | " << initial[1] << endl;
			// string name = "ERROR_1-SurroundedTwoTimes_" + to_string(insert_counter) + "_"+ to_string(p.x()) + "_" + to_string(p.y());
			// plotContourSurround(name, data,contour,opt);


			findMostRightP(contour,p);
			s = p;
			p = p+v_right;

			findSecondP(data,p,s);
			if(p.x() == meta.grid[0]-1){
				stop = true;
			} else {
				// wholeContour.insert(contour.begin(),contour.end());
				// cout << "After surrounding the contour two times. Size:" << contour.size() << " Coordinates: " << p << " with initial " << initial[0] << " | " << initial[1] << endl;
				// cout << "found new coordinates." << endl;
				contour.clear();
				contour.insert(p);
				direction = 0;
				insert_counter =1;
			}
			

			// throw std::string("Contour::trackContour error #1: Surrounded the contour two times.");





			 
			// findMostRightP(contour,p);
			// s = p;
			// p = p + v_right;
			// contour.clear();
			// findInitialP(data,p,s/*,initial*/);
			// contour.insert(p);
			// insert_counter = 1;			
		}
	if(insert_counter >= meta.grid[0] * meta.grid[1]){
		stop = true;
	}

	}while(stop == false);
	
	wholeContour.insert(contour.begin(),contour.end());
	smooth(wholeContour);
	return wholeContour;
}

void Contour::smooth(c_set &c){
	
}