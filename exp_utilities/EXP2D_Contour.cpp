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

inline void Contour::findInitialP(RealGrid &data,Coordinate<int32_t> &p,Coordinate<int32_t> &s, Coordinate<int32_t> *initial){

	Coordinate<int32_t> tmp = p;


	for(int x = tmp.x(); x < data.width(); x++){

		if(data(0,x,tmp.y(),0) > 0){
			p = data.make_coord(x,tmp.y(),0);
			s = p + v_left;
			initial[0] = p;
			initial[1] = s;			
			break;
		}
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
	// c_set::iterator it;
	// std::pair<c_set::iterator,bool> ret;

	Coordinate<int32_t> s;
	Coordinate<int32_t> p = data.make_coord(0,data.height()/2,0);
	Coordinate<int32_t> initial[2];

	findInitialP(data,p,s,initial);
	contour.insert(p);
	
	int32_t counter = 0;
	int direction = 0;
	int insert_counter = 0;
	bool singlepoint = false;
	bool stop = false;
	cout << endl;
	cout << "Test1: " << endl;
	string name2 = "Test1";
	// string name3 = "Test2";
	// plotdatatopng(name3,data,opt);
	plotContourSurround(name2,data,contour,opt);

	do{
		cout << "\r" << flush;
		cout << initial[0] << " | " << initial[1] << " | " << p << " | " << s << " | " << insert_counter << " ";
		
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
		counter = 0;

		if(singlepoint == true){
			cout << endl << "Found single point, continuing the search." << endl;
			string name = "singlepoint" + to_string(insert_counter) + "_"+ to_string(p.x()) + "_" + to_string(p.y());
			plotContourSurround(name, data,contour,opt);
			contour.clear();			
			s = p;
			p = p + v_right;			
			findInitialP(data,p,s,initial);
			contour.insert(p);
			insert_counter = 1;
			singlepoint = false;
			break;

		}

		if(insert_counter >= 2 * contour.size()){
			cout << "Surrounded the contour two times. Searching new contour." << endl;
			string name = "twotimes" + to_string(insert_counter) + "_"+ to_string(p.x()) + "_" + to_string(p.y());
			plotContourSurround(name, data,contour,opt);
			findMostRightP(contour,p);
			s = p;
			p = p + v_right;
			contour.clear();
			findInitialP(data,p,s,initial);
			contour.insert(p);
			insert_counter = 1;			
		}

		int size_condition = (data.width()/2 - p.x()) * 2 * M_PI * 0.9; // Circumference of a circle going through p, 90%
		if(insert_counter > size_condition){
			if((initial[0] == p) && (initial[1] == s)){
				cout << endl << "Found initial conditions with big enough contour. Size: " << contour.size() << endl;
				stop = true;
			}
		}else if((initial[0] == p) && (initial[1] == s)){
			cout << endl << "Found initial conditions with small contour. Size:" << contour.size() << " Searching new contour." << endl;
			string name = "toosmall" + to_string(insert_counter) + "_"+ to_string(p.x()) + "_"+ to_string(p.y());
			plotContourSurround(name, data,contour,opt);
			findMostRightP(contour,p);
			s = p;
			p = p + v_right;
			contour.clear();
			findInitialP(data,p,s,initial);
			contour.insert(p);
			insert_counter = 1;	
		}

	}while(stop == false);
	

	return contour;
}