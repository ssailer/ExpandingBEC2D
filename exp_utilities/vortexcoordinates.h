//*****Vortex Initial Coordinates*****

//Spacing (6% across and 3% up relative to the lattice grid size)
const int across=6,up=3;
int half_across=across/2;

//Central Vortex (at the origin)
int x_1=opt.grid[1]/2;
int y_1=opt.grid[2]/2; 

//Ring 1 
int x_2=(100-across)*opt.grid[1]/200,x_3=(100+across)*opt.grid[1]/200,x_4=(100+half_across)*opt.grid[1]/200,x_5=(100-half_across)*opt.grid[1]/200,x_6=(100-half_across)*opt.grid[1]/200,x_7=(100+half_across)*opt.grid[1]/200;
int y_2=opt.grid[2]/2,y_3=opt.grid[2]/2,y_4=(100+up)*opt.grid[2]/200,y_5=(100+up)*opt.grid[2]/200,y_6=(100-up)*opt.grid[2]/200,y_7=(100-up)*opt.grid[2]/200; 

//Ring 2
int x_8=(100-2*across)*opt.grid[1]/200,x_9=(100-across-half_across)*opt.grid[1]/200,x_10=(100-across)*opt.grid[1]/200,x_11=opt.grid[1]/2,x_12=(100+across)*opt.grid[1]/200,x_13=(100+across+half_across)*opt.grid[1]/200,x_14=(100+2*across)*opt.grid[1]/200,x_15=(100+across+half_across)*opt.grid[1]/200,x_16=(100+across)*opt.grid[1]/200,x_17=opt.grid[1]/2,x_18=(100-across)*opt.grid[1]/200,x_19=(100-across-half_across)*opt.grid[1]/200;
int y_8=opt.grid[2]/2,y_9=(100+up)*opt.grid[2]/200,y_10=(100+2*up)*opt.grid[2]/200,y_11=(100+2*up)*opt.grid[2]/200,y_12=(100+2*up)*opt.grid[1]/200,y_13=(100+up)*opt.grid[2]/200,y_14=opt.grid[2]/2,y_15=(100-up)*opt.grid[2]/200,y_16=(100-2*up)*opt.grid[2]/200,y_17=(100-2*up)*opt.grid[2]/200,y_18=(100-2*up)*opt.grid[2]/200,y_19=(100-up)*opt.grid[2]/200;

//Ring 3
int x_20=(100-3*across)*opt.grid[1]/200,x_21=(100-2*across-half_across)*opt.grid[1]/200,x_22=(100-2*across)*opt.grid[1]/200,x_23=(100-across-half_across)*opt.grid[1]/200,x_24=(100-half_across)*opt.grid[1]/200,x_25=(100+half_across)*opt.grid[1]/200,x_26=(100+across+half_across)*opt.grid[1]/200,x_27=(100+2*across)*opt.grid[1]/200,x_28=(100+2*across+half_across)*opt.grid[1]/200,x_29=(100+3*across)*opt.grid[1]/200,x_30=(100+2*across+half_across)*opt.grid[1]/200,x_31=(100+2*across)*opt.grid[1]/200,x_32=(100+across+half_across)*opt.grid[1]/200,x_33=(100+half_across)*opt.grid[1]/200,x_34=(100-half_across)*opt.grid[1]/200,x_35=(100-across-half_across)*opt.grid[1]/200,x_36=(100-2*across)*opt.grid[1]/200,x_37=(100-2*across-half_across)*opt.grid[1]/200;
int y_20=opt.grid[2]/2,y_21=(100+up)*opt.grid[2]/200,y_22=(100+2*up)*opt.grid[2]/200,y_23=(100+3*up)*opt.grid[2]/200,y_24=(100+3*up)*opt.grid[2]/200,y_25=(100+3*up)*opt.grid[2]/200,y_26=(100+3*up)*opt.grid[2]/200,y_27=(100+2*up)*opt.grid[2]/200,y_28=(100+up)*opt.grid[2]/200,y_29=opt.grid[2]/2,y_30=(100-up)*opt.grid[2]/200,y_31=(100-2*up)*opt.grid[2]/200,y_32=(100-3*up)*opt.grid[2]/200,y_33=(100-3*up)*opt.grid[2]/200,y_34=(100-3*up)*opt.grid[2]/200,y_35=(100-3*up)*opt.grid[2]/200,y_36=(100-2*up)*opt.grid[2]/200,y_37=(100-up)*opt.grid[2]/200;  

//Ring 4
int x_38=(100-4*across)*opt.grid[1]/200,x_39=(100-3*across-half_across)*opt.grid[1]/200,x_40=(100-3*across)*opt.grid[1]/200,x_41=(100-2*across-half_across)*opt.grid[1]/200,x_42=(100-2*across)*opt.grid[1]/200,x_43=(100-across)*opt.grid[1]/200,x_44=opt.grid[1]/2,x_45=(100+across)*opt.grid[1]/200,x_46=(100+2*across)*opt.grid[1]/200,x_47=(100+2*across+half_across)*opt.grid[1]/200,x_48=(100+3*across)*opt.grid[1]/200,x_49=(100+3*across+half_across)*opt.grid[1]/200,x_50=(100+4*across)*opt.grid[1]/200,x_51=(100+3*across+half_across)*opt.grid[1]/200,x_52=(100+3*across)*opt.grid[1]/200,x_53=(100+2*across+half_across)*opt.grid[1]/200,x_54=(100+2*across)*opt.grid[1]/200,x_55=(100+across)*opt.grid[1]/200,x_56=opt.grid[1]/2,x_57=(100-across)*opt.grid[1]/200,x_58=(100-2*across)*opt.grid[1]/200,x_59=(100-2*across-half_across)*opt.grid[1]/200,x_60=(100-3*across)*opt.grid[1]/200,x_61=(100-3*across-half_across)*opt.grid[1]/200;
int y_38=opt.grid[2]/2,y_39=(100+up)*opt.grid[2]/200,y_40=(100+2*up)*opt.grid[2]/200,y_41=(100+3*up)*opt.grid[2]/200,y_42=(100+4*up)*opt.grid[2]/200,y_43=(100+4*up)*opt.grid[2]/200,y_44=(100+4*up)*opt.grid[2]/200,y_45=(100+4*up)*opt.grid[2]/200,y_46=(100+4*up)*opt.grid[2]/200,y_47=(100+3*up)*opt.grid[2]/200,y_48=(100+2*up)*opt.grid[2]/200,y_49=(100+up)*opt.grid[2]/200,y_50=opt.grid[2]/2,y_51=(100-up)*opt.grid[2]/200,y_52=(100-2*up)*opt.grid[2]/200,y_53=(100-3*up)*opt.grid[2]/200,y_54=(100-4*up)*opt.grid[2]/200,y_55=(100-4*up)*opt.grid[2]/200,y_56=(100-4*up)*opt.grid[2]/200,y_57=(100-4*up)*opt.grid[2]/200,y_58=(100-4*up)*opt.grid[2]/200,y_59=(100-3*up)*opt.grid[2]/200,y_60=(100-2*up)*opt.grid[2]/200,y_61=(100-up)*opt.grid[2]/200;


 double vortex(int a, int b, int x, int y) //Vortex with phase [0,2*pi)          
{
        if(atan2(b-y,a-x)<0){ return 2*pi+atan2(b-y,a-x); } //atan2 is defined from [-pi,pi) so it needs to be changed to [0,2*pi)
	else{ return atan2(b-y,a-x); }        
}

void add_vortex() // computes the phasefield by adding all vortices together and saves the final psi to use for timeevolution
{	
	ComplexGrid Phase = ComplexGrid(opt.grid[0],opt.grid[1],opt.grid[2],opt.grid[3]);
	ComplexGrid c = ComplexGrid(opt.grid[0],opt.grid[2],opt.grid[2],opt.grid[3]);

   	for(i=0;i<opt.grid[1];i++)
 	{
      		for(j=0;j<opt.grid[2];j++)
       		{
 		        Phase(0,i,j,0) = run->phase_save(real(run->pPsi),i,j)+vortex(i,j,x_1,y_1)
			  /*Ring 1*/
				+vortex(i,j,x_2,y_2)+vortex(i,j,x_3,y_3)+vortex(i,j,x_4,y_4)+vortex(i,j,x_5,y_5)+vortex(i,j,x_6,y_6)+vortex(i,j,x_7,y_7)
 			  /*Ring 2*/
 				+vortex(i,j,x_8,y_8)+vortex(i,j,x_9,y_9)+vortex(i,j,x_10,y_10)+vortex(i,j,x_11,y_11)+vortex(i,j,x_12,y_12)+vortex(i,j,x_13,y_13)+vortex(i,j,x_14,y_14)+vortex(i,j,x_15,y_15)+vortex(i,j,x_16,y_16)+vortex(i,j,x_17,y_17)+vortex(i,j,x_18,y_18)+vortex(i,j,x_19,y_19)
 			  /*Ring 3*/
 				+vortex(i,j,x_20,y_20)+vortex(i,j,x_21,y_21)+vortex(i,j,x_22,y_22)+vortex(i,j,x_23,y_23)+vortex(i,j,x_24,y_24)+vortex(i,j,x_25,y_25)+vortex(i,j,x_26,y_26)+vortex(i,j,x_27,y_27)+vortex(i,j,x_28,y_28)+vortex(i,j,x_29,y_29)+vortex(i,j,x_30,y_30)+vortex(i,j,x_31,y_31)+vortex(i,j,x_32,y_32)+vortex(i,j,x_33,y_33)+vortex(i,j,x_34,y_34)+vortex(i,j,x_35,y_35)+vortex(i,j,x_36,y_36)+vortex(i,j,x_37,y_37) 
 			  /*Ring 4*/
 			  	+vortex(i,j,x_38,y_38)+vortex(i,j,x_39,y_39)+vortex(i,j,x_40,y_40)+vortex(i,j,x_41,y_41)+vortex(i,j,x_42,y_42)+vortex(i,j,x_43,y_43)+vortex(i,j,x_44,y_44)+vortex(i,j,x_45,y_45)+vortex(i,j,x_46,y_46)+vortex(i,j,x_47,y_47)+vortex(i,j,x_48,y_48)+vortex(i,j,x_49,y_49)+vortex(i,j,x_50,y_50)+vortex(i,j,x_51,y_51)+vortex(i,j,x_52,y_52)+vortex(i,j,x_53,y_53)+vortex(i,j,x_54,y_54)+vortex(i,j,x_55,y_55)+vortex(i,j,x_56,y_56)+vortex(i,j,x_57,y_57)+vortex(i,j,x_58,y_58)+vortex(i,j,x_59,y_59)+vortex(i,j,x_60,y_60)+vortex(i,j,x_61,y_61);

 			c(0,i,j,0)=polar(abs(run->pPsi->at(0,i,j,0)),Phase(0,i,j,0)); // compute psi by using the initial psi^2 and adding the phase 

       		}
    	}

    for(i=0;i<opt.grid[1];i++){for(j=0;j<opt.grid[2];j++){ run->pPsi->at(0,i,j,0) = c(0,i,j,0); } }
 }
