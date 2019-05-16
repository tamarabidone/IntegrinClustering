package tryout;

import java.util.Random;

public class Ligand_3D {

	double time = 0, dt;
	double Bead_PositionL[];  
	double Bead_ForceL[];	

	double initialp[];	//Initial positions

	double Length_L, k, Tension;	

	
	double zetaL;	//Drag coefficient

	double DomainSize;

	double temperature;
int numBoundIntegrins;
public int boundYN;
	
    Random rand = new Random();
    

    double zLigand;


    public void start(double[] initialpos){
    	initialp = initialpos;
	    Bead_PositionL=new double[3];
	      
	    Bead_PositionL[0]=initialpos[0];
	    Bead_PositionL[1]=initialpos[1];
	    Bead_PositionL[2]=initialpos[2];   
	    

    }
    
    /*Update forces*/
	public void step(){		
	    Bead_ForceL=new double[3];
		time+=dt;


		Thermal_Force_Update();
		Spring_Force_Update();
		


		move();
		Boundary_Update();
		
	
				
	    
    }
    
  
	

	protected void Thermal_Force_Update(){


		  	 
			double thermal_Force=Math.sqrt(2*0.0000138*temperature*zetaL/dt);
	
			Bead_ForceL[0]=thermal_Force*rand.nextGaussian();
			Bead_ForceL[1]=thermal_Force*rand.nextGaussian();
			Bead_ForceL[2]=thermal_Force*rand.nextGaussian();
			
		

	}	
	
	protected void Spring_Force_Update(){
			Bead_ForceL[2]+=Spring_Force(Bead_PositionL[2]);
		}
	protected void Boundary_Update(){
		for (int j=0; j<2; j++){
		if (Bead_PositionL[j]>DomainSize/2)
		Bead_PositionL[j]=Bead_PositionL[j]-DomainSize;
		if (Bead_PositionL[j]<-DomainSize/2)
			Bead_PositionL[j]=Bead_PositionL[j]+DomainSize;
		}
	}
	


	



	protected void move(){
	
			for(int j=0; j<3;j++){
				Bead_PositionL[j]=Bead_PositionL[j]-Bead_ForceL[j]*dt/zetaL;

			
		}
	}
	
	
	//FORCES********************************************************
	protected double Spring_Force(double z){		
		
		double f0 = 100*(z-zLigand);	//harmonic restrain for the lipid bilayer.
		return f0;
	} 
    
	
	
	
	
	
	





}
