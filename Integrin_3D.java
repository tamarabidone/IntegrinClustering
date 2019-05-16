package tryout;

import java.util.*;



public class Integrin_3D{
	
	
	double time = 0, dt;
	double Bead_Position[];  
	double Bead_Force[];	
double v, cut_off, ka, kd, affinityI, affinityL;
public int active;
public int clustered;
int ligandBoundYN, boundB, boundF;
public int ligandNum;
	double initialp[];	//Initial positions

	double zetaI, zetaI_initial;	//Drag coefficient

	double DomainSize;

	double temperature;

	
    Random rand = new Random();
    

    double zIntegrin;


    public void start(double[] initialpos){
    	initialp = initialpos;
	    Bead_Position=new double[3];
	      
	    Bead_Position[0]=initialpos[0];
	    Bead_Position[1]=initialpos[1];
	    Bead_Position[2]=initialpos[2];   
	    

    }
    

    
    /*Update forces*/
	public void step(){		
	    Bead_Force=new double[3];
		time+=dt;


		Thermal_Force_Update();
		Spring_Force_Update();
		


		move();
		Boundary_Update();
		
		clustered=0;
				
	    
    }
    
	
	
	
	

	protected void Thermal_Force_Update(){


		  	 
			double thermal_Force=Math.sqrt(2*0.0000138*temperature*zetaI/dt);
	
			Bead_Force[0]=thermal_Force*rand.nextGaussian();
			Bead_Force[1]=thermal_Force*rand.nextGaussian();
			Bead_Force[2]=thermal_Force*rand.nextGaussian();
			
		

	}	
	
	protected void Spring_Force_Update(){
			Bead_Force[2]+=Spring_Force(Bead_Position[2]);
		}
	protected void Boundary_Update(){
		for (int j=0; j<2; j++){
		if (Bead_Position[j]>DomainSize/2)
		Bead_Position[j]=Bead_Position[j]-DomainSize;
		if (Bead_Position[j]<-DomainSize/2)
			Bead_Position[j]=Bead_Position[j]+DomainSize;
		}
	}
	


	



	protected void move(){
	
			for(int j=0; j<3;j++){
				Bead_Position[j]=Bead_Position[j]-Bead_Force[j]*dt/zetaI;

			
		}
	}
	
	
	//FORCES********************************************************
	protected double Spring_Force(double z){		
		
		double f0 = 100*(z-zIntegrin);	//harmonic restrain for the lipid bilayer.
		return f0;
	}	
	
   
}
