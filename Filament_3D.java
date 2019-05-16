package tryout;

import java.util.*;



public class Filament_3D{
	
	double DomainSize;
	double time = 0;
	double Bead_Position[][];  
	double Bead_Force[][];	
	double Bead_radius[];
	double initialp[][];	//Initial positions


	int N=30;
	
	double RodLength_R=0.015;	//Simulated segment
	double PersL=10;	//Persistence length
	double sigmatheta;	//sigmatheta is used for generating the initial WLC
	double zeta=1;//4*Math.PI*0.3* RodLength_R/(0.84 + Math.log( RodLength_R/2/0.0035));;	//Drag coefficient
	double k, V, dt;	//Spring constant between the semiflexible polymers
	double kapa;	//Flexual rigidity
	double bEndGrowLength = 0;

	double pEndGrowLength = 0;
	double pEndGrowRate=0;	//Pointed end growth rate
	double temperature=300;

	double filamentlength=0.2;
	
    Random rand = new Random();
    
    //temporary vectors
    double[] x2,x1,x0;
	double[] deltas;
    double[][] basis,t;
    
    double breaktime=0;
    String severingData;


    public void start(double[][] initialpos, double theta){
    	initialp = initialpos;
	    Bead_Position=new double[N][3];
	    x2=new double[3];
	    x1=new double[3];
	    x0=new double[3];
	    t=new double[N][3];

	    
	    Bead_Position[0][0]=initialpos[0][0];
	    Bead_Position[0][1]=initialpos[0][1];
	    Bead_Position[0][2]=initialpos[0][2];

	   
	    
	    for (int count=1; count<N; count++){
	    	 Bead_Position[count][0]=Bead_Position[count-1][0]+RodLength_R*Math.cos(theta);
	 	    Bead_Position[count][1]=Bead_Position[count-1][1]+RodLength_R*Math.sin(theta);
	 	    Bead_Position[count][2]=initialpos[0][2];
	 	   
	    }
	    /*Bead_radius = new double[N];
	    
	    sigmatheta=Math.sqrt(RodLength_R/PersL);
	    
	    
	    for(int i=0;i<3;i++){
	    	t[0][i]=Bead_Position[1][i]-Bead_Position[0][i];

		    
	    }

	    
	    //this part calculates the positions for the WLC beads
	    for(int i=1;i<N-1;i++){
	    	if(i==1){		    	
	    		for(int j=0;j<3;j++){
	    		x2[j]=Bead_Position[i][j];
	    		x1[j]=Bead_Position[i-1][j];
	    	}
	    		x0[0] = x1[0];
	    		x0[1] = x1[1]-1;
	    		x0[2] = x1[2];
	    		
	    	}

	    	else{
		    	for(int j=0;j<3;j++){
		    		x2[j]=Bead_Position[i][j];
		    		x1[j]=Bead_Position[i-1][j];
		    		x0[j]=Bead_Position[i-2][j];	
		    	}
	    	}
	    	
	    	deltas = vectorNextBead(x0, x1, x2);
            
	    	Bead_Position[i+1][0]=Bead_Position[i][0]+deltas[0];
	    	Bead_Position[i+1][1]=Bead_Position[i][1]+deltas[1];
	    	Bead_Position[i+1][2]=0.2;//Bead_Position[i][2]+deltas[2];
	    	
	    	
	    	t[i][0]=deltas[0];
	    	t[i][1]=deltas[1];
	    	t[i][2]=deltas[2];
	    	

	    }
	    
	    */

    }
    	
	

	
	public void step(){		
	    Bead_Force=new double[N][3];



	
		Flow_Force_Update();
		


		move();
		Boundary_Update();
		
		
				
	    
    }


	protected void Spring_Force_Update(){
		for (int i = 0; i<N; i++){
		Bead_Force[i][2]+=Spring_Force(Bead_Position[i][2]);
	}
	}
	
	protected void Flow_Force_Update(){
		for (int i = 0; i<N; i++){
		Bead_Force[i][1]+=-zeta*V;
	}
	}
protected void Boundary_Update(){
	for (int j=0; j<2; j++){
		for (int i = 0; i<N; i++){
	
	if (Bead_Position[i][j]>DomainSize/2)
	Bead_Position[i][j]=Bead_Position[i][j]-DomainSize;
	if (Bead_Position[i][j]<-DomainSize/2)
		Bead_Position[i][j]=Bead_Position[i][j]+DomainSize;
		}
	}
}
	protected void move(){
		for (int i = 0; i<N; i++){
			for(int j=0; j<3;j++){
				Bead_Position[i][j]=Bead_Position[i][j]+Bead_Force[i][j]*dt/zeta;

			}	
		}
	}


//FORCES********************************************************
protected double Spring_Force(double z){		
	
	double f0 = 100*(z-0.2);	//harmonic restrain for the lipid bilayer.
	return f0;
} 
    private double[] vectorNextBead(double[] x02, double[] x12, double[] x22) {
    	

	    double[] tangent=new double[3];
	    double[] ntangent=new double[3];
	    double[] binormal=new double[3];
	    double[] nbinormal=new double[3];
	    double[] nnormal=new double[3];
	    double[][] basis=new double[3][3];
	    double[] changes=new double[3];
	    double[] heredeltas = new double[3];

	    
    	tangent=vectorDifference(x22, x12);
        ntangent=normalizeVector(tangent);

        binormal=crossProduct(vectorDifference(x12,x02), tangent);
        nbinormal= normalizeVector(binormal);
    
        nnormal = crossProduct(nbinormal, ntangent);
        
        basis = createBasis(ntangent, nnormal, nbinormal);
        
    	double Angle_Theta;
    	double Angle_Phi;
        
        Angle_Theta = GenerateTheta();
        Angle_Phi = GenerateRandomAngle();
        
        
        changes[0] = RodLength_R*Math.cos(Angle_Theta);
        changes[1] = RodLength_R*Math.sin(Angle_Theta)*Math.cos(Angle_Phi);
        changes[2] = RodLength_R*Math.sin(Angle_Theta)*Math.sin(Angle_Phi);
        
        heredeltas = matrixProduct(basis, changes);
        
        return heredeltas;
		
	}



	
	
	
	
	
	//CALCULATIONS**************************************************
	
    protected double GenerateTheta(){
    	double r=rand.nextDouble();
    	double randtheta=Math.sqrt(-2*sigmatheta*sigmatheta*Math.log(1-r));
    	return randtheta;
    }
    protected double GenerateRandomAngle(){
    	double r=rand.nextDouble();
    	return r*2*Math.PI;
    }
    
    //MATRIX CALCULATIONS*******************************************
    protected static double[] vectorDifference( double[] a, double[] b){
        double[] c = new double[3];
        for(int i = 0; i<3; i++)
            c[i] = a[i] - b[i];      
        return c;
    }
    protected static double[] vectorSum( double[] a, double[] b){
        double[] c = new double[3];
        for(int i = 0; i<3; i++)
            c[i] = a[i] + b[i];
        return c;
    }
    protected static double vectorproduct(double[] a, double[] b){
    	double c=0;
    	for(int i=0;i<3;i++){
    		c+=a[i]*b[i];
    	}
    	return c;
    	
    }
    protected static double[] normalizeVector(double[] a){
        double mag = magnitude(a);  
        double[] c = new double[3];
        for(int i = 0; i<3; i++)
            c[i] = a[i]/mag;
        return c;
    }
    protected static double magnitude(double[] a){
        double mag = 0;
        for(int i = 0; i<3; i++)
            mag += Math.pow(a[i],2);
        return Math.sqrt(mag);
    }
    protected static double[] crossProduct(double[] a, double[] b){
    	double[] c = new double[3];
        for(int i = 0; i<3; i++)
            c[i] = a[(i+1)%3]*b[(i+2)%3] - a[(i+2)%3]*b[(i+1)%3]; 
        return c;       
    }
    protected static double[] matrixProduct(double[][] basis, double[] a){
        double[] b = new double[3];        
        for(int i = 0; i<3; i++){           
            for(int j = 0; j<3; j++)
                b[i] += basis[i][j]*a[j];         
        }    
        return b;     
    }
    protected static double[][] createBasis(double[] t, double[] n, double[] b){
        double[][] basis = new double[3][3];  
        for(int i = 0; i<3; i++){
            basis[i][0] = t[i];
            basis[i][1] = n[i];
            basis[i][2] = b[i];
        } 
        return basis;
    }



}
