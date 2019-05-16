package tryout;

import java.util.ArrayList;

public class PairwiseInteractions {
	
	public void TruncatedHarmonic (ArrayList <Integrin_3D> sIntegrin,double cut_off){
		

		 for (int i1=0; i1<sIntegrin.size()-1; i1++){
				if (sIntegrin.get(i1).active==1){
		 	for (int i2=i1+1; i2<sIntegrin.size(); i2++){
		 		if (sIntegrin.get(i2).active==1)
		 		{
		 		double dist=Distance(sIntegrin.get(i1), sIntegrin.get(i2));
		 		double k=200;
		 		double x=dist-0.002;
		 		double force=k*x;
		 		
		 		if (dist<sIntegrin.get(i1).cut_off){
		 			sIntegrin.get(i1).clustered=1;
		 			sIntegrin.get(i2).clustered=1;
		 			
		 			sIntegrin.get(i2).Bead_Position=move(sIntegrin.get(i2).Bead_Position, sIntegrin.get(i1).Bead_Position, force, dist, sIntegrin.get(i2).dt, sIntegrin.get(i2).zetaI);
					sIntegrin.get(i1).Bead_Position=move(sIntegrin.get(i1).Bead_Position, sIntegrin.get(i2).Bead_Position, force, dist, sIntegrin.get(i1).dt, sIntegrin.get(i1).zetaI);											
		 		
					sIntegrin.get(i1).zetaI =  400;//4.114*1E-3/0.00004;
		 			sIntegrin.get(i2).zetaI =  400;//4.114*1E-3/0.00004;
		 			
		 		}		
		 		}		
		 		}
		 	}	
		 }
		
	}
	
	public void Gaussian (ArrayList <Integrin_3D> sIntegrin,double cut_off){
		 for (int i1=0; i1<sIntegrin.size()-1; i1++){
				if (sIntegrin.get(i1).active==1){
		 	for (int i2=i1+1; i2<sIntegrin.size(); i2++){
		 		if (sIntegrin.get(i2).active==1)
		 		{
		 		double dist=Distance(sIntegrin.get(i1), sIntegrin.get(i2));
		 		double v0=0.1;//0.01;
		 		double sigma=0.011;
		 		double x=dist-0.002;
		 		double force=(-v0*Math.exp(-(x*x)/(2*sigma*sigma))*(-x/(sigma*sigma)));
		 		
		 		if (dist<sIntegrin.get(i1).cut_off){
		 			sIntegrin.get(i1).clustered=1;
		 			sIntegrin.get(i2).clustered=1;

		 			
		 			sIntegrin.get(i2).Bead_Position=move(sIntegrin.get(i2).Bead_Position, sIntegrin.get(i1).Bead_Position, force, dist, sIntegrin.get(i2).dt, sIntegrin.get(i2).zetaI);
					sIntegrin.get(i1).Bead_Position=move(sIntegrin.get(i1).Bead_Position, sIntegrin.get(i2).Bead_Position, force, dist, sIntegrin.get(i1).dt, sIntegrin.get(i1).zetaI);											
		 			sIntegrin.get(i1).zetaI =  400;//4.114*1E-3/0.04;
		 			sIntegrin.get(i2).zetaI =  400;//4.114*1E-3/0.04;		
		 		}
		 		}
		 	}	
		 	}
		 }
		
	}
	
	public void TruncatedHarmonicShifted (ArrayList <Integrin_3D> sIntegrin,double cut_off){
		

		 for (int i1=0; i1<sIntegrin.size()-1; i1++){
			 if (sIntegrin.get(i1).active==1){
		 	for (int i2=i1+1; i2<sIntegrin.size(); i2++){
		 		if (sIntegrin.get(i2).active==1)
		 		{
		 		double dist=Distance(sIntegrin.get(i1), sIntegrin.get(i2));
		 		double k=1;
		 		double x=dist-0.005;
		 		double force=k*x-k*(sIntegrin.get(i1).cut_off-0.01);		 		
		 		if ( dist<sIntegrin.get(i1).cut_off){
		 			sIntegrin.get(i1).clustered=1;
		 			sIntegrin.get(i2).clustered=1;
		 			sIntegrin.get(i1).zetaI = 400;// 4.114*1E-3/0.04;
		 			sIntegrin.get(i2).zetaI =  400;//4.114*1E-3/0.04;		 			
		 			sIntegrin.get(i2).Bead_Position=move(sIntegrin.get(i2).Bead_Position, sIntegrin.get(i1).Bead_Position, force, dist, sIntegrin.get(i2).dt, sIntegrin.get(i2).zetaI);
					sIntegrin.get(i1).Bead_Position=move(sIntegrin.get(i1).Bead_Position, sIntegrin.get(i2).Bead_Position, force, dist, sIntegrin.get(i1).dt, sIntegrin.get(i1).zetaI);					
		 		}		
		 		}
		 	}	
			 }
		 }
		
	}
	
	
	
	private double Distance(Integrin_3D integrin_3d1, Integrin_3D integrin_3d2) {
		double D=Math.sqrt(Math.pow(integrin_3d1.Bead_Position[0]-integrin_3d2.Bead_Position[0],2)+Math.pow(integrin_3d1.Bead_Position[1]-integrin_3d2.Bead_Position[1],2)+Math.pow(integrin_3d1.Bead_Position[2]-integrin_3d2.Bead_Position[2],2));
		return D;
	}
	
	private double[] move(double[] beadA, double[] beadB, double force, double r, double dt, double zeta) {		
		double temp[] = new double[3];		
		temp[0] = beadA[0]-force*(beadA[0]-beadB[0])/r*dt/zeta;
		temp[1] = beadA[1]-force*(beadA[1]-beadB[1])/r*dt/zeta;
		temp[2] = beadA[2]-force*(beadA[2]-beadB[2])/r*dt/zeta;	
		return temp;
	}

}
