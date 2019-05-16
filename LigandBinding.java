package tryout;

import java.util.ArrayList;
import java.util.Random;

public class LigandBinding {
	public void BindLigand (ArrayList<Integrin_3D> sIntegrin,ArrayList <Ligand_3D> sLigand, double cut_off){
		

		 for (int i1=0; i1<sIntegrin.size(); i1++){
			 //sIntegrin.get(i1).zetaI=4.11*1E-3/0.3;
		 	for (int l=0; l<sLigand.size(); l++){
		 		//sLigand.get(l).boundYN=0;
		 		double dist=Distance(sIntegrin.get(i1), sLigand.get(l));
		 		
		 		double x=dist-(Math.abs(sLigand.get(l).zLigand));
		 		double force=sLigand.get(l).k*x;

		 		if ((sIntegrin.get(i1).ligandBoundYN==0 && sLigand.get(l).boundYN==0 && sIntegrin.get(i1).active==1 && dist<cut_off) || (sIntegrin.get(i1).ligandBoundYN==1 && sLigand.get(l).boundYN==1 && sIntegrin.get(i1).ligandNum==l)){
		 			sIntegrin.get(i1).ligandBoundYN=1;
		 			sIntegrin.get(i1).ligandNum=l;
		 			sLigand.get(l).boundYN=1;
		 			sLigand.get(l).Bead_PositionL=move(sLigand.get(l).Bead_PositionL, sIntegrin.get(i1).Bead_Position, force, dist, sLigand.get(l).dt, sLigand.get(l).zetaL);
					sIntegrin.get(i1).Bead_Position=move(sIntegrin.get(i1).Bead_Position, sLigand.get(l).Bead_PositionL, force, dist, sIntegrin.get(i1).dt, sIntegrin.get(i1).zetaI);											
					sIntegrin.get(i1).zetaI=400;	
					
		 					
		 		}
		 	}	
		 }
		 
		 
		
		
	}
	
	
	public void UnbindLigand (ArrayList <Integrin_3D> sIntegrin,ArrayList <Ligand_3D> sLigand){

		 for (int i=0; i<sIntegrin.size(); i++){
			 Random rand = new Random();
			 if (sIntegrin.get(i).active==1 && sIntegrin.get(i).ligandBoundYN==1 && rand.nextDouble()<Math.exp(-sIntegrin.get(i).affinityL)){
				
				 sIntegrin.get(i).ligandBoundYN=0;
				
				 sLigand.get(sIntegrin.get(i).ligandNum).boundYN=0;
				 sLigand.get(sIntegrin.get(i).ligandNum).Tension=0;
				
				 sIntegrin.get(i).ligandNum=-1;
				 if (sIntegrin.get(i).clustered==0)
					 
				 sIntegrin.get(i).zetaI=sIntegrin.get(i).zetaI_initial;
			 }
			 
		 }
		 
		}
	
	private double Distance(Integrin_3D i, Ligand_3D l) {
		double D=Math.sqrt(Math.pow(i.Bead_Position[0]-l.Bead_PositionL[0],2)+Math.pow(i.Bead_Position[1]-l.Bead_PositionL[1],2)+Math.pow(i.Bead_Position[2]-l.Bead_PositionL[2],2));
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
