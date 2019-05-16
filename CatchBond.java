package tryout;

import java.util.ArrayList;
import java.util.Random;

public class CatchBond {

	public int ApplyCatchB (int numero,ArrayList <Integrin_3D> sIntegrin, ArrayList <Ligand_3D> sLigand){
		int num=0;
		 for (int i1=0; i1<sIntegrin.size()/1-1; i1++){
				if (sIntegrin.get(i1).ligandBoundYN==1 )
				{
					 Random rand = new Random();
					 sIntegrin.get(i1).Bead_Position[0]=sIntegrin.get(i1).Bead_Position[0]+ sIntegrin.get(i1).v*sIntegrin.get(i1).dt;
		 		double dist=Distance(sIntegrin.get(i1), sLigand.get(sIntegrin.get(i1).ligandNum));
		 		
		 		double eq=Math.abs(sLigand.get(sIntegrin.get(i1).ligandNum).zLigand);
		 		
		 		double Tension =50*(sLigand.get(sIntegrin.get(i1).ligandNum).k*(dist-eq));
		 		sLigand.get(sIntegrin.get(i1).ligandNum).Tension=Tension;
		 	
		 		//if (i1==4)
			 //	System.out.println(i1 + "   dist   "+dist +"    Ligand     "+ (sIntegrin.get(i1).ligandNum)+" tension    "+ Tension);
		 		//double k01=3, lambda1=0.25, lambda2=0.25, k02=0.042;
		 		//double k01=0.8,  k02=0.000078;
		 		double lambda1=0.8, lambda2=0.088;
		 		double k01=0.22,  k02=0.000078;
		 	//kong double kub=0.9*Math.exp(-0.005*Tension+3E-7*Math.exp(0.4*Tension));
		 		
		 		////DDDdouble kub=0.4*Math.exp(-0.004*Tension+4E-7*Math.exp(0.2*Tension));
		 		
		 		double kub=1.6*Math.exp(-0.1*Tension+5E-2*Math.exp(0.05*Tension));
/*	if (Math.abs(Tension)<30){
		 			kub=k01*Math.exp(lambda1*Math.abs((Tension)/(4.11*1E-3)));
		 		}
		 		else {
		 			kub=k02*Math.exp(-lambda2*Math.abs(Tension)/(4.11*1E-3));
		 		}*/
		 		//System.out.println(Tension);
	
		 		if (rand.nextDouble() <kub*sIntegrin.get(i1).dt)
		 		{
		 			System.out.println(Tension+"    "+kub);
		 			//System.out.println("unbound Ligand");
		 			 num=num+1;
		 			// System.out.println(i1+"    "+numero);
		 			sIntegrin.get(i1).ligandBoundYN=0;
		 			sIntegrin.get(i1).boundF=-1;
		 			sIntegrin.get(i1).boundB=-1;
		 			 sLigand.get(sIntegrin.get(i1).ligandNum).boundYN=0;
		 			sLigand.get(sIntegrin.get(i1).ligandNum).Tension=0;
		 			sIntegrin.get(i1).ligandNum=-1;
		 			 if (sIntegrin.get(i1).clustered==0)
						 
						 sIntegrin.get(i1).zetaI=sIntegrin.get(i1).zetaI_initial;
		 			
		 					
		 		}
				}
				else
				{
		 			//System.out.println(Tension+"    "+kub);
		 			//System.out.println("unbound Ligand");
		 			
		 			// System.out.println(i1+"    "+numero);
		 			sIntegrin.get(i1).ligandBoundYN=0;
		 			sIntegrin.get(i1).boundF=-1;
		 			sIntegrin.get(i1).boundB=-1;
		 			// sLigand.get(sIntegrin.get(i1).ligandNum).boundYN=0;
		 			//sLigand.get(sIntegrin.get(i1).ligandNum).Tension=0;
		 			sIntegrin.get(i1).ligandNum=-1;
		 			 if (sIntegrin.get(i1).clustered==0)
						 
						 sIntegrin.get(i1).zetaI=sIntegrin.get(i1).zetaI_initial;
		 			
		 					
		 		}
		 }
		 		
		 return num;	
		 	
		
	}
	
	private double Distance(Integrin_3D integrin_3d1, Ligand_3D lig) {
		double D=Math.sqrt(Math.pow(integrin_3d1.Bead_Position[0]-lig.Bead_PositionL[0],2)+Math.pow(integrin_3d1.Bead_Position[1]-lig.Bead_PositionL[1],2)+Math.pow(integrin_3d1.Bead_Position[2]-lig.Bead_PositionL[2],2));
		return D;
	}
	
}
