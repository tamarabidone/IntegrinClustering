package tryout;

import java.util.ArrayList;
import java.util.Random;

public class Activation {
	public void active_deactive (ArrayList <Integrin_3D> sIntegrin, ArrayList <Ligand_3D> sLigand){
		

		 for (int i=0; i<sIntegrin.size(); i++){
			 Random rand = new Random();
			 
			 if (sIntegrin.get(i).active==0 && rand.nextDouble()<sIntegrin.get(i).ka*sIntegrin.get(i).dt)
				 sIntegrin.get(i).active=1;
			 if (sIntegrin.get(i).active==1 && rand.nextDouble()<sIntegrin.get(i).kd*sIntegrin.get(i).dt){
				 sIntegrin.get(i).active=0;
				 sIntegrin.get(i).clustered=0;
				 sIntegrin.get(i).ligandBoundYN=0;
				 sIntegrin.get(i).boundF=-1;
				 sIntegrin.get(i).boundB=-1;
				 if (sIntegrin.get(i).ligandNum>-1){
				 sLigand.get(sIntegrin.get(i).ligandNum).boundYN=0;
				 sLigand.get(sIntegrin.get(i).ligandNum).Tension=0;
				 sIntegrin.get(i).ligandNum=-1;}
				 sIntegrin.get(i).zetaI =  sIntegrin.get(i).zetaI_initial;
			 }
		 	
	
		 }
		
	}
}
