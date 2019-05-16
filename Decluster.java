package tryout;

import java.util.ArrayList;
import java.util.Random;

public class Decluster {
	
	public void declustering (ArrayList <Integrin_3D> sIntegrin){

	 for (int i=0; i<sIntegrin.size(); i++){
		 Random rand = new Random();
		 if (sIntegrin.get(i).clustered==1 && rand.nextDouble()<Math.exp(-sIntegrin.get(i).affinityI)){
			 sIntegrin.get(i).clustered=0;
			 
			 if (sIntegrin.get(i).ligandBoundYN==0){
			 sIntegrin.get(i).zetaI = sIntegrin.get(i).zetaI_initial;
			 //sLigand.get(sIntegrin.get(i).ligandNum).Tension=0;
			 }
		 }
		 
	 }
	 
	}
}
