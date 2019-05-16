package tryout;

import java.util.ArrayList;
import java.util.Random;


public class Data {

	public int BoundActin (ArrayList <Integrin_3D> sIntegrin){
int num=0;
		 for (int i=0; i<sIntegrin.size(); i++){
			
			 if (sIntegrin.get(i).boundF>-1 ){
				 num=num+1;
			 }
			 
		 }
		 return num;
		 
		}
	public double Tension (ArrayList <Ligand_3D> sLigand){
double num=0;
double T=0;
double avT=0;
		 for (int i=0; i<sLigand.size(); i++){
			
			 if (sLigand.get(i).boundYN>-1){
				 T=T+Math.abs(sLigand.get(i).Tension);
				 num=num+1;
				// System.out.println(T);
			 }
			 
		 }
		 
		 avT=T/num;
		 return T;
		 
		}
	
	public double TensionSign (ArrayList <Ligand_3D> sLigand){
double num=0;
double T=0;
double avT=0;
		 for (int i=0; i<sLigand.size(); i++){
			
			 if (sLigand.get(i).Tension!=0){
				 T=T+(sLigand.get(i).Tension);
				 num=num+1;
				// System.out.println(T);
			 }
			 
		 }
		 
		 avT=T/num;
		 return T;
		 
		}
	
	public int Lnum (ArrayList <Ligand_3D> sLigand){
int num=0;
		 for (int i=0; i<sLigand.size(); i++){
			
			 if (sLigand.get(i).boundYN==1 ){
				 num=num+1;
			 }
			 
		 }
		 return num;
		 
		}
	
	

	public int Cnum (ArrayList <Integrin_3D> sIntegrin){
int num=0;
		 for (int i=0; i<sIntegrin.size(); i++){
			
			 if (sIntegrin.get(i).clustered==1 ){
				 num=num+1;
			 }
			 
		 }
		 return num;
		 
		}
	
	public int ActiveNum (ArrayList <Integrin_3D> sIntegrin){
int num=0;
		 for (int i=0; i<sIntegrin.size(); i++){
			
			 if (sIntegrin.get(i).active==1 ){
				 num=num+1;
			 }
			 
		 }
		 return num;
		 
		}
	
	public int Lbound (ArrayList <Integrin_3D> sIntegrin){
int num=0;
		 for (int i=0; i<sIntegrin.size(); i++){
			
			 if (sIntegrin.get(i).ligandNum>-1 ){
				 num=num+1;
			 }
			 
		 }
		 return num;
		 
		}
}
