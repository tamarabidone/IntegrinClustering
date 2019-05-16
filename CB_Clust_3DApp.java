package tryout;
import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.JOptionPane;








public class CB_Clust_3DApp {

	/*An array of filament instances*/
	ArrayList<Integrin_3D> sIntegrin=new ArrayList<Integrin_3D>();
	ArrayList<Ligand_3D> sLigand=new ArrayList<Ligand_3D>();
	ArrayList<Filament_3D> sFilament=new ArrayList<Filament_3D>();
	

	
	double displayInterval;
	double eta;	//Viscosity
	double steptime;

	int nOfIntegrins, nOfLigands, nOfFil;

	
    Random rand = new Random();
   
	tangentPlot tP;	//Plots tangent correlation on the contour of a filament

	/*Plotting and displaying frames*/

	
	
	
	public void initialize(){


		eta = 0.3;
		steptime = 0.1;
		nOfIntegrins = 100;
		nOfLigands = 1000;
		nOfFil = 0;
	 double DomainSize=1;
	 double dt=1E-4;
				
		for(int i=0;i<nOfFil;i++){
			sFilament.add(new Filament_3D());
		
			sFilament.get(i).DomainSize = DomainSize;
			sFilament.get(i).dt=dt;
			sFilament.get(i).V=700*1E-3;
			/*Initial positions for the first beads*/
			double[][] initialpos = new double[2][3];		
			initialpos[0][0] = -sFilament.get(i).DomainSize/2+Math.random()*sFilament.get(i).DomainSize;
			initialpos[0][1] = -sFilament.get(i).DomainSize/2+Math.random()*sFilament.get(i).DomainSize;
			initialpos[0][2] = 0.02;
			sFilament.get(i).start(initialpos, Math.random()*360);
			
		}
		for(int i=0;i<nOfIntegrins;i++){		
			sIntegrin.add(new Integrin_3D());

			sIntegrin.get(i).DomainSize = DomainSize;
			sIntegrin.get(i).dt=dt;//time step	
			sIntegrin.get(i).clustered=0;
			sIntegrin.get(i).temperature = 300;
			sIntegrin.get(i).zIntegrin = 0;			
			sIntegrin.get(i).zetaI_initial =  4.114*1E-3/0.3;//0.3;//Math.PI*(eta*(1+Math.random()*1000))*sIntegrin.get(i).Length_R/(0.84 + Math.log(sIntegrin.get(i).Length_R/(2*0.010)));
			sIntegrin.get(i).zetaI =sIntegrin.get(i).zetaI_initial;
			sIntegrin.get(i).ka=0.5;//0.5;
		    sIntegrin.get(i).kd=1;//0.1;	   	    	
			sIntegrin.get(i).cut_off=0.03;
			sIntegrin.get(i).ligandBoundYN=0;
			sIntegrin.get(i).ligandNum=-1;
			sIntegrin.get(i).boundF=-1;
			sIntegrin.get(i).boundB=-1;
			double[]initialpos = new double[3];			
			initialpos[0]=-sIntegrin.get(i).DomainSize/2+Math.random()*sIntegrin.get(i).DomainSize;
			initialpos[1]=-sIntegrin.get(i).DomainSize/2+Math.random()*sIntegrin.get(i).DomainSize;
			initialpos[2] = sIntegrin.get(i).zIntegrin;
			sIntegrin.get(i).start(initialpos);		
		/*	if (i<25){
		    sIntegrin.get(i).affinityI=9; 		    	sIntegrin.get(i).affinityL=9;
			}
			else{*/		 	
			    sIntegrin.get(i).affinityI=1;
			    sIntegrin.get(i).affinityL=111;		  				
			//}
			sIntegrin.get(i).v=2;		
		}
		
		for(int i=0;i<nOfLigands;i++){		
			sLigand.add(new Ligand_3D()); 	sLigand.get(i).DomainSize = DomainSize;
			sLigand.get(i).Length_L = 0.01;	sLigand.get(i).dt=dt;//time step	
			sLigand.get(i).temperature = 300;	sLigand.get(i).zLigand = -0.02; sLigand.get(i).k=1;
			sLigand.get(i).Tension=0; 	sLigand.get(i).numBoundIntegrins=0;
			sLigand.get(i).zetaL =1000;//4*Math.PI*eta/100*sLigand.get(i).Length_L/(0.84 + Math.log(sLigand.get(i).Length_L/2/0.0035)); //1000;//
			sLigand.get(i).boundYN =0; 			double[]initialpos = new double[3];	
			initialpos[0]=-sLigand.get(i).DomainSize/2+Math.random()*sLigand.get(i).DomainSize;
			initialpos[1]=-sLigand.get(i).DomainSize/2+Math.random()*sLigand.get(i).DomainSize;
			initialpos[2] = sLigand.get(i).zLigand;			sLigand.get(i).start(initialpos);
		    		}
		}
	
	
	int I1=1, L1=111, F=103;    

	public void doStep() {

		/*performs a steptime run of the simulation*/
		for (int dv=0; dv < (steptime<sIntegrin.get(0).dt?1:steptime)/sIntegrin.get(0).dt;dv++) 
			doStep2();
		}

	
	
	
    /*The "REAL" step*/
	public void doStep2() {		
		int num=1;
		int numero=0;
			Activation activate=new Activation ();
			activate.active_deactive(sIntegrin, sLigand);			
			
			LigandBinding lig=new LigandBinding();
		 	lig.BindLigand(sIntegrin,sLigand, Math.abs(sLigand.get(0).zLigand));			
		
			PairwiseInteractions pair=new PairwiseInteractions();
		   pair.Gaussian (sIntegrin, sIntegrin.get(0).cut_off);	   
		
		
		
		Decluster dec=new Decluster ();
	dec.declustering(sIntegrin);		
		
		lig.UnbindLigand(sIntegrin, sLigand);		
	
		CatchBond catchb=new CatchBond();

		int numCatchUnbound=catchb.ApplyCatchB(numero, sIntegrin, sLigand);	
		
		Data d=new Data();		
		int LigandNum=1000;
		
		//ligands1000
		String BoundLClusteredIActiveI= "ImplicitFlow_B_G_"+ Double.toString(I1)+"_"+Double.toString(L1)+"_"+Double.toString(F)+"_"+Integer.toString(num);
		String IntegrinCoordinates= "ImplicitFlow_B_I_"+ Double.toString(I1)+"_"+Double.toString(L1)+"_"+Double.toString(F)+"_"+Integer.toString(num);
		String LigandCoordinates= "ImplicitFlow_B_L_"+ Double.toString(I1)+"_"+Double.toString(L1)+"_"+Double.toString(F)+"_"+Integer.toString(num);
		if (sIntegrin.get(1).time%1<0.0001) {
			{
				try
				{
						
				// System.out.println(j+"    "+sIntegrin.get(j).Bead_Position[0]);
				    FileWriter fw = new FileWriter(BoundLClusteredIActiveI,true); //the true will append the new data
				    fw.write(String.format("%.3f", sIntegrin.get(0).time)+"    "+"    "+d.Lnum(sLigand)+"     "+ d.ActiveNum(sIntegrin)+"       "+ d.Cnum(sIntegrin)+"       "+ d.Lbound(sIntegrin)+"        "+d.BoundActin(sIntegrin) +"    "+String.format("%.4f",d.TensionSign(sLigand))+"    "+String.format("%.4f",d.Tension(sLigand))+"   "+numCatchUnbound +" \n");//appends the string to the file
					
				    fw.close();
				    
				    FileWriter fw2 = new FileWriter(IntegrinCoordinates,true); //the true will append the new data
				    for (int i=0; i<nOfIntegrins; i++){  
				    	fw2.write(String.format("%.3f", sIntegrin.get(0).time)+"       "+i+"       "+sIntegrin.get(i).Bead_Position[0]+"        "+ sIntegrin.get(i).Bead_Position[1]+"        "+sIntegrin.get(i).active +"        "+sIntegrin.get(i).clustered +"       "+sIntegrin.get(i).ligandBoundYN+"         "+sIntegrin.get(i).ligandNum +" \n");//appends the string to the file
					   
				    }
				    	fw2.close();
				    	
				    	FileWriter fw3 = new FileWriter(LigandCoordinates,true); //the true will append the new data
					    for (int i=0; i<nOfLigands; i++){  
					    	fw3.write(String.format("%.3f", sIntegrin.get(0).time)+"       "+i+"       "+sLigand.get(i).Bead_PositionL[0]+"        "+ sLigand.get(i).Bead_PositionL[1]+"        "+sLigand.get(i).boundYN +"        "+sLigand.get(i).Tension +" \n");//appends the string to the file
						   
					    }
					    	fw3.close();
				    
				
				}
				catch(IOException ioe)
				{
				    System.err.println("IOException: " + ioe.getMessage());
				}
				
				}
		}
		for(int j=0;j<nOfIntegrins;j++){		
			sIntegrin.get(j).step();
		}


		for(int j=0;j<nOfLigands;j++){		
			sLigand.get(j).step();
			
		}
		

		
		for(int j=0;j<nOfFil;j++){		
			sFilament.get(j).step();
		}
	
		
		

		
		
			
	}
	

	


	public static String getSaveFileName(Frame parent){
		   FileDialog fd = new FileDialog(parent,"Save File",FileDialog.SAVE);
		   fd.setVisible(true);
		   String fname = fd.getFile();
		   String dirname = fd.getDirectory();
		   String fullname = dirname +  fname;
		   return fullname;
		}
/*	public void writeData(String value){
		String line = value;
	 
		try {
		RMSWriter.write(line);			 
		RMSWriter.newLine();
		RMSWriter.flush();
		}
		catch (Exception e){e.printStackTrace();
		
	 	}
	}*/

	


    public static void main(String[] args) {
        CB_Clust_3DApp simulation = new CB_Clust_3DApp  ();
         simulation.initialize();
          
         for (int i=0;i<(int)(500/simulation.steptime);i++){
             simulation.doStep();
         }
          
         System.exit(0);
     }
	
    protected static double magnitude(double[] a){
        double mag = 0;
        for(int i = 0; i<3; i++)
            mag += Math.pow(a[i],2);
        return Math.sqrt(mag);
    }
	
}

class tangentPlot{
	double[] taverage;
	double[] tstd;
	int numofds;
}
