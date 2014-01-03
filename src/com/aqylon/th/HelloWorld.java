package com.aqylon.th;
import com.aqylon.aircooledCondenser.AircooledCondenser;
import com.aqylon.thermodynamics.physics.ThFluid;
import com.aqylon.thermodynamics.physics.ThState;

public class HelloWorld {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Hello from Aqylon!");
	
		double di =0.34;  //di internal diameter of tube (m)
		double Do =0.35;  //Do external diameter of tube without fins (m) 
		double Db =0.40 ;  //Db external diameter of tube with fins (m) // ask
		double b =0.001;   //b fin width (m) // ask
		double fs =0.003;  //fs fin spacing (m) 
		double L  =5;  //L condenser tubes length (m)
		double W = 2 ;  //W condenser width (m)
		int ntf = 10 ; // ntf number of tube horizontal rows // ask
		int np = 2 ; // np number of passes
		int nt = 2 ; // nt number of tubes per pass
		double airTotalMassFlow =50 ;//airTotalMassFlow air mass flow at inlet (kg/s)
		double airInletTemperature  = 290;//airInletTemperature air temperature at inlet (K)
		
		AircooledCondenser condenseur =new AircooledCondenser(di, Do, Db, b, fs, L, W,ntf,np,nt,airTotalMassFlow, airInletTemperature);
		
		ThFluid fluid=ThFluid.Water;
		ThState fluidInletState=new ThState(fluid);
		fluidInletState=fluidInletState.create();
		fluidInletState.setQuality(0.1);
		
		ThState outlet=condenseur.computeOutlet( fluidInletState, 4);
		System.out.println( outlet.quality );
		
		condenseur.printFlowPatternMap(fluid);
		
	}

}
 