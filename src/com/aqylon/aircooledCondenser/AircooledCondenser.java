package com.aqylon.aircooledCondenser;

import com.aqylon.thermodynamics.physics.ThFluid;
import com.aqylon.thermodynamics.physics.ThState;


/**
 * 
 * @author Vincent Lamblot & Tristan Agaësse
 *
 */
public class AircooledCondenser {

	/**
	 * General parameters
	 */
	public final static double g = 9.81; // Gravity of Earth (m/s^2)

	/**
	 * Tube geometric parameters
	 */
	public double di,Do,Db,Ao,Ad,Af,b;
	public int Nm; 

	/**
	 * Condenser geometric parameters
	 */
	public double L,W;
	public int nf, ntf, nt, np;

	
	/**
	 * Air flow parameters
	 */
	public double airTotalMassFlow, airNodeMassFlow, airInletTemperature;
	public final static double cpAir = 1.004E3; // Air thermal capacity at 20°C in J/(kg.K)
	public final static double PrAir = 0.713; // Prandtl number for air flow
	public final static double muAir = 1.983E-5; // Air dynamic viscosity at 20°C in kg/(m.s)
	public final static double kAir = 0.0257; // Air thermal conductivity at 20°C in W/(m.K)

	/**
	 * Fins thermal properties
	 */
	public final static double kFins = 205; // Aluminium thermal conductivity at 20°C in W/(m.K)

	/**
	 * Fluid flow parameters
	 */
	private double fluidTotalMassFlow, fluidNodeMassFlow;
	private static double sigma = 1.0; //FIXME



	/**
	 * Computation parameters
	 */
	public int N;
	public double dA;
	public double dx = 1.0E-3; // Computation spatial step (m)

	/**
	 * Geometric parameters as shown in Fig.8.8
	 * E. CAO - "Heat transfer in process engineering" - Ch.8 "Finned tubes" - Sec. 8-2-2 "Heat Transfer in Air Coolers" - pp.236-241
	 * @param di internal diameter of tube (m)
	 * @param Do external diameter of tube without fins (m) 
	 * @param Db external diameter of tube with fins (m) // ask
	 * @param b fin width (m) // ask
	 * @param fs fin spacing (m) 
	 * @param L condenser tubes length (m)
	 * @param W condenser tubes width (m)
	 * @param ntf number of tube horizontal rows // ask
	 * @param airTotalMassFlow air mass flow at inlet (kg/s)
	 * @param airInletTemperature air temperature at inlet (K)
	 */
	public AircooledCondenser(double di, double Do, double Db, double b, double fs, double L, int ntf, double airTotalMassFlow, double airInletTemperature){

		this.di = di;
		this.Do = Do;
		this.Db = Db;
		this.b = b;
		this.Nm = (int)(1/fs);
		this.L = L;
		this.ntf = ntf;

		Ao = Math.PI*Do; // Plain tube area per unit tube length (m^2)
		Ad = Ao*(1-b*Nm); // Exposerd area of tube (not covered by fins) per unit tube length (m^2)
		Af = Math.PI*2*Nm*(Db*Db-Do*Do)/4; // Fin area per unit tube length (m^2)

		np = 2; // Number of passes
		nt = 2; // Number of tube(s) per pass

		N = (int)(L/dx); // First size of the computation arrays
		dx = L/N;
		dA = Ao*dx; // Reference exchange area for each node (m^2)
		nf = np*nt; // Second size of the computation arrays


		this.airTotalMassFlow = airTotalMassFlow; // In kg/s
		this.airNodeMassFlow = airTotalMassFlow/ntf; // In kg/s   ????????????????
		this.airInletTemperature = airInletTemperature; // In K
		
	}

	public ThState computeOutlet(ThState fluidInletState, double fluidTotalMassFlow){
		
		this.fluidTotalMassFlow = fluidTotalMassFlow;
		this.fluidNodeMassFlow = fluidTotalMassFlow/(nt*ntf);

		NumericalSolver solver=new NumericalSolver(this, fluidInletState, fluidTotalMassFlow);
    
		ThState[] statesToMix=solver.computeOutletStateOfEachPipe();
		
    ThState outletState = mixSeveralStates(statesToMix);
		return outletState;
	}
	
	
	
  public ThState mixSeveralStates(ThState[] statesToMix){
    
    ThFluid fluid=statesToMix[0].fluid;
    ThState mixedState =new ThState(fluid);
    
    double quality = 0;
    for(int i=0; i<statesToMix.length ;i++ ){
      quality=statesToMix[i].quality;
    }
    quality=quality/statesToMix.length ;
        
    mixedState.setQuality( quality );
    
    return mixedState ;
  }
  
  
  
  
  
  public void printFlowPatternMap(ThFluid fluid){
    
    int nPoints=50;
    double[] mWavyArray=new double[nPoints];
    double[] mStratArray=new double[nPoints];
    
    ThState fluidState = new ThState(fluid);
    double airNodeInletTemperature=300;
    
    for(int i=0;i<nPoints+1;i++){
       double x=i/nPoints;
       fluidState.setQuality(x);
       HeatTransferLocalUnit transferUnit = new HeatTransferLocalUnit(this, airNodeInletTemperature, fluidState, this.fluidNodeMassFlow);
       transferUnit.computeFlowPatternMapBoundaries() ;
    
       mWavyArray[i]=transferUnit.mWavy ;
       mStratArray[i]=transferUnit.mStrat ;
    }
   
    
  }
  
  
}
