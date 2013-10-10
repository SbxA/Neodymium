package com.aqylon.aircooledCondenser;

import com.aqylon.thermodynamics.physics.ThFluid;
import com.aqylon.thermodynamics.physics.ThState;


/**
 * 
 * @author Vincent Lamblot
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
	private TransferUnit node;

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
	 * Working data arrays
	 */
	private double[] airTemperatures;
	private double[] deltaT;
	private ThState[] fluidNode;

	/**
	 * Computation parameters
	 */
	private int N;
	public double dA;
	private double dx = 1.0E-3; // Computation spatial step (m)

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

		airTemperatures = new double[N*(nf+1)]; // In K		
		deltaT = new double[N*nf]; // In K
		fluidNode = new ThState[(N+1)*nf];

		this.airTotalMassFlow = airTotalMassFlow; // In kg/s
		this.airNodeMassFlow = airTotalMassFlow/ntf; // In kg/s
		this.airInletTemperature = airInletTemperature; // In K

		for(int i=0; i < airTemperatures.length; i++){
			airTemperatures[i] = this.airInletTemperature;
		}

	}

	public ThState computeOutlet(ThState fluidInletState, double fluidTotalMassFlow){
		//TODO
		this.fluidTotalMassFlow = fluidTotalMassFlow;
		this.fluidNodeMassFlow = fluidTotalMassFlow/nf;
		ThState finalMixture = new ThState(ThFluid.Water);

		for(int i=0; i<100; i++){
			// use unit.airOutletTemperature
		}
		return null;
	}



	/*
	 * Utilities
	 */

	/**
	 * 
	 * @param i index of integration node along the tube
	 * @param j index of tube
	 * @return air temperature at node (i,j)
	 */
	private double getAirTemperature(int i, int j){
		return airTemperatures[i+nf*j];
	}


	/**
	 * 
	 * @param i index of integration node along the tube
	 * @param j index of tube
	 * @param temperature air temperature to set un node (i,j)
	 */
	private void setAirTemperature(double temperature, int i, int j){
		airTemperatures[i+nf*j] = temperature;
	}


	/**
	 * 
	 * @param i index of integration node along the tube
	 * @param j index of tube
	 * @return air temperature difference at node (i,j)
	 */
	private double getDeltaT(int i, int j){
		return deltaT[i+nf*j];
	}


	/**
	 * 
	 * @param i index of integration node along the tube
	 * @param j index of tube
	 * @param temperature air temperature difference to set un node (i,j)
	 */
	private void setDeltaT(double dt, int i, int j){
		deltaT[i+nf*j] = dt;
	}

	/**
	 * 
	 */
	private void updateAirTemperatureMap(){
		for(int j=nf; j>=0; j--){
			for(int i=0; i<N; i++){
				setAirTemperature(getAirTemperature(i, j)+getDeltaT(i, j), i, j);
			}
		}
	}

}
