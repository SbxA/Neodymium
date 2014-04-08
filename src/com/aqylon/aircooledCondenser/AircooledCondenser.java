package com.aqylon.aircooledCondenser;

import java.text.DecimalFormat;
import java.text.Format;

import com.aqylon.thermodynamics.physics.ThFluid;
import com.aqylon.thermodynamics.physics.ThState;

/**
 * 
 * @author Vincent Lamblot 
 *
 */
public class AircooledCondenser {

	/**
	 * Standard format 0.##
	 */
	public final static Format std = new DecimalFormat("0.##");

	/**
	 * Gravity of Earth (m/s^2)
	 */
	public final static double g = 9.81;
	
	

	/*
	 * Parameters as shown in Fig.8.8
	 * E. CAO - "Heat transfer in process engineering" - Ch.8 "Finned tubes" - Sec. 8-2-2 "Heat Transfer in Air Coolers" - pp.236-241
	 */

	/**
	 * Internal diameter (in m)
	 */
	public double di;

	/**
	 * External diameter of tube without fins (m) 
	 */
	public double Do;

	/**
	 * External diameter of tube with fins (m)
	 */
	public double Db ;

	/**
	 * Plain tube area per unit tube length (m^2)
	 */
	public double Ao;

	/**
	 * Exposed area of tube (not covered by fins) per unit tube length (m^2)
	 */
	public double Ad;

	/**
	 * Fin area per unit tube length (m^2)
	 */
	public double Af;

	/**
	 * Condenser tubes length (m)
	 */
	public double L;

	/**
	 * Condenser width (m)
	 */
	public double W;

	/**
	 * Number of tube horizontal rows
	 */
	public int ntf;

	/**
	 * Number of passes
	 */
	public int np ;

	/**
	 * Number of tube(s) per pass
	 */
	public int nt ;

	/**
	 * Number of lines in the computation arrays
	 */
	public int nf;

	// Air flow parameters

	/**
	 * Air mass flow at inlet (kg/s)
	 */
	public double airTotalMassFlow  ;

	/**
	 * Air temperature at inlet (K)
	 */
	public double airInletTemperature  ;

	/**
	 * Air mass flow in each node (kg/s)
	 */
	public double  airNodeMassFlow ;

	/**
	 * Air thermal capacity in J/(kg.K)
	 */
	public final static double cpAir = 1.004E3; // Air thermal capacity at 20°C 

	/**
	 * Prandtl number for air flow
	 */
	public final static double PrAir = 0.713; 

	/**
	 * Air dynamic viscosity in kg/(m.s)
	 */
	public final static double muAir = 1.983E-5; // Air dynamic viscosity at 20°C 

	/**
	 * Air thermal conductivity in W/(m.K)
	 */
	public final static double kAir = 0.0257; // Air thermal conductivity at 20°C 

	/**
	 * fin width (m)
	 */
	public double b;

	/**
	 * fin spacing (m) 
	 */
	public double fs;

	/**
	 * Number of fins per unit length of tube (m^-1)
	 */
	public int Nm; 

	/**
	 * Fins thermal conductivity in W/(m.K)
	 */
	public final static double kFins = 205; // Aluminium thermal conductivity at 20°C 

	/**
	 * Total fluid mass flow arriving in the condenser
	 */
	private double fluidTotalMassFlow; 

	/**
	 * Fluid mass flow in one tube
	 */
	double fluidNodeMassFlow;

	/**
	 * Number of computational subdivisions steps along the pipes
	 */
	public int N=10; 

	/**
	 * Area of a pipe subdivision step : exchange area for each node (m^2)
	 */
	public double dA;

	/**
	 * Length of a pipe subdivion step
	 */
	public double dx ; 

	/**
	 * @param di internal diameter of tube (m)
	 * @param Do external diameter of tube without fins (m) 
	 * @param Db external diameter of tube with fins (m)
	 * @param b fin width (m)
	 * @param fs fin spacing
	 * @param L condenser tubes length
	 * @param W condenser width
	 * @param ntf number of tube horizontal rows 
	 * @param np number of passes
	 * @param nt number of tubes per pass
	 * @param airTotalMassFlow air mass flow at inlet (kg/s)
	 * @param airInletTemperature air temperature at inlet (K)
	 */
	public AircooledCondenser(double di, double Do, double Db, double b, double fs, double L, double W, int ntf, int np, int nt, double airTotalMassFlow, double airInletTemperature){

		this.di = di;
		this.Do = Do;
		this.Db = Db;
		this.b = b;
		this.Nm = (int)(1/fs);
		this.L = L;
		this.W=W;
		this.ntf = ntf;
		this.airTotalMassFlow = airTotalMassFlow;  
		this.airInletTemperature = airInletTemperature;

		Ao = Math.PI*this.Do; 
		Ad = Ao*(1-this.b*Nm); 
		Af = Math.PI*2*Nm*(this.Db*this.Db-this.Do*this.Do)/4; 

		this.np = np; 
		this.nt = nt;

		dx = L/N;
		dA = Ao*dx; 
		nf = np*nt; 
		airNodeMassFlow = this.airTotalMassFlow*(dx/this.L)*(this.W/nf);
	}

	public ThState computeOutlet(ThState fluidInletState, double fluidTotalMassFlow) throws Exception{

		this.fluidTotalMassFlow = fluidTotalMassFlow;
		fluidNodeMassFlow = this.fluidTotalMassFlow/(nt*ntf);

		NumericalSolver solver = new NumericalSolver(this, fluidInletState, fluidTotalMassFlow);

		ThState[] statesToMix = solver.computeOutletStateOfEachPipe();

		ThState outletState = mixSeveralStates(statesToMix);

		return outletState;
	}

	private ThState mixSeveralStates(ThState[] statesToMix){

		ThFluid fluid=statesToMix[0].fluid;
		ThState mixedState =new ThState(fluid);
		mixedState=mixedState.create();

		double quality = 0;
		for(int i=0; i<statesToMix.length ;i++ ){
			quality=statesToMix[i].quality;
		}
		quality=quality/statesToMix.length ;

		mixedState.setQuality( quality );

		return mixedState ;
	}

	public void printFlowPatternMap(ThFluid fluid){
		//compute the flow pattern map boundaries for several qualities. Useful to debug the flow pattern map.
		int nPoints=50;
		double[] mWavyArray=new double[nPoints];
		double[] mStratArray=new double[nPoints];
		double quality;

		double airNodeInletTemperature=300;

		for(int i=1;i<nPoints-1;i++){
			quality=(double) i/nPoints;
			ThState fluidState = new ThState(fluid);
			fluidState=fluidState.create();
			fluidState.setQuality(quality);

			LocalHeatTransferUnit transferUnit = new LocalHeatTransferUnit(airNodeInletTemperature, fluidState);
			transferUnit.computeFlowPatternMapBoundaries() ;

			mWavyArray[i]=transferUnit.mWavy ;
			mStratArray[i]=transferUnit.mStrat ;

		}
	}

	/**
	 * Compute the local biphasic heat exchanges based on the model developped in
	 * E. CAO - "Heat transfer in process engineering" - Ch.8 "Finned tubes" - Sec. 8-2-2 "Heat Transfer in Air Coolers" - pp.236-241
	 */
	public class LocalHeatTransferUnit {

		/**
		 * Temperature of the air arriving in this node
		 */
		private double airNodeInletTemperature ;
		/**
		 * Thermodynamical state of the fluid arriving in this node
		 */
		private ThState fluidNodeInletState;
		/**
		 * Computed temperature gain of air at this node due to heat exchange
		 */
		private double airNodeOutletDeltaT;


		public double totalTransferCoefficient, enthalpyTransfered, fluidTransferCoefficient, airTransferCoefficient;

		/**
		 * Mass velocity transition between annular and stratified-wavy flows patterns
		 */
		public double mWavy;
		/**
		 * Mass velocity transition between stratified-wavy and fully stratified flows patterns
		 */
		public double mStrat;
		/**
		 * Adimensionned height of liquid in the pipe (case of a stratified flow)
		 */
		public double hLd;
		/**
		 * Local void fraction
		 */
		public double epsilon;
		/**
		 * Stratified angle around upper perimeter of the tube to stratified liquid level
		 */
		public double thetaStrat;
		/**
		 * Area of a tube cross section (m²)
		 */
		public double A;
		/**
		 * Cross-sectional area occupied by the liquid phase
		 */
		public double AL;
		/**
		 * Falling film angle around the top perimeter of the pipe
		 */
		public double theta;

		/**
		 * Fluid quality
		 */
		private double x;

		/**
		 * Fluid saturated liquid state density
		 */
		private double rhoL;

		/**
		 * Fluid saturated vapor state density
		 */
		private double rhoG;

		/**
		 * Fluid saturated liquid state dynamic viscosity
		 */
		private double muL;
		/**
		 * Fluid saturated liquid state isobaric heat capacity
		 */
		private double cpL;
		/**
		 * Fluid liquid-vapor surface tension
		 */
		private double sigma;
		/**
		 * Fluid saturated liquid state thermal conductivity
		 */
		private double kL;

		/**
		 * Local perimeter-averaged condensing heat transfer coefficient
		 */
		double alpha;
		/**
		 *  Convective heat transfer coefficient 
		 */
		double alphaC=0;
		/**
		 * Film condensation heat transfer
		 */
		double alphaF=0;
		/**
		 * Interfacial roughness correction factor for turbulent liquid film condensation heat transfer 
		 */

		double fi;
		/**
		 * Liquid film Reynolds number 
		 */
		double ReL;
		/**
		 * Liquid film Prandtl number 
		 */
		double PrL;
		/**
		 * Liquid film thickness 
		 */
		double delta;

		/**
		 *  Local surfacic heat transfer between interior and exterior of a pipe
		 */
		double q;
		/**
		 * Ratio of gas and liquid velocities
		 */
		double velocitiesRatio;


		/**
		 * @param condenser
		 * @param airNodeInletTemperature
		 * @param fluidNodeInletState
		 * @param fluidNodeMassFlow
		 */
		public LocalHeatTransferUnit(double airNodeInletTemperature, ThState fluidNodeInletState){

			this.fluidNodeInletState = fluidNodeInletState;
			this.airNodeInletTemperature = airNodeInletTemperature;

			x = fluidNodeInletState.quality;
			rhoL = fluidNodeInletState.saturatedLiquidState.density;
			rhoG = fluidNodeInletState.saturatedVaporState.density;
			muL = fluidNodeInletState.saturatedLiquidState.dynamicViscosity;
			cpL = fluidNodeInletState.saturatedLiquidState.isobaricHeatCapacity;
			kL = fluidNodeInletState.saturatedLiquidState.thermalConductivity;
			sigma = fluidNodeInletState.saturatedLiquidState.surfaceTension;

			A=Math.PI*di*di/4;
		}

		public void computeTransfer() throws Exception{

			airTransferCoefficient = computeAirTransferCoefficient();

			// Start automatic computation
			this.enthalpyTransfered = 1000.0; // In J, 1st hypothesis
			q=enthalpyTransfered/dA;
			double enthalpyTransferedBuffer = 0.0;

			while(Math.abs(enthalpyTransfered-enthalpyTransferedBuffer)/enthalpyTransfered > 1.0E-5){
				enthalpyTransferedBuffer = enthalpyTransfered;
				fluidTransferCoefficient = computeFluidTransferCoefficient();
				totalTransferCoefficient = 1/(1/fluidTransferCoefficient+1/airTransferCoefficient);
				enthalpyTransfered=totalTransferCoefficient*dA*(fluidNodeInletState.temperature-airNodeInletTemperature); // Chosen as positive in the fluid to air way.
				q=enthalpyTransfered/dA;
			}

			System.out.println("Air heat transfer coefficient : "+std.format(airTransferCoefficient)+" W.m-2.K-1");
			System.out.println("Fluid heat transfer coefficient : "+std.format(fluidTransferCoefficient)+" W.m-2.K-1");
			System.out.println("Total heat transfer coefficient : "+std.format(totalTransferCoefficient)+" W.m-2.K-1");

			airNodeOutletDeltaT = enthalpyTransfered/(cpAir*airNodeMassFlow);
		}


		/**
		 * alpha as defined in eq. 8.1.23.
		 * E. CAO - "Heat transfer in process engineering" - Ch.8 "Finned tubes" - Sec. 8-2-2 "Heat Transfer in Air Coolers" - pp.236-241
		 * Computing based on the Thome-El Hajal-Calvini flow pattern method
		 * Wolverine Tube, INC. - Ch.8.1 "Condensation inside Horizontal Tubes" - "Heat transfer model" pp.8.8-8.12
		 * @return alpha
		 * @throws Exception 
		 */
		private double computeFluidTransferCoefficient() throws Exception{

			this.epsilon=computeEpsilon(x);
			this.thetaStrat=computeThetaStrat(epsilon);
			this.AL = (1-epsilon)*A;
			
			FlowPattern pattern = getPattern();

			switch(pattern){
			case AnnularFlow:
				theta = 0.0;
				break;
			case StratifiedFlow:
				theta = thetaStrat;
				break;
			case StratifiedWavyFlow:
				theta = thetaStrat*Math.pow((mWavy-fluidNodeMassFlow)/(mWavy-mStrat), 0.5); // eq. (8.1.30)
				break;
			default:
				throw new Exception("Flow pattern not defined. Do no know how to solve this case !!");
			}


			if(x==0){theta=0.0;}  //ill defined when x=0

			velocitiesRatio = rhoL*x*(1-epsilon)/(rhoG*epsilon*(1-x)); // eq. (8.1.36/37)
			if(x==0){velocitiesRatio=0;} //ill defined when x=0


			if(epsilon>0.5){
				delta=(di-Math.sqrt(di*di-8*AL/(2*Math.PI-theta)))/2; //modification of equation 8.1.35 to make it homogeneous. TODO : check the new formula //FIXME ask LC

				//quadraticEquation eq8135 = new quadraticEquation(4.0, -4*di*di, 8*AL/(2*Math.PI-theta)+Math.pow(di, 4)-di*di); 
				//delta = eq8135.solution; // eq. (8.1.35)
			}
			else{
				delta = di/2;
			}

		
			if (pattern==FlowPattern.StratifiedFlow){
<<<<<<< HEAD
				fi = 1+Math.pow(velocitiesRatio, 0.5)*Math.pow((rhoL-rhoG)*g*delta*delta/sigma, 0.25)*(fluidNodeMassFlow/mStrat); // eq. (8.1.41)
			}
			else{
				fi = 1+Math.pow(velocitiesRatio, 0.5)*Math.pow((rhoL-rhoG)*g*delta*delta/sigma, 0.25); // eq. (8.1.40)
=======
				fi = fi*(fluidNodeMassFlow/mStrat); // eq. (8.1.41) //FIXME to check
>>>>>>> branch 'master' of https://github.com/Aqylon/Neodymium.git
			}

			ReL = 4*fluidNodeMassFlow*(1-x)*delta/((1-epsilon)*muL); // eq. (8.1.33)
			PrL = cpL*muL/kL; // eq. (8.1.34)

			alphaC = 0.003*Math.pow(ReL, 0.74)*Math.pow(PrL, 0.5)*kL*fi/delta; // eq. (8.1.32) 

			alphaF = 0.655*Math.pow(rhoL*(rhoL-rhoG)*g*hLd*Math.pow(kL, 3)/(muL*di*q), 1/3); // eq. (8.1.43) //FIXME to check hLd au lieu de HLg ?

			alpha = (alphaF*theta+(2*Math.PI-theta)*alphaC)/(2*Math.PI); // eq. (8.1.23)

			System.out.println("Epsilon : "+std.format(epsilon));
			System.out.println("Thetastrat : "+std.format(Math.toDegrees(airTransferCoefficient))+"°");
			System.out.println("Pattern : "+pattern);
			System.out.println("ReL : "+std.format(ReL));
			System.out.println("PrL : "+std.format(PrL));
			System.out.println("AlphaC : "+std.format(alphaC)+" W.m-2.K-1");
			System.out.println("AlphaF : "+std.format(alphaF)+" W.m-2.K-1");
			System.out.println("Alpha : "+std.format(alpha)+" W.m-2.K-1\n\n");
			
			return alpha;	
		}


		/**
		 * heat transfer for the air side as defined in eq. (8-2-9) 
		 * E. CAO - "Heat transfer in process engineering" - Ch.8 "Finned tubes" - Sec. 8-2-2 "Heat Transfer in Air Coolers" - pp.236-241
		 * @return heat transfer film coefficient for the air side corrected for external fouling and fin efficiency
		 */
		private double computeAirTransferCoefficient(){
			double hf, hfop, hfp; // p for "prime"
			double Rfo = 0.0; // External fouling resistance (m^2.K/W) 

			// Computing hf : heat transfer film coefficient for the air side 
			double as = W*L-ntf*L*(Do+Nm*(Db-Do)*b); // eq. (8-2-17)
			double projectedPermimeter = 2*(Db-Do)*Nm+2*(1-b*Nm); // eq. (8-2-15)
			double De = 2*(Af+Ad)/(Math.PI*projectedPermimeter); // eq. (8-2-14)
			double Res = De*airTotalMassFlow/(as*muAir); // eq. (8-2-16)
			hf = 0.0959*Math.pow(Res, 0.718)*kAir/(De*Math.pow(PrAir, -1/3));

			// Computing hfp : heat transfer film coefficient for the air side corrected for external fouling
			hfp = 1/((1/hf)+Rfo); // eq. (8-2-9)

			// Computing hfop : heat transfer film coefficient for the air side corrected for external fouling and fin efficiency
			double m = Math.sqrt(2*hfp/(kFins*b));
			double H = (Db-Do)/2;
			double Y = (H+b/2)*(1+0.35*Math.log(Db/Do));
			double omega = Math.tanh(m*Y)/(m*Y); // fig. 8-10
			hfop = (omega*Af+Ad)*hfp/Ao; // eq. (8-2-10)

			return hfop;
		}


		private FlowPattern getPattern() throws Exception{

			double fluidNodeMassVelocity = 4*fluidNodeMassFlow/(Math.PI*di*di);
			computeFlowPatternMapBoundaries();

			if(fluidNodeMassVelocity<mStrat || x==0 ) return FlowPattern.StratifiedFlow;
			else if(fluidNodeMassVelocity>mStrat && fluidNodeMassVelocity<mWavy) return FlowPattern.StratifiedWavyFlow;
			else if(fluidNodeMassVelocity>mWavy) return FlowPattern.AnnularFlow;
			else throw new Exception("Flow pattern can not be defined.");
		}

		private void computeFlowPatternMapBoundaries(){

			// Determine the local flow pattern using Thome-El Hajal map
			// Boundary between annular/intermittent & stratified-wavy flow
			double x_mWavyMin = get_x_mWavyMin();

			if(x<x_mWavyMin) mWavy = mWavy(x);
			else mWavy = mWavy(x_mWavyMin);

			// Boundary between stratified-wavy & fully-stratified flow
			double AGd, ALd;
			AGd = A*epsilon/(di*di); // eq. (12.4.21)
			ALd = A*(1-epsilon)/(di*di); // eq. (12.4.20)
			mStrat = Math.pow(266.3*266.3*ALd*AGd*AGd*rhoG*(rhoL-rhoG)*muL*g/(x*x*(1-x)*Math.pow(Math.PI, 3)), 1/3)+20*x; // eq. (12.4.17)
		}

		private double computeEpsilon(double x){
			// Determine epsilon - void fraction
			if(x==0.0) return 0.0;
			if(x==1.0) return 1.0;
			
			double epsilonH, epsilonR;
			epsilonH = 1.0/(1+rhoG*((1.0-x)/x)/rhoL); // eq. (8.1.28)
			epsilonR = (x/rhoG)*1/((1.0+0.12*(1-x))*((x/rhoG)+(1.0-x)/rhoL)+(1.18*(1-x)*Math.pow(g*sigma*(rhoL-rhoG),0.25)/(fluidNodeMassFlow*Math.pow(rhoL, 0.5)))); // eq. (8.1.29)
			
			System.out.println("epsilonH : "+std.format(epsilonH));
			System.out.println("epsilonR : "+std.format(epsilonR));

			return (epsilonH-epsilonR)/Math.log(epsilonH/epsilonR); // eq. (8.1.27)
		}


		private double computeThetaStrat(double epsilon){
			return thetaStrat=2*Math.PI-2*(Math.PI*(1-epsilon)+Math.pow(3*Math.PI/2, 1/3)*(1-2*(1-epsilon)+Math.pow(1-epsilon, 1/3)-Math.pow(epsilon, 1/3))-(1-epsilon)*epsilon*(1-2*(1-epsilon))*(1+4*(Math.pow(1-epsilon, 2)+Math.pow(epsilon, 2)))/200); // eq. (8.1.31)
		}



		private double computeHld(double thetaStrat){
			return 0.5*(1-Math.cos((2*Math.PI-thetaStrat)/2)); // eq.(12.4.22)
		}

		private double mWavy(double x){
			if(x==0){
				x=0.01; //I changed the value of x because I had troubles with infinite values for x=0. 
			}

			double qDNB; // Heat flux of departure from nucleate boiling

			double WeOnFrL, F1, F2, AGd;
			double epsilonLocal = computeEpsilon(x);
			double thetaStratLocal = computeThetaStrat(epsilonLocal);

			double hLdLocal = computeHld( thetaStratLocal);
			AGd = A*epsilonLocal/(di*di); // eq. (12.4.21)
			WeOnFrL = g*di*di*rhoL/sigma; //eq. (12.4.6)

			double hLG=this.fluidNodeInletState.saturatedVaporState.enthalpy-this.fluidNodeInletState.saturatedLiquidState.enthalpy ;//FIXME
			qDNB = 0.131*Math.pow(rhoG, 0.5)*hLG*Math.pow(g*(rhoL-rhoG)*sigma, 0.25) ;// eq. (12.4.9)


			F1 = 646.0*Math.pow(q/qDNB, 2)+64.8*(q/qDNB); // eq. (12.4.8a)
			F2 = 18.8*(q/qDNB)+1.023; // eq. (12.4.8b)

			double b=-75*Math.exp(-Math.pow(x*x-0.97, 2)/(x*(1-x)));
			double c=((Math.PI*Math.PI*Math.pow(1-x, -F1)/(25*hLdLocal*hLdLocal))*Math.pow(WeOnFrL, -F2)+1);
			return Math.pow(16*Math.pow(AGd, 3)*g*di*rhoL*rhoG*c/(x*x*Math.PI*Math.PI*Math.pow(1-Math.pow(2*hLdLocal-1, 2), 0.5)), 0.5)+50+b; // eq. (12.4.1/18)

		}


		/**
		 * Find the minimum of the mWawy function using a trial and error method.
		 * H.H. Rosenbrock, « An automatic method for finding the greatest or least value of a function », Comput. J., 3, 175 (1960)
		 * @return
		 */
		private double get_x_mWavyMin(){
			double precision = 1.0E-5 ;
			double step=0.1 , dilatation=3 , contraction=0.4;

			double currentPoint=0.6 ;
			double mWavyCurrent=mWavy(currentPoint);
			double trialPoint=currentPoint+step;
			double mWavyTrial=mWavy(trialPoint);   

			while(Math.abs( (mWavyCurrent-mWavyTrial)/mWavyCurrent )>precision){
				if(mWavyTrial<mWavyCurrent){
					currentPoint=trialPoint;
					mWavyCurrent=mWavyTrial;

					step=dilatation*step;
					trialPoint=currentPoint+step;
					trialPoint=Math.min(trialPoint,0.99);trialPoint=Math.max(trialPoint,0.01);
					mWavyTrial=mWavy(trialPoint);
				}
				else{
					step=-contraction*step;
					trialPoint=currentPoint+step;
					trialPoint=Math.min(trialPoint,0.99);trialPoint=Math.max(trialPoint,0.01);
					mWavyTrial=mWavy(trialPoint);  
				}
			}
			return currentPoint;
		}


		public double getAirNodeOutletDeltaT(){
			return airNodeOutletDeltaT;
		}

		/*    public Hashtable getDebugPhysicalProperties(){

    }*/

	}

}
