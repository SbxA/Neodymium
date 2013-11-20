package com.aqylon.aircooledCondenser;

import com.aqylon.thermodynamics.physics.ThState;
import com.aqylon.thermodynamics.physics.ThFluid;
import com.aqylon.utilities.quadraticEquation.quadraticEquation;

import java.util.Hashtable;

/**
 * 
 * @author Vincent Lamblot & Tristan Agaësse
 *
 */
public class HeatTransferLocalUnit {

	@SuppressWarnings("unused")
	private AircooledCondenser condenser;

	private double airNodeOutletDeltaT, airNodeInletTemperature, fluidNodeMassFlow;
	private ThState fluidNodeInletState ;


	public double totalTransferCoefficient, enthalpyTransfered, fluidTransferCoefficient, airTransferCoefficient;
	public double mWavy, mStrat, hLd, epsilon, thetaStrat, A, theta; 
	private double x ,rhoL, rhoG ,muL , cpL , sigma, kL ;

	/*
	 * Parameters inherited from condenser
	 */
	private double dA, cpAir, g, 
	di, W, L, ntf, Do, Nm, Af, b, Db, Ad, airNodeMassFlow, muAir, 
	kAir, PrAir, kFins, Ao;

  double ReL, PrL, delta, fi, velocitiesRatio, q=enthalpyTransfered;
  double alfaC=0,alfaF=0,alfa;
  double AL;
	
	/**
	 * 
	 * @param condenser
	 * @param airNodeInletTemperature
	 * @param fluidNodeInletState
	 */
	public HeatTransferLocalUnit(AircooledCondenser condenser, double airNodeInletTemperature, ThState fluidNodeInletState, double fluidNodeMassFlow){

		this.fluidNodeInletState = fluidNodeInletState;
		this.fluidNodeMassFlow = fluidNodeMassFlow;
		this.airNodeInletTemperature=airNodeInletTemperature;
		    
		this.condenser = condenser;
		
		dA = condenser.dA;
		cpAir = AircooledCondenser.cpAir;
		g = AircooledCondenser.g;
		di = condenser.di;
		W = condenser.W;
		L = condenser.L;
		ntf = condenser.ntf;
		Do = condenser.Do;
		Nm= condenser.Nm;
		Af = condenser.Af;
		b = condenser.b;
		Db= condenser.Db;
		Ad = condenser.Ad;
		airNodeMassFlow = condenser.airNodeMassFlow;
		muAir = AircooledCondenser.muAir; 
		kAir = AircooledCondenser.kAir;
		PrAir = AircooledCondenser.PrAir;
		kFins = AircooledCondenser.kFins;
		Ao = condenser.Ao;


		x = fluidNodeInletState.quality;
		rhoL = fluidNodeInletState.saturatedLiquidState.density;
		rhoG = fluidNodeInletState.saturatedVaporState.density;
		muL = fluidNodeInletState.saturatedLiquidState.dynamicViscosity;
		cpL = fluidNodeInletState.saturatedLiquidState.isobaricHeatCapacity;
		kL = fluidNodeInletState.saturatedLiquidState.thermalConductivity;
		sigma = fluidNodeInletState.saturatedLiquidState.surfaceTension;

	}

	
	public void computeTransfer(){
	  
	   airTransferCoefficient = computeAirTransferCoefficient();

	    // Start automatically computation

	    this.enthalpyTransfered = 1000.0; // In J, 1st hypothesis
	    double enthalpyTransferedBuffer = 0.0;

	    while(Math.abs(enthalpyTransfered-enthalpyTransferedBuffer)/enthalpyTransfered > 1.0E-5){
	      enthalpyTransferedBuffer = enthalpyTransfered;
	      fluidTransferCoefficient = computeFluidTransferCoefficient();
	      totalTransferCoefficient = 1/(1/fluidTransferCoefficient+1/airTransferCoefficient);
	      enthalpyTransfered = totalTransferCoefficient*dA*(fluidNodeInletState.temperature-airNodeInletTemperature); // Chosen as positive in the fluid to air way. 
	    }

	    airNodeOutletDeltaT = enthalpyTransfered/(cpAir*airNodeMassFlow);
	  
	}
	
	
	
	/**
	 * hio' as defined in eq. (8-2-11)
	 * E. CAO - "Heat transfer in process engineering" - Ch.8 "Finned tubes" - Sec. 8-2-2 "Heat Transfer in Air Coolers" - pp.236-241
	 * Computing based on the Thome-El Hajal-Calvini flow pattern method
	 * Wolverine Tube, INC. - Ch.8.1 "Condensation inside Horizontal Tubes" - "Heat transfer model" pp.8.8-8.12
	 * @param fluidNodeInletState fluid thermodynamical properties at node inlet
	 * @return hio'
	 */
	private double computeFluidTransferCoefficient(){
		
		
		computeFlowPatternMapBoundaries();
		
		
		FlowPattern pattern = findPattern();

		// Compute hio = alfa
		switch(pattern){

		case AnnularFlow:
		  theta = 0.0;
		  break;

		case StratifiedFlow:
		  theta = thetaStrat*Math.pow((mWavy-fluidNodeMassFlow)/(mWavy-mStrat), 0.5); // eq. (8.1.30)
		  break;

		case StratifiedWavyFlow:
		  theta = thetaStrat;
		  break;

		default:
		  throw new RuntimeException("Flow pattern not defined. Do no know how to solve this case !!");
		}

		
		velocitiesRatio = rhoL*x*(1-epsilon)/(rhoG*epsilon*(1-x)); // eq. (8.1.36/37)

		if(epsilon>0.5){
			quadraticEquation eq8135 = new quadraticEquation(4.0, -4*di*di, 8*AL/(2*Math.PI-theta)+Math.pow(di, 4)-di*di); 
			delta = eq8135.solution; // eq. (8.1.35)
		}
		else{
			delta = di/2;
		}

		fi = 1+Math.pow(velocitiesRatio, 0.5)*Math.pow((rhoL-rhoG)*g*delta*delta/sigma, 0.25); // eq. (8.1.40)

		if (pattern==FlowPattern.StratifiedFlow){
			fi = fi*(fluidNodeMassFlow/mStrat); // eq. (8.1.41)
		}

		ReL = 4*fluidNodeMassFlow*(1-x)*delta/((1-epsilon)*muL); // eq. (8.1.33)
		PrL = cpL*muL/kL; // eq. (8.1.34)

		alfaC = 0.003*Math.pow(ReL, 0.74)*Math.pow(PrL, 0.5)*kL*fi/delta; // eq. (8.1.32)

		alfaF = 0.655*Math.pow(rhoL*(rhoL-rhoG)*g*hLd*Math.pow(kL, 3)/(muL*di*q), 1/3); // eq. (8.1.43)

		alfa = (alfaF*theta+(2*Math.PI-theta)*alfaC)/2*Math.PI; // eq. (8.1.23)

		return alfa;
	}


	/**
	 * hfo' as defined in eq. (8-2-9) 
	 * E. CAO - "Heat transfer in process engineering" - Ch.8 "Finned tubes" - Sec. 8-2-2 "Heat Transfer in Air Coolers" - pp.236-241
	 * @return hfo'
	 */
	private double computeAirTransferCoefficient(){
		double hf, hfop, hfp; // p for '
		double Rfo = 0.0; // External fouling resistance (m^2.K/W) 

		// Computing hf
		double as = W*L-ntf*L*(Do+Nm*(Db-Do)*b); // eq. (8-2-17)
		double projectedPermimeter = 2*(Db-Do)*Nm+2*(1-b*Nm); // eq. (8-2-15)
		double De = 2*(Af+Ad)/(Math.PI*projectedPermimeter); // eq. (8-2-14)
		double Res = De*airNodeMassFlow/(as*muAir); // eq. (8-2-16)
		hf = 0.0959*Math.pow(Res, 0.718)*kAir/(De*Math.pow(PrAir, 1/3));

		// Computing hfp
		hfp = 1/((1/hf)+Rfo); // eq. (8-2-9)

		// Computing hfop
		double m = Math.sqrt(2*hfp/(kFins*b));
		double H = (Db-Do)/2;
		double Y = (H+b/2)*(1+0.35*Math.log(Db/Do));
		double omega = Math.tan(m*Y)/(m*Y); // fig. 8-10
		hfop = (omega*Af+Ad)*hf/Ao; // eq. (8-2-10)

		return hfop;
	}


	public FlowPattern findPattern(){
	  
    double fluidNodeMassVelocity = 4*fluidNodeMassFlow/(Math.PI*di*di);
	  
    if(fluidNodeMassVelocity<mStrat){
      return FlowPattern.StratifiedFlow;

    }
    else if(fluidNodeMassVelocity>mStrat && fluidNodeMassVelocity<mWavy){
      return FlowPattern.StratifiedWavyFlow;

    }
    else if(fluidNodeMassVelocity>mWavy){
      return FlowPattern.AnnularFlow;
    }
    else{
      throw new RuntimeException("Flow pattern can not be defined.");
    }
	  
	}
	
	
	public void computeFlowPatternMapBoundaries(){
	  
	  // Determine epsilon - void fraction
    double epsilonH, epsilonR;
    epsilonH = 1/(1+rhoG*((1-x)/x)/rhoL); // eq. (8.1.28)
    epsilonR = (x/rhoG)*1/((1+0.12*(1-x))*((x/rhoG)+(1-x)/rhoL)+(1.18*(1-x)*Math.pow(g*sigma*(rhoL-rhoG),0.25)/(fluidNodeMassFlow*Math.pow(rhoL, 0.5)))); // eq. (8.1.29)
    this.epsilon = (epsilonH-epsilonR)/Math.log(epsilonH/epsilonR); // eq. (8.1.27)



    // Determine the local flow pattern using Thome-El Hajal map
    this.thetaStrat=2*Math.PI-2*(Math.PI*(1-epsilon)+Math.pow(3*Math.PI/2, 1/3)*(1-2*(1-epsilon)+Math.pow(1-epsilon, 1/3)-Math.pow(epsilon, 1/3))-(1-epsilon)*epsilon*(1-2*(1-epsilon))*(1+4*(Math.pow(1-epsilon, 2)+Math.pow(epsilon, 2)))/200); // eq. (8.1.31)
    AL = Math.pow(di, 2)*((2*Math.PI-thetaStrat)-Math.sin(2*Math.PI-thetaStrat))/8; // eq. (8.1.24)
    this.A = AL/(1-epsilon);
    
    
    // Boundary between annular/intermittent & stratified-wavy flow
    double x_mWavyMin = get_x_mWavyMin();
    
    if(x<x_mWavyMin){
      mWavy = mWavy(x);
    }
    else{
      mWavy = mWavy(x_mWavyMin);
    }


    // Boundary between stratified-wavy & fully-stratified flow
    double AGd, ALd;
    AGd = A*epsilon/(di*di); // eq. (12.4.21)
    ALd = A*(1-epsilon)/(di*di); // eq. (12.4.20)
    mStrat = Math.pow(266.3*266.3*ALd*AGd*AGd*rhoG*(rhoL-rhoG)*muL*g/(x*x*(1-x)*Math.pow(Math.PI, 3)), 1/3)+20*x; // eq. (12.4.17)
    
    
	}
	
  private double mWavy(double x){
    double WeOnFrL, qDNB, F1, F2, AGd, q=enthalpyTransfered;

    hLd = 0.5*(1-Math.cos((2*Math.PI-thetaStrat)/2)); // eq.(12.4.22)
    AGd = A*epsilon/(di*di); // eq. (12.4.21)
    WeOnFrL = g*di*di*rhoL/sigma; //eq. (12.4.6)
    qDNB = 0.131*Math.pow(rhoG, 0.5)*hLd*Math.pow(g*(rhoL-rhoG)*sigma, 0.25); // eq. (12.4.9)
    F1 = 646.0*Math.pow(q/qDNB, 2)+64.8*(q/qDNB); // eq. (12.4.8a)
    F2 = 18.8*(q/qDNB)+1.023; // eq. (12.4.8b)

    return Math.pow(16*Math.pow(AGd, 3)*g*di*rhoL*rhoG*((Math.PI*Math.PI*Math.pow(1-x, -F1)/(25*hLd*hLd))*Math.pow(WeOnFrL, -F2)+1)/(x*x*Math.PI*Math.PI*Math.pow(1-Math.pow(2*hLd-1, 2), 0.5)), 0.5)+50-75*Math.exp(-Math.pow(x*x-0.97, 2)/(x*(1-x))); // eq. (12.4.1/18)
  }

  /**
   * 
   * @param x
   * @return
   */
  private double mWavy_derivative(double x){
    double dx = 1.0E-5;

    return (mWavy(x+dx)-mWavy(x))/dx;
  }

  /**
   * 
   * @return
   */
  private double get_x_mWavyMin(){
    double x = 0.1, error = 1.0E-5, alfa = 0.5;

    while(mWavy_derivative(x)>error){
      x = x - alfa*mWavy_derivative(x);
    }
    return x;
  }


	public double getAirNodeOutletDeltaT(){
		return airNodeOutletDeltaT;
	}
	
	public Hashtable getDebugPhysicalProperties(){
	  
	}
	
}
