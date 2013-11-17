package com.aqylon.aircooledCondenser;

import com.aqylon.thermodynamics.physics.ThState;
import com.aqylon.thermodynamics.physics.ThFluid;
import com.aqylon.utilities.quadraticEquation.quadraticEquation;


/**
 * 
 * @author Vincent Lamblot
 *
 */
public class HeatTransferLocalUnit {

	@SuppressWarnings("unused")
	private AircooledCondenser condenser;

	private double airNodeOutletDeltaT, fluidNodeMassFlow;
	private ThState fluidNodeInletState, fluidNodeOutletState ;


	public double totalTransferCoefficient, enthalpyTransfered, fluidTransferCoefficient, airTransferCoefficient;
	private double mWavy, mStrat, hLd, epsilon, thetaStrat, A; 
	private double x ,rhoL, rhoG ,muL , cpL , sigma, kL ;


	/*
	 * Parameters inherited from condenser
	 */
	private double dA, cpAir, g, 
	di, W, L, ntf, Do, Nm, Af, b, Db, Ad, airNodeMassFlow, muAir, 
	kAir, PrAir, kFins, Ao;

	/**
	 * 
	 * @param condenser
	 * @param airNodeInletTemperature
	 * @param fluidNodeInletState
	 */
	public HeatTransferLocalUnit(AircooledCondenser condenser, double airNodeInletTemperature, ThState fluidNodeInletState, double fluidNodeMassFlow){

		this.fluidNodeInletState = fluidNodeInletState;
		this.fluidNodeMassFlow = fluidNodeMassFlow;
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

		airNodeOutletDeltaT = enthalpyTransfered/cpAir;

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
		// Determine epsilon - void fraction
		double epsilon, epsilonH, epsilonR;
		epsilonH = 1/(1+rhoG*((1-x)/x)/rhoL); // eq. (8.1.28)
		epsilonR = (x/rhoG)*1/((1+0.12*(1-x))*((x/rhoG)+(1-x)/rhoL)+(1.18*(1-x)*Math.pow(g*sigma*(rhoL-rhoG),0.25)/(fluidNodeMassFlow*Math.pow(rhoL, 0.5))));
		epsilon = (epsilonH-epsilonR)/Math.log(epsilonH/epsilonR); // eq. (8.1.27)
		this.epsilon = epsilon;

		double theta, thetaStrat; 
		double ReL, PrL, delta, fi, velocitiesRatio, q=enthalpyTransfered;
		double alfaC=0,alfaF=0,alfa;
		double AL;

		// Determine the local flow pattern using Thome-El Hajal map
		thetaStrat=2*Math.PI-2*(Math.PI*(1-epsilon)+Math.pow(3*Math.PI/2, 1/3)*(1-2*(1-epsilon)+Math.pow(1-epsilon, 1/3)-Math.pow(epsilon, 1/3))-(1-epsilon)*epsilon*(1-2*(1-epsilon))*(1+4*(Math.pow(1-epsilon, 2)+Math.pow(epsilon, 2)))/200); // eq. (8.1.31)
		this.thetaStrat = thetaStrat;
		AL = Math.pow(di, 2)*((2*Math.PI-thetaStrat)-Math.sin(2*Math.PI-thetaStrat))/8; // eq. (8.1.24)
		this.A = AL/(1-epsilon);
		
		FlowPattern pattern = determineFlowPattern();

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




	/*
	 * Utilities
	 */


	public double getAirNodeOutletDeltaT(){
		return airNodeOutletDeltaT;
	}
}
