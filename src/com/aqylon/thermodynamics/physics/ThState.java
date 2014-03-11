package com.aqylon.thermodynamics.physics;


/**
 * 
 * @author Pierre Convert & Cie
 * 
 * properties are accessed in two different ways :
 *    - direct access to the (public) fields for fluid properties calculation engine
 *    - through getter/setter for use in cycle modeling
 */

public class ThState{

	public boolean isEnthalpyDefined = false, isQualityDefined = false, isStateDefined = false;
	
	public ThFluid fluid;

	/**
	 * Phase (one of SubcooledLiquid, SaturatedLiquid, Saturated, SaturatedVapor, SuperheatedVapor, Critical, SuperCritical)
	 */
	public ThPhase phase = ThPhase.Saturated;

	/**
	 * Temperature (in K)
	 */
	public double temperature;

	/**
	 * Pressure (in Pa)
	 */
	public double pressure;

	/**
	 * Vapor quality (in %)
	 */
	public double quality;

	/**
	 * Density (in kg.m-3)
	 */
	public double density;

	/**
	 * Enthalpy (in J.kg-1)
	 */
	public double enthalpy;

	/**
	 * Isobaric Specific Heat Capacity (J/K)
	 */
	public double isobaricHeatCapacity;

	/**
	 * Thermal conductivity (W/(m.K))
	 */
	public double thermalConductivity;

	/**
	 * Dynamic Viscosity (Pa/s)
	 */
	public double dynamicViscosity;

	/**
	 * Surface tension (N/m)
	 */
	public double surfaceTension;
	
	/**
	 * Saturated liquid state sharing (T,p) with this state
	 */
	public ThState saturatedLiquidState;

	/**
	 * Saturated vapor state sharing (T,p) with this state
	 */
	public ThState saturatedVaporState;

	public ThState(ThFluid fluid){
		this.fluid = fluid;

	}
	
	
	public ThState create(){
		
		ThState s = new ThState(fluid);
		
		switch(fluid){
		
		case R245fa:
			s.temperature = 36.0 + 273.15;
			s.pressure = 0.9894E5;
			s.saturatedLiquidState = new ThState(fluid);
			s.saturatedVaporState = new ThState(fluid);

			// Saturated liquid state properties
			s.saturatedLiquidState.temperature = s.temperature;
			s.saturatedLiquidState.pressure = s.pressure;
			s.saturatedLiquidState.phase = ThPhase.SaturatedLiquid;
			s.saturatedLiquidState.quality = 0.0;
			s.saturatedLiquidState.density = 1.343E3;
			s.saturatedLiquidState.enthalpy = 234.96E3;
			s.saturatedLiquidState.isobaricHeatCapacity = 1.0473E3;
			s.saturatedLiquidState.thermalConductivity = 0.08;
			s.saturatedLiquidState.dynamicViscosity = 4.57E-4;
			s.saturatedLiquidState.surfaceTension = 11.773e-3;

			// Saturated vapor state properties
			s.saturatedVaporState.temperature = s.temperature;
			s.saturatedVaporState.pressure = s.pressure;
			s.saturatedLiquidState.phase = ThPhase.SaturatedVapor;
			s.saturatedVaporState.quality = 1.0;
			s.saturatedVaporState.density = 7.49;
			s.saturatedVaporState.enthalpy = 366.41E3;
			s.saturatedVaporState.isobaricHeatCapacity = 0.8979E3;
			s.saturatedVaporState.thermalConductivity = 0.012;
			s.saturatedVaporState.dynamicViscosity = 1.14E-5;
			break;
			
		case SES36:
			s.temperature = 36.0 + 273.15;
			s.pressure = 0.9894E5;
			s.saturatedLiquidState = new ThState(fluid);
			s.saturatedVaporState = new ThState(fluid);

			// Saturated liquid state properties
			s.saturatedLiquidState.temperature = s.temperature;
			s.saturatedLiquidState.pressure = s.pressure;
			s.saturatedLiquidState.phase = ThPhase.SaturatedLiquid;
			s.saturatedLiquidState.quality = 0.0;
			s.saturatedLiquidState.density = 1.343E3;
			s.saturatedLiquidState.enthalpy = 234.96E3;
			s.saturatedLiquidState.isobaricHeatCapacity = 1.0473E3;
			s.saturatedLiquidState.thermalConductivity = 0.08;
			s.saturatedLiquidState.dynamicViscosity = 4.57E-4;
			s.saturatedLiquidState.surfaceTension = 11.773e-3;

			// Saturated vapor state properties
			s.saturatedVaporState.temperature = s.temperature;
			s.saturatedVaporState.pressure = s.pressure;
			s.saturatedLiquidState.phase = ThPhase.SaturatedVapor;
			s.saturatedVaporState.quality = 1.0;
			s.saturatedVaporState.density = 7.49;
			s.saturatedVaporState.enthalpy = 366.41E3;
			s.saturatedVaporState.isobaricHeatCapacity = 0.8979E3;
			s.saturatedVaporState.thermalConductivity = 0.012;
			s.saturatedVaporState.dynamicViscosity = 1.14E-5;
			break;
			
		case Water: // source : http://www.peacesoftware.de/einigewerte/calc_dampf.php5
			s.temperature = 100.0 + 273.15;
			s.pressure = 0.9894E5;
			s.saturatedLiquidState = new ThState(fluid);
			s.saturatedVaporState = new ThState(fluid);

			// Saturated liquid state properties
			s.saturatedLiquidState.temperature = s.temperature;
			s.saturatedLiquidState.pressure = s.pressure;
			s.saturatedLiquidState.phase = ThPhase.SaturatedLiquid;
			s.saturatedLiquidState.quality = 0.0;
			s.saturatedLiquidState.density = 1.343E3;
			s.saturatedLiquidState.enthalpy = 419E3;
			s.saturatedLiquidState.isobaricHeatCapacity = 1.0473E3;
			s.saturatedLiquidState.thermalConductivity = 0.08;
			s.saturatedLiquidState.dynamicViscosity = 4.57E-4;
			s.saturatedLiquidState.surfaceTension = 11.773e-3;

			// Saturated vapor state properties
			s.saturatedVaporState.temperature = s.temperature;
			s.saturatedVaporState.pressure = s.pressure;
			s.saturatedVaporState.phase = ThPhase.SaturatedVapor;
			s.saturatedVaporState.quality = 1.0;
			s.saturatedVaporState.density = 7.49;
			s.saturatedVaporState.enthalpy = 2676E3;
			s.saturatedVaporState.isobaricHeatCapacity = 0.8979E3;
			s.saturatedVaporState.thermalConductivity = 0.012;
			s.saturatedVaporState.dynamicViscosity = 1.14E-5;
			break;
			
		default:
			break;
		}
		return s;	
	}
	
	
	
	/*
	 * Setters
	 */
	
	/**
	 * 
	 * @param enthalpy
	 */
	public void setEnthalpy(double enthalpy){
		if(isEnthalpyDefined){
			throw new RuntimeException("Enthalpy has already been defined");
		}
		if(isStateDefined){
			throw new RuntimeException("State is already fully defined");
		}
		
		this.enthalpy = enthalpy;
		isEnthalpyDefined = true;
		
		double Q, h, hSatL,hSatV;

		h=enthalpy;
		//Enthalpy (in J.kg-1)
		hSatL=saturatedLiquidState.enthalpy;
		hSatV=saturatedVaporState.enthalpy;
		//Vapor quality (in %)
		Q=(h-hSatL)/(hSatV-hSatL);
		quality=Q;
		
		isQualityDefined = true;
		isStateDefined = true;
	}
	
	/**
	 * 
	 * @param quality
	 */
	public void setQuality(double quality){
		if(isEnthalpyDefined){
			throw new RuntimeException("Quality has already been defined");
		}
		if(isStateDefined){
			throw new RuntimeException("State is already fully defined");
		}
		
		this.quality = quality;
		isEnthalpyDefined = true;
		
		double Q, h, hSatL,hSatV;
		
		Q=quality;
		//Enthalpy (in J.kg-1)
		hSatL=saturatedLiquidState.enthalpy;
		hSatV=saturatedVaporState.enthalpy;
		h=hSatL+Q*(hSatV-hSatL);
		enthalpy=h;
		
		isEnthalpyDefined = true;
		isStateDefined = true;
	}
	
}
