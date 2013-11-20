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
		
		
		AircooledCondenser condenseur =new AircooledCondenser(1, 1,1, 1,1,1, 1,2, 300);
		
		ThFluid fluid=ThFluid.Water;
		ThState fluidInletState=new ThState(fluid);
		fluidInletState=fluidInletState.create();
		ThState outlet=condenseur.computeOutlet( fluidInletState, 4);

	}

}
 