package com.aqylon.aircooledCondenser;

import com.aqylon.aircooledCondenser.AircooledCondenser.LocalHeatTransferUnit;
import com.aqylon.thermodynamics.physics.ThFluid;
import com.aqylon.thermodynamics.physics.ThState;

public class Test {

    public static void main(String[] args) throws Exception {
        final double di = 20.8E-3;  //di internal diameter of tube (m)
        final double Do = 26.7E-3;  //Do external diameter of tube without fins (m) 
        final double Db = Do*1.5 ;  //Db external diameter of tube with fins (m) // ask
        final double b = 0.5E-3;   //b fin width (m) // ask
        final double fs = 2.4E-3;  //fs fin spacing (m) 
        final double L = 10.2;  //L condenser tubes length (m)
        final double W = Math.sqrt(2.5*2.5+2.3*2.3/4);  //W condenser width (m)
        final int ntf = 168/2; // ntf number of tube horizontal rows // ask
        final int np = 2; // np number of passes
        final int nt = 2; // nt number of tubes per pass
        final double airTotalMassFlow = 354854*1.2/3600 ;//airTotalMassFlow air mass flow at inlet (kg/s)
        final double airInletTemperature = 20+273;//airInletTemperature air temperature at inlet (K)
        final double fluidTotalMassFlow = 20.;

        AircooledCondenser condenser = new AircooledCondenser(di, Do, Db, b, fs, L, W,ntf,np,nt,airTotalMassFlow, airInletTemperature);
        condenser.fluidNodeMassFlow = fluidTotalMassFlow/(nt*ntf);

        ThFluid fluid = ThFluid.Water;
        ThState inletState = new ThState(fluid);
        inletState = inletState.create();
        inletState.setQuality(0.8);
        
        LocalHeatTransferUnit unit = condenser.new LocalHeatTransferUnit(airInletTemperature, inletState);        
    
        unit.computeTransfer();

        
    }

}