package com.aqylon.aircooledCondenser;

public class FlowPatternMap {
  
  public FlowPatternMap(){
    
  }
  
  public FlowPattern findPattern(){
    
  //TODO to check by FD

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
    // Defining flow pattern
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
  
  public void printMap(){
    
    
  }
  
  
  
  /**
   * 
   * @param x
   * @return
   */
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
  
}
