package com.aqylon.aircooledCondenser;

import com.aqylon.thermodynamics.physics.ThFluid;
import com.aqylon.thermodynamics.physics.ThState;
import java.util.Hashtable;

public class NumericalSolver {

  /**
   * Geometry of the condenser
   */
  private AircooledCondenser condenser;

  /**
   * Working data arrays
   */
  private double[] airTemperatures;
  private double[] deltaT;
  private ThState[] innerFluidStates;
  
  /**
   * Discretisation parameters and array sizes
   */
  private int nPipe ; //number of pipes in a vertical slice (count one for two pipes linked in one pass) ; nPipe=nt
  private int N ; //number of discretisation steps for the condenser length 
  private int nPass ; //number of pass nPass=np
  private int nf ; //nPipe* nPass
  
  /**
   * Options to control the class behavior  
   */
  String typeEnergyBalance ; // energy terms to take into account in the energy balance in function computeStateOfNextFluidNode
  String debugMode ; // to control if you want to get lots of physical properties in arrays for debug
  
  private Hashtable debugPhysicalProperties ; // lots of physical properties in arrays for debug
  
  /**
   * Inlet state of the fluid
   */
  private ThState fluidInletState;
  private double fluidNodeMassFlow ;
  
  /**
  * FlowPatternMap of the pipes
  */
  private FlowPatternMap flowPatternMap ;
 
  /**
   * Grid to visualize computed properties with paraview
   */
  private RectilinearGrid grid;
  
  
  public NumericalSolver(AircooledCondenser condenser, ThState fluidInletState, double fluidTotalMassFlow){
    
    typeEnergyBalance="NoConvection";
    
    nPipe=condenser.nt;
    nf=condenser.nf;
    N=condenser.N;
    nPass=condenser.np;
    
    this.fluidInletState=fluidInletState; 
    fluidNodeMassFlow=fluidTotalMassFlow/ nPipe;
    
    deltaT = new double[N*nf]; // In K
    innerFluidStates = new ThState[(N+1)*nf];
    airTemperatures = new double[N*(nf+1)]; // In K   
    
    for(int i=0; i < airTemperatures.length; i++){
      airTemperatures[i] = condenser.airInletTemperature;
    }
    
    this.flowPatternMap = new FlowPatternMap(   ) ;
  }
  
    
  public ThState[] computeOutletStateOfEachPipe(){
    
    computeAllStatesWithConvergence();
    
    ThState[] statesToMix= new ThState[nPipe];
    
    for(int iPipe=0 ; iPipe < nPipe ; iPipe++){
      statesToMix[iPipe] = this.innerFluidStates[ getPipeOutletNode(iPipe) ];
    }
    
    return statesToMix ;
  }
  
  
  public void computeAllStatesWithConvergence(){
    
    boolean converged=false ;
    
    while(!converged){
      
      computeAllInnerStatesAlongPipes();
      
      //Calcul des critères de convergence  
      converged=true ;                                // TO DO
      
      updateAirTemperatureMap();

    }
  }
  
  
  private void computeAllInnerStatesAlongPipes(){
        
    for(int iPipe=0; iPipe < nPipe; iPipe++){
      for(int iPass=0 ; iPass < nPass ; iPass++ ){
        for(int iNode=0 ; iNode < N  ; iNode ++){
        
          int currentFluidNode = getFluidNode(iPipe, iPass, iNode);
          int currentDeltaTNode = getDeltaTNode(iPipe, iPass, iNode);
          int currentAirNode = getAirNode(iPipe, iPass, iNode);
          
          HeatTransferLocalUnit currentTransferUnit = new HeatTransferLocalUnit(this.condenser, this.flowPatternMap, this.airTemperatures[currentAirNode], this.innerFluidStates[currentFluidNode], fluidNodeMassFlow);
          this.deltaT[currentDeltaTNode] = currentTransferUnit.getAirNodeOutletDeltaT();
          
          ThState  fluidNodeOutletState = computeStateOfNextFluidNode(currentFluidNode, currentTransferUnit);
        
          int nextFluidNode = getFluidNode(iPipe, iPass, iNode+1);
          this.innerFluidStates[nextFluidNode] = fluidNodeOutletState ;
        }
        
        // transition from one pass to the other
        int previousPassFluidNode= getFluidNode(iPipe, iPass, N);
        int nextPassFluidNode= getFluidNode(iPipe, iPass+1, 0);
        this.innerFluidStates[nextPassFluidNode]=this.innerFluidStates[previousPassFluidNode];
      }
    }
  }

  
  private ThState computeStateOfNextFluidNode(int currentFluidNode, HeatTransferLocalUnit currentTransferUnit){
        
    ThState currentState=this.innerFluidStates[currentFluidNode];  
    
    ThState fluidNodeOutletState = new ThState(currentState.fluid);
    
    if(typeEnergyBalance.equals("NoConvection")){
      fluidNodeOutletState.pressure = currentState.pressure;
      double enthalpy=currentState.enthalpy + currentTransferUnit.enthalpyTransfered ;
      fluidNodeOutletState.setEnthalpy(enthalpy);
    }
    
     
    return fluidNodeOutletState ;
  }
  

  private void updateAirTemperatureMap(){
    for(int iAirFlow=0 ; iAirFlow<N ; iAirFlow ++){
      
      for(int iNode=0 ; iNode < nf ; iNode++){
        
        int currentAirNode=getAirFlowNode(iAirFlow,iNode);
        int nextAirNode=getAirFlowNode(iAirFlow,iNode+1);
        int currentDeltaTNode=getAirFlowDeltaTNode(iAirFlow,iNode);
                
        airTemperatures[nextAirNode]=airTemperatures[currentAirNode]+this.deltaT[currentDeltaTNode];
      }
    }
  }
  
  
  /* 
   * Functions to access properties inside arrays...
   */
  
  /*
   * ...first when you are counting nodes along pipes
   */
  
  private int getFluidNode(int iPipe, int iPass, int iNode){
    
    int node = iPipe*(N+1)+iPass*(nPipe*(N+1)) ;  //last column is outside the condenser on the right 
    
    if(iPass%2==1){
      node=node+iNode;
    }
    else{
      node=node+N-iNode;
    }
    return node ;
  }
  
  private int getAirNode(int iPipe, int iPass, int iNode){
    
    int node = N + iPipe*N+iPass*(nPipe*N) ;//first row is upside the condenser, it's the air output
    
    if(iPass%2==1){
      node=node+iNode;
    }
    else{
      node=node+N-iNode;
    }
    return node;
  }
  
  private int getDeltaTNode(int iPipe, int iPass, int iNode){
    
    int node = iPipe*N+iPass*(nPipe*N) ;//first row is upside the condenser, it's the air output
    
    if(iPass%2==1){
      node=node+iNode;
    }
    else{
      node=node+N-iNode;
    }
    return node ;
  }  
  
  
  /*
   * ... second when you are counting nodes along air flow  ...
   */
  
  private int getAirFlowNode(int iAirFlow, int iNode){
        
    int node = (nf+1-iNode)*N + iAirFlow ;
    return node;
  }  
  
  
  private int getAirFlowDeltaTNode(int iAirFlow, int iNode){
    
    int node = (nf-iNode)*N + iAirFlow ;
    return node;
  }
  
  /*
   * ... then when you're looking for outlet states.
   */
  
  private int getPipeOutletNode(int iPipe){
    
    int node=(nPass-1)*nPipe*(N+1);
    
    if(nPass%2==0){
      node=node+iPipe*(N+1);
    }
    else{
      node=node+(iPipe+1)*(N+1)-1; 
    }
    return node ;
  }
  

  
}
