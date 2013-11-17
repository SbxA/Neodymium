package com.aqylon.aircooledCondenser;
import java.util.Enumeration;
import java.util.Hashtable;

public class RectilinearGrid {

  int nodeNumber ;
  Hashtable  physicalProperties ;
  
  
  public RectilinearGrid(int nodeNumber){
    this.nodeNumber=nodeNumber;
    this.physicalProperties= new Hashtable();
  }
  
  public void addField(double[] property, String name){
    assert property.length == nodeNumber :
      "This property doesn't have the same length as the RectilinearGrid.";
    physicalProperties.put(name,property);
  }
    
  public void writeVTKFile(){
    Enumeration e = physicalProperties.elements();
    
    while(e.hasMoreElements()){
      e.nextElement();
    }
  }
  
  
}
