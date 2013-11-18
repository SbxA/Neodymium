package com.aqylon.aircooledCondenser;
import java.util.Enumeration;
import java.util.Hashtable;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;


/**
 * RectilinearGrid
 * Class used to visualize physical properties on a cartesian grid with Paraview 
 */

public class RectilinearGrid {

  int Nx, Ny, Nz ;
  Hashtable  physicalProperties ;
  
  
  
  
  public RectilinearGrid(int Nx, int Ny){
    this.Nx=Nx;
    this.Ny=Ny;
    this.Nz=1;
    this.physicalProperties= new Hashtable();
  }
  
  
  
  
  public void addField(double[][] property, String name){
    assert property.length == Nx && property[0].length==Ny :
      "This property doesn't have the same dimensions as the RectilinearGrid.";
    physicalProperties.put(name,property);
  }
    
  
  
  
  public void writeVTKFile(String name) throws IOException{
    Enumeration enumProperties = physicalProperties.elements();
    
    try {
      File file = new File(name);
 
      // if file doesnt exists, then create it
      if (!file.exists()) {
        file.createNewFile();
      }
 
      FileWriter fw = new FileWriter(file.getAbsoluteFile());
      BufferedWriter bw = new BufferedWriter(fw);

      //Write the VTK RectilinearGrid headline
      bw.write("DATASET RECTILINEAR_GRID"); 
      bw.newLine();
      bw.write("DIMENSIONS "+ Nx +" "+ Ny +" "+ Nz); 
      bw.newLine();

      bw.write("X_COORDINATES "+     ); 
      bw.newLine();
      for(int i=0; i<Nx ; i++){
        bw.write(i);
      }
      bw.newLine();

      bw.write("Y_COORDINATES "+     );
      bw.newLine();
      for(int i=0; i<Ny ; i++){
        bw.write(i);
      }
      bw.newLine();

      bw.write("Z_COORDINATES "+     );
      bw.newLine();
      for(int i=0; i<Nz ; i++){
        bw.write(i);
      }
      bw.newLine();

      //Write the properties
      while(enumProperties.hasMoreElements()){
        array =   enumProperties.nextElement();
        
        bw.write("CELL_DATA " +   );
        bw.newLine();
        bw.write("SCALARS "+dataName+ "float "+ 1    );
        bw.newLine();
        LOOKUP_ TABLE tableName 
        bw.newLine();
        
        
        for(int iLine=0 ; iLine <array.length  ; iLine++){
          for(int iColumn=0 ; iColumn < array[iLine].length ; iColumn++){
            bw.write(array[iLine][iColumn] +" ");
          }
          bw.newLine();
        }
        bw.newLine();
      }
      
      bw.close();
 
      System.out.println("Paraview file writen.");
 
    } catch (IOException ex) {
      ex.printStackTrace();
    }

  }
  
  
}
