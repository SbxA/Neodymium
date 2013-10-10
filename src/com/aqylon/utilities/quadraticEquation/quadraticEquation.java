package com.aqylon.utilities.quadraticEquation;

public class quadraticEquation {

	public double r1, r2, solution;
	private double  a, b, c, delta;

	/**
	 * @author Vincent Lamblot
	 * Solves a quadratic equation for physics use only (positive roots are chosen) 
	 * a*x^2+b*x+c (x is the unknown)
	 * r1 & r2 are the roots
	 * solution is the chosen root (the positive one)
	 * @param a
	 * @param b
	 * @param c
	 */
	public quadraticEquation(double a, double b, double c){

		this.a = a;
		this.b = b;
		this.c = c;
		solve();
	}

	private void solve(){
		delta = b*b-4*a*c;
		if(delta>0.0){
			r1 = (-b+Math.sqrt(delta))/(2*a);
			r2 = (-b-Math.sqrt(delta))/(2*a);

			if(r1>0.0 && r2>0.0){
				throw new RuntimeException("Both roots are positive. Do not know which one to choose.");
			}
			else if(r1<0.0 && r2<0.0){
				throw new RuntimeException("Both roots are negative.");
			}
			else if(r1>0.0 && r2<0.0){
				solution = r1;
			}
			else if(r1<0.0 && r2>0.0){
				solution = r2;
			}
		}
		else if(delta==0){	
			r1 = -b/(2*a);
			r2 = r1;
		}

		else if(delta<0.0){
			throw new RuntimeException("Discriminant is negative. Do not know how to solve this case !!");
		}
		else{
			throw new RuntimeException("Unexpected error.");
		}
	}
}