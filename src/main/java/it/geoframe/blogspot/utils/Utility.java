package it.geoframe.blogspot.utils;

public class Utility {

	/**
	 * Get the conversion factor, from mm/dt to m3/s.
	 * 
	 * @param a area
	 * @param dt time in minutes
	 * @return
	 */
	
	public static final double getCOnversionToM3SCoeff(double a, double dt ) {
		 return a * Math.pow(10, 3) / (dt * 60);
	}
	
}
