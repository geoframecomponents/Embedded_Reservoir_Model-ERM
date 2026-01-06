/*
 * GNU GPL v3 License
 *
 * Copyright 2023 Daniele Andreis
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package it.geoframe.blogspot.rungekutta;

/**
 * 
 * @author Giuseppe Formetta, Daniele Andreis
 *
 */
public class Utils {

	/**
	 * Evaluate the Runge-Kutta 4 mean;
	 * 
	 * @param k1
	 * @param k2
	 * @param k3
	 * @param k4
	 * @param index
	 * @return
	 */
	public static final double getRKMean(double[] k1, double[] k2, double[] k3, double[] k4, int index) {
		return (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / 6;
	}

	/**
	 * Evaluate the Runge-Kutta 4 mean;
	 * 
	 * @param k1
	 * @param k2
	 * @param k3
	 * @param k4
	 * @param index
	 * @param den   denominator if is not 6.
	 * @return
	 */
	public static final double getRKMean(double[] k1, double[] k2, double[] k3, double[] k4, int index, double den) {
		return (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / den;
	}
}
