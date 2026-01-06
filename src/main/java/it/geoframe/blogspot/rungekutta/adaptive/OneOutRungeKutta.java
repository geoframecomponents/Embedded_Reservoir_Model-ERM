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
package it.geoframe.blogspot.rungekutta.adaptive;

/**
 * 
 * Extension for single-out reservoir.
 * 
 * @author Giuseppe Formetta, Daniele Andreis
 *
 */
public class OneOutRungeKutta extends AdaptiveRungeKutta4 {
	private double d;
	private double c;
	private double maxStorage;

	public OneOutRungeKutta(double c, double d, double maxStorage) {
		this.c = c;
		this.d = d;
		this.maxStorage = maxStorage;
	}

	public double[] computeFunction(double storageN, double in) {
		if (storageN < 0) {
			storageN = 0;
		}
		double out = computeOut(storageN, in);
		double fun = in - out;
		return new double[] { fun, out };
	}

	private double computeOut(double storageN, double in) {
		double out = c * Math.pow( storageN, d);
		out = out + Math.max(0, storageN - maxStorage + in - out);
		out = Math.min(storageN + in, out);
		return out;
	}

	@Override
	protected int getOutDimension() {
		return 2;
	}
}
