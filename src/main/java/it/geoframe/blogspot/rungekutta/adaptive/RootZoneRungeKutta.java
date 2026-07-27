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

import static org.hortonmachine.gears.libs.modules.HMConstants.isNovalue;

/**
 * Extension for root zone.
 * 
 * 
 * @author Giuseppe Formetta, Daniele Andreis
 *
 */
public class RootZoneRungeKutta extends AdaptiveRungeKutta4 {
	double sRootZoneMax;
	double pBSoil;
	double rain;
	double coeff;
	double exponent;

	public RootZoneRungeKutta(double coeff, double exponent, double sRootZoneMax, double pBSoil) {
		this.sRootZoneMax = sRootZoneMax;
		this.pBSoil = pBSoil;
		this.coeff = coeff;
		this.exponent = exponent;
	}

	@Override
	protected int getOutDimension() {
		return 6;
	}

	public double[] run(double storageStart, double in, double out, double rkiter) {
		this.rain = in;
		double[] result = super.run(storageStart, out, rkiter);
		double S_new = result[0];
		double AET = result[3];
		double recharge = result[2];
		double quick = result[5];
		result[getOutDimension()] = Math.abs(S_new - storageStart - this.rain + AET + recharge + quick);
		return result;
	}

	// compute dS/dt
	public double[] computeFunction(double sN, double etpnet) {
		if (sN < 0) {
			sN = 0;
		}
		double alpha = alpha(sN, rain);
		double[] o = actualInputs(alpha);
		double actualInputs = o[0];
		double quick = o[1];
		double aet = computeAET(sN, actualInputs, etpnet);
		double recharge = computeR(sN, actualInputs, aet);
		double fun = actualInputs - aet - recharge;
		return new double[] { fun, actualInputs, recharge, aet, alpha, quick };
	}

	// compute alpha according to Hymod
	private double alpha(double sN, double pVal) {
		double pCmax = sRootZoneMax * (pBSoil + 1);
		double coeff1 = 1.0 - ((pBSoil + 1.0) * (sN) / pCmax);
		double exp = 1.0 / (pBSoil + 1.0);
		double ctPrev = pCmax * (1.0 - Math.pow(coeff1, exp));
		double uT1 = Math.max((pVal - pCmax + ctPrev), 0.0);
		double dummy = Math.min(((ctPrev + pVal - uT1) / pCmax), 1.0);
		double coeff2 = (1.0 - dummy);
		double exp2 = (pBSoil + 1.0);
		double xn = (pCmax / (pBSoil + 1.0)) * (1.0 - (Math.pow(coeff2, exp2)));
		double uT2 = Math.max(pVal - uT1 - (xn - sN), 0);
		double alpha = (uT1 + uT2) / pVal;
		if (isNovalue(alpha) || alpha > 1)
			alpha = 1;
		return alpha;
	}

	// compute actual inputspublic
	private double[] actualInputs( double alfa) {
		return new double[] { (1 - alfa) * rain, alfa * rain };
	}

	// compute groundwater recharge
	private double computeR(double sn, double in, double et) {
		double out = coeff * Math.pow(sn, exponent);
		out = Math.min(sn + in - et, out + Math.max(0, sn - sRootZoneMax + in - et - out));
		return out;

	}

	// compute AET
	private double computeAET(double sN, double in, double etpnet) {
		return Math.min(sN + in, etpnet * Math.min(1, 1.33 * Math.min(1, sN / sRootZoneMax)));
	}

}
