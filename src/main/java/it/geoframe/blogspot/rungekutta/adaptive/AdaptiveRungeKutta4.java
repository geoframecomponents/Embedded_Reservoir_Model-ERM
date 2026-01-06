

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

import it.geoframe.blogspot.rungekutta.Utils;

/**
 * 
 * 
 * Runge Kutta with adaptive step
 * 
 * 
 * from:
 * https://www.glowscript.org/#/user/wlane/folder/Runge-Kutta/program/Runge-Kutta-Adaptive-Step
 * 
 * 
 * @author Giuseppe Formetta, Daniele Andreis
 *
 */
public abstract class AdaptiveRungeKutta4 {
	// the increments of the step.
	private static final double[] oneStepCoefficent = new double[] { 0.5, 0.5, 1 };
	private static final double[] halfStepCoefficent = new double[] { 0.25, 0.25, 0.5 };
	private static final double[] doubleStepCoefficent = new double[] { 1.0, 1.0, 2.0 };

	// control value
	static final double DT_MIN = 0.005;
	double dtMax = 0.1;
	// relative value!!!
	double dSMax = 0.1;
	// relative value!!!
	double dSMin = 0.01;
	double dSToll = 0.01;
	double[] output;

	protected AdaptiveRungeKutta4() {
	}

	// RK4
	public double[] run(double storageStart, double in, double rkiter) {
		double dt = Math.min(1.0 / rkiter,1.0);
		double dtMin =Math.max(DT_MIN,dt);
		double t = 0;
		output = new double[getOutDimension() + 1];
		output[0] = storageStart;
		dSToll = 0.0001 * storageStart;
		while (t < 1.0) {
			double[] k1 = computeFunction(output[0], in);
			double[] oneStepValue = this.computeValue(output[0], dt, k1, oneStepCoefficent, in);
			double oneStepSlope = oneStepValue[0];

			if (Math.abs(dt * oneStepSlope) < dSToll) {
				updateOutput(dt, in, oneStepValue);
				t = t + dt;
				dt = checkDt(t, dtMin);
				continue;
			}
			double[] halfStepValue = this.computeValue(output[0], dt, k1, halfStepCoefficent, in);
			double halfStepSlope = halfStepValue[0];
			if (Math.abs(dt * oneStepSlope) > dSToll
					&& Math.abs(oneStepSlope - halfStepSlope) / Math.abs(oneStepSlope) > dSMax) {
				updateOutput(dt / 2, in, halfStepValue);
				t = t + dt / 2;
				dt = checkDt(t, dt / 2);
				continue;
			}
			double[] doubleStepValue = this.computeValue(output[0], dt, k1, doubleStepCoefficent, in);
			double sN1DoubleStep = doubleStepValue[0];

			if (dt * 2 < dtMax && t + 2.0 * dt <= 1.0 && Math.abs(dt * oneStepSlope) > dSToll
					&& Math.abs(oneStepSlope - sN1DoubleStep) / Math.abs(oneStepSlope) < dSMin) {
				dt = dt * 2;
				updateOutput(dt, in, doubleStepValue);

			} else {
				updateOutput(dt, in, oneStepValue);
			}
			t = t + dt;
			dt = checkDt(t, dt);

		}
		return output;

	}

	private double checkDt(double t, double dt) {
		if (t + dt > 1.0) {
			return 1.0 - t;
		}
       if(dt<DT_MIN) {
    	   return DT_MIN;
       }
		return dt;
	}

	private void updateOutput(double dt, double in, double[] slopesValue) {
		double storageStart = output[0];
		double totalOutput = 0;
		output[0] = output[0] + dt * slopesValue[0];
		int i = 0;
		for (i = 1; i < output.length - 1; i++) {
			double tmpValue = dt * slopesValue[i];
			output[i] = output[i] + tmpValue;
			totalOutput = +tmpValue;
		}
		output[i] = output[i]+Math.abs(storageStart - output[0] + dt * in - totalOutput);
	}

	private double[] computeValue(double storagePreviousStep, double dt, double[] k1, double[] params, double in) {
		double[] k2 = computeFunction(storagePreviousStep + params[0] * dt * k1[0], in);
		double[] k3 = computeFunction(storagePreviousStep + params[1] * dt * k2[0], in);
		double[] k4 = computeFunction(storagePreviousStep + params[2] * dt * k3[0], in);
		double[] result = new double[k1.length];
		for (int i = 0; i < k1.length; i = i + 1) {
			result[i] = Utils.getRKMean(k1, k2, k3, k4, i);
		}

		return result;
	}

	protected abstract double[] computeFunction(double storage, double in);

	protected abstract int getOutDimension();
}
