/*
 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
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

package org.geoframe.erm.groundWater;

import static org.geoframe.erm.utility.Utils.getRKMean;
import static org.hortonmachine.gears.libs.modules.HMConstants.isNovalue;

import java.util.HashMap;
import java.util.Set;
import java.util.Map.Entry;

import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;

/**
 * The Class WaterBudget solves the water budget equation for the groudwater
 * layer. The input s the recharge from the root zone and the output is the
 * discharge, modeled with a non linear reservoir model.
 * 
 * @author Marialaura Bancheri, Riccardo Busti
 */
public class WaterBudgetGround {

	@Description("Input recharge Hashmap")
	@In
	public HashMap<Integer, double[]> inHMRechargeValues;

	@Description("Input CI Hashmap")
	@In
	public HashMap<Integer, double[]> initialConditionS_i;

	@Description("Time Step simulation")
	@In
	public double tTimestep;

	@Description("Coefficient of the non-linear Reservoir model ")
	@In
	public double e;

	@Description("Exponent of non-linear reservoir")
	@In
	public double f;

	@Description("The area of the HRUs in km2")
	@In
	public double A;

	@Description("s_GroundWaterMax")
	@In
	public double s_GroundWaterMax;

	@Description("RK iterations")
	@In
	public double RKiter = 100;

	// @Description("ODE solver model: dp853, Eulero ")
	// @In
	// public String solver_model;

	@Description("The output HashMap with the Water Storage")
	@Out
	public HashMap<Integer, double[]> outHMStorage = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the discharge")
	@Out
	public HashMap<Integer, double[]> outHMDischarge = new HashMap<Integer, double[]>();

	@Description("The output HashMap with the discharge in mm")
	@Out
	public HashMap<Integer, double[]> outHMDischarge_mm = new HashMap<Integer, double[]>();

	@Description("The output HashMap with error")
	@Out
	public HashMap<Integer, double[]> outHMError = new HashMap<Integer, double[]>();

	int step;
	private HashMap<Integer, double[]>ciMap= new HashMap<Integer, double[]>();

	/**
	 * Process: reading of the data, computation of the storage and outflows
	 *
	 * @throws Exception the exception
	 */
	@Execute
	public void process() throws Exception {

		// reading the ID of all the stations
		Set<Entry<Integer, double[]>> entrySet = inHMRechargeValues.entrySet();

		// iterate over the station
		for (Entry<Integer, double[]> entry : entrySet) {
			Integer ID = entry.getKey();

			/** Input data reading */
			double recharge = inHMRechargeValues.get(ID)[0];
			if (isNovalue(recharge))
				recharge = 0;

			if (step == 0) {
				double CI;
				if (initialConditionS_i != null) {
					CI = initialConditionS_i.get(ID)[0];
					if (isNovalue(CI))
						CI = 0.01 * s_GroundWaterMax;
				} else {
					CI = 0.01 * s_GroundWaterMax;
				}
				ciMap.put(ID, new double[] { CI });
			}
			double CI = ciMap.get(ID)[0];

			WaterBudgetGroundStepResult r = calculateWaterBudgetGround(recharge, CI, e, f, s_GroundWaterMax, A, tTimestep);

			// save results
			storeResult_series(ID, r.waterStorage, r.discharge_mm, r.discharge, r.error);

			// update storage
			ciMap.put(ID, new double[] { r.waterStorage });
		}
		step++;
	}

	public static WaterBudgetGroundStepResult calculateWaterBudgetGround(double recharge, double CI,  double e, double f, double s_GroundWaterMax, double A, double tTimestepInMinutes) {
		double m3s = A * Math.pow(10, 3) / (tTimestepInMinutes * 60);

		// solve S at t^n+1
		double[] out = RK4(CI, recharge, e, f, s_GroundWaterMax);
		double waterStorage = out[0];
		if (waterStorage < 0)
			waterStorage = 0;
		double error = out[1];

		// update variables at t^n+1
		double deep_mm = out[2];
		double deep = deep_mm * m3s;
		
		WaterBudgetGroundStepResult r = new WaterBudgetGroundStepResult(
				waterStorage,
				deep,
				deep_mm,
				error);
		return r;
	}

	// compute dS/dt
	public static double[] computeFunction(double Sn, double recharge, double e, double f, double s_GroundWaterMax) {
		if (Sn < 0) {
			Sn = 0;
		}
		double deep = computeDeep(Sn, recharge, e, f, s_GroundWaterMax);
		double fun = recharge - deep;
		return new double[] { fun, deep };
	}

	// compute deep discharge
	public static double computeDeep(double Sn, double recharge, double e, double f, double s_GroundWaterMax) {
		double out = e * Math.pow(Math.min(1, Sn / s_GroundWaterMax), f);
		out = out + Math.max(0, Sn - s_GroundWaterMax + recharge - out);
		return Math.min(Sn + recharge, out);
	}

	// RK4
	public static double[] RK4(double Sn, double recharge, double e, double f, double s_GroundWaterMax) {

		double balance = 0;
		double[] k1 = computeFunction(Sn, recharge, e, f, s_GroundWaterMax);
		double[] k2 = computeFunction(Sn + 0.5 * k1[0], recharge, e, f, s_GroundWaterMax);
		double[] k3 = computeFunction(Sn + 0.5 * k2[0], recharge, e, f, s_GroundWaterMax);
		double[] k4 = computeFunction(Sn + k3[0], recharge, e, f, s_GroundWaterMax);
		double Sn1 = Sn + getRKMean(k1, k2, k3, k4, 0);
		double deltaDeep = getRKMean(k1, k2, k3, k4, 1);
		balance = balance + Math.abs(Sn - Sn1 + recharge - deltaDeep);
		return new double[] { Sn1, balance, deltaDeep };

	}

	private void storeResult_series(int ID, double S, double d_mm, double d, double err) {
		outHMStorage.put(ID, new double[] { S });
		outHMDischarge.put(ID, new double[] { d });
		outHMDischarge_mm.put(ID, new double[] { d_mm });
		outHMError.put(ID, new double[] { err });
	}
	
	public record WaterBudgetGroundStepResult(
			double waterStorage,
			double discharge,
			double discharge_mm,
			double error) {
	}

}