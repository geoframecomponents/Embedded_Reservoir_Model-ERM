/*
 * GNU GPL v3 License
 *
 * Copyright 2021 Niccolo` Tubini, Giuseppe Formetta
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
package it.geoframe.blogspot.simplebucket;

import static org.hortonmachine.gears.libs.modules.HMConstants.isNovalue;

import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import it.geoframe.blogspot.rungekutta.adaptive.AdaptiveRungeKutta4;
import it.geoframe.blogspot.rungekutta.adaptive.OneOutRungeKutta;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;
import oms3.annotations.Unit;


/**
 * The Class WaterBudget solves the water budget equation for the runoff layer.
 * The input s the recharge from the root zone and the output is the discharge,
 * modeled with a non linear reservoir model.
 * 
 * @author Marialaura Bancheri, Riccardo Busti
 */
public class WaterBudget {

	@Description("Input recharge Hashmap")
	@In
	public HashMap<Integer, double[]> inHMRechargeValues;

	@Description("Input CI Hashmap")
	@In
	public HashMap<Integer, double[]> initialConditionS_i;

	@Description("Time Step simulation")
	@Unit("minutes")
	@In
	public double tTimestep;

	@Description("Coefficient of the non-linear Reservoir model ")
	@In
	public double c;

	@Description("Exponent of non-linear reservoir")
	@In
	public double d;

	@Description("The area of the HRUs in km2")
	@In
	public double A;

	@Description("s_RunoffMax")
	@In
	public double s_RunoffMax;

	@Description("RK iterations")
	@In
	public static double RKiter = 10;

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

	 AdaptiveRungeKutta4 rk = null;
 double m3s = 0;

	
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
				rk = new OneOutRungeKutta(c, d, s_RunoffMax);
				m3s = A * Math.pow(10, 3) / (tTimestep * 60);

				if (initialConditionS_i != null) {
					CI = initialConditionS_i.get(ID)[0];
					if (isNovalue(CI))
						CI = 0.5 * s_RunoffMax;
				} else {
					CI = 0.5 * s_RunoffMax;
				}
				ciMap.put(ID, new double[] { CI });
			}
			double CI = ciMap.get(ID)[0];

			WaterBudgetStepResult r = calculateWaterBudget(rk,m3s,   recharge, CI);

			// save results
			storeResult_series(ID, r.waterStorage, r.runoff_mm, r.runoff, r.error);

			// update storage
			ciMap.put(ID, new double[] { r.waterStorage });
		}
		step++;
	}

	public static WaterBudgetStepResult calculateWaterBudget(AdaptiveRungeKutta4 rk, double m3s2,double recharge, double CI) {
		// solve S at t^n+1
		double[] output = rk.run(CI, recharge, RKiter);;
		double waterStorage = output[0];
		if (waterStorage < 0)
			waterStorage = 0;
		double error = output[2];

		// update variables at t^n+1
		double runoff_mm = output[1];
		double runoff = runoff_mm * m3s2;
		WaterBudgetStepResult r = new WaterBudgetStepResult(
				waterStorage,
				runoff,
				runoff_mm,
				error
		);
		return r;
	}





	private void storeResult_series(int ID, double S, double r_mm, double r, double err) {

		outHMStorage.put(ID, new double[] { S });
		outHMDischarge.put(ID, new double[] { r });
		outHMDischarge_mm.put(ID, new double[] { r_mm });
		outHMError.put(ID, new double[] { err });

	}
	
	public record WaterBudgetStepResult(
			double waterStorage,
			double runoff,
			double runoff_mm,
			double error
	) {
	}
}
